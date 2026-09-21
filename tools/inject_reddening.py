"""Closure test of the law-mode fit: inject a known G23 reddening into the
(dereddened) template calibrators and refit with A_V, R_V and T_eff free, as
the science fields are fitted.  Measures the R_V-T_eff-A_V coupling directly.

  python tools/inject_reddening.py dwarfs 1.0 3.05 free [plx10] [max_stars]
  python tools/inject_reddening.py giants 1.0 3.05 locked

variants: free   - T_eff free, log g / [Fe/H] priors as in the field (DESI / APOGEE)
          locked - T_eff locked to the calibration value (colour / ASPCAP), A_V, R_V free
          nocorr - free, template corrections off
plx10: inflate the parallax errors to S/N 10 (the field stars at 1-2 kpc), so the
       radius prior is as weak as in the field.

Output: inject_<sample>_av<A>_rv<R>_<variant>[_plx10].csv in the calibrator workspace.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from build_template_corrections import stars_table, workspace  # noqa: E402

from dustline import calib, extinction, fit as fit_mod  # noqa: E402

BANDS = ["2MASS_H", "2MASS_J", "2MASS_Ks", "PS1_g", "PS1_i", "PS1_r", "PS1_y", "PS1_z",
         "WISE_W1", "WISE_W2"]


def load_sample(sample: str):
    ws = workspace(sample)
    stars, _ = stars_table(ws)
    xp = dict(np.load(ws.path("xp_sampled.npz")))
    av0 = pd.read_csv(ws.path("xp_stars_av0.csv"))
    truth = av0[["source_id", "av_map", "imodel", "teff_spec", "logg_spec", "feh_spec",
                 "teff_best", "logg_best", "mh_best"]].rename(
        columns={"teff_spec": "teff_true", "logg_spec": "logg_true", "feh_spec": "feh_true"})
    stars = stars.merge(truth, on="source_id")
    keep = np.isin(xp["source_id"], stars.source_id.values)
    xp = {k: (v[keep] if k != "wave_nm" else v) for k, v in xp.items()}
    ids = list(xp["source_id"])
    stars = stars.set_index("source_id").loc[ids].reset_index()
    xp, stars = calib.deredden(xp, stars, BANDS, stars.av_map.values)
    return ws, stars, xp


def inject(ws, stars, xp, av_inj: float, rv_inj: float):
    """Redden the XP spectra with G23(rv_inj) at av_inj and the photometry with the
    band extinctions of each star's calibration best-fit model (uncorrected grid,
    the same products the science fit uses for the law)."""
    grid0 = fit_mod.build_grid(ws, BANDS, mh=None, corrections=None)
    ir = int(np.argmin(np.abs(grid0["rv"] - rv_inj)))
    ia = int(np.argmin(np.abs(grid0["av"] - av_inj)))
    assert abs(grid0["rv"][ir] - rv_inj) < 1e-6 and abs(grid0["av"][ia] - av_inj) < 1e-6, \
        "inject on grid values so the recovery is not quantised"
    ext = extinction.curve(np.asarray(xp["wave_nm"]) * 10.0, rv_inj)
    scale = 10.0 ** (-0.4 * av_inj * ext)
    out = dict(xp)
    out["flux"] = xp["flux"] * scale[None, :]
    out["flux_err"] = xp["flux_err"] * scale[None, :]
    st = stars.copy()
    im = st.imodel.values.astype(int)
    for j, b in enumerate(BANDS):
        st[f"mag_{b}"] = st[f"mag_{b}"] + grid0["A_band"][ir, ia, im, j]
    st["av_inj"], st["rv_inj"] = av_inj, rv_inj
    return st, out


def main():
    sample, av_inj, rv_inj = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    variant = sys.argv[4] if len(sys.argv) > 4 else "free"
    plx10 = "plx10" in sys.argv[5:]
    nums = [a for a in sys.argv[5:] if a.isdigit()]
    max_stars = int(nums[0]) if nums else 0
    ws, stars, xp = load_sample(sample)
    trange = [a for a in sys.argv[5:] if a.startswith("teff=")]     # teff=3000,4500
    if trange:
        lo, hi = (float(x) for x in trange[0][5:].split(","))
        keep = (stars.teff_true >= lo) & (stars.teff_true < hi)
        stars = stars[keep].reset_index(drop=True)
        xp = {k: (v[keep.values] if k != "wave_nm" else v) for k, v in xp.items()}
    if "synth" in sys.argv[5:]:
        # machinery check: replace each XP spectrum by its calibration model (corrected
        # template x C_best) plus Gaussian noise at the reported errors
        av0 = pd.read_csv(ws.path("xp_stars_av0.csv")).set_index("source_id")
        gridc = fit_mod.build_grid(ws, BANDS, mh=None, corrections="default")
        C = av0.C_best.reindex(stars.source_id).values
        m = gridc["m_xp"][stars.imodel.values.astype(int)] * C[:, None]
        rng = np.random.default_rng(1)
        xp = dict(xp, flux=(m + rng.normal(size=m.shape) * np.nan_to_num(xp["flux_err"], nan=0.0)
                            ).astype(xp["flux"].dtype))
        variant += "_synth"
    if "flip" in sys.argv[5:]:
        # mirror each star's residual about its calibration model: obs -> model^2 / obs
        av0 = pd.read_csv(ws.path("xp_stars_av0.csv")).set_index("source_id")
        gridc = fit_mod.build_grid(ws, BANDS, mh=None, corrections="default")
        C = av0.C_best.reindex(stars.source_id).values
        m = gridc["m_xp"][stars.imodel.values.astype(int)] * C[:, None]
        f = xp["flux"].astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            flipped = np.where(f > 0, m * m / f, f)
        xp = dict(xp, flux=flipped.astype(xp["flux"].dtype))
        variant += "_flip"
    stars, xp = inject(ws, stars, xp, av_inj, rv_inj)
    if variant.startswith("locked"):
        stars["teff_spec"] = stars.teff_true
        stars["teff_spec_err"] = calib.TEFF_PRIOR_SIGMA
    else:
        stars["teff_spec"] = np.nan          # law mode: T_eff from the data + radius prior
    if sample == "giants":
        stars["logg_spec_err"] = 0.1
    if plx10:
        stars["parallax_error"] = np.maximum(stars.parallax_error, stars.parallax / 10.0)
    stars["teff_desi"] = stars.teff_true
    drop = [a[5:] for a in sys.argv[5:] if a.startswith("drop=")]        # drop=PS1,WISE / drop=all
    for pref in (drop[0].split(",") if drop else []):
        for b in BANDS:
            if pref == "all" or b.startswith(pref):
                stars[f"mag_{b}"] = np.nan
    corr = None if variant == "nocorr" else "default"
    fit = fit_mod.fit_stars(ws, stars, xp, BANDS, offsets=None, corrections=corr, max_stars=max_stars)
    fit = fit.merge(stars[["source_id", "av_map", "av_inj", "rv_inj", "teff_true", "logg_true",
                           "feh_true"]], on="source_id", how="left")
    name = (f"inject_{sample}_av{av_inj:.2f}_rv{rv_inj:.2f}_{variant}" + ("_plx10" if plx10 else "")
            + (f"_T{trange[0][5:].replace(',', '-')}" if trange else "") + (f"_n{max_stars}" if max_stars else "")
            + (f"_drop{drop[0].replace(',', '')}" if drop else ""))
    fit.round(5).to_csv(ws.path(name + ".csv"), index=False)
    print(f"wrote {ws.path(name + '.csv')}: {len(fit)} stars")
    summarize(fit)


def summarize(fit: pd.DataFrame):
    d = fit.copy()
    d["dT"] = d.teff - d.teff_true
    d["dA"] = d.av - d.av_inj
    d["dR"] = d.rv - d.rv_inj
    print(f"all N {len(d)}: R_V {d.rv.median():.2f} (MAD {1.4826 * (d.rv - d.rv.median()).abs().median():.2f}),"
          f" dA_V {d.dA.median():+.3f}, dT_eff {d.dT.median():+.0f} K, mean R_V {d.rv.mean():.2f}")
    edges = np.array([3000, 4000, 4500, 5000, 5500, 6000, 6500, 7500, 12000])
    print("by true T_eff:")
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = d[(d.teff_true >= lo) & (d.teff_true < hi)]
        if len(s) < 5:
            continue
        print(f"  {lo:5d}-{hi:5d} N {len(s):4d}  R_V {s.rv.median():.2f}  dA_V {s.dA.median():+.3f}"
              f"  dT_eff {s.dT.median():+5.0f} K  (16-84: {s.dT.quantile(.16):+5.0f} {s.dT.quantile(.84):+5.0f})")
    print("by fitted A_V bin (the field diagnostic):")
    for q in np.linspace(0, 1, 5)[:-1]:
        lo, hi = d.av.quantile(q), d.av.quantile(q + 0.25)
        s = d[(d.av >= lo) & (d.av <= hi)]
        print(f"  A_V {lo:.2f}-{hi:.2f} N {len(s):4d}  R_V {s.rv.median():.2f}  dT_eff {s.dT.median():+5.0f} K")
    print("by fitted - true T_eff:")
    for lo, hi in ((-2000, -150), (-150, -50), (-50, 50), (50, 150), (150, 2000)):
        s = d[(d.dT >= lo) & (d.dT < hi)]
        if len(s) < 5:
            continue
        print(f"  dT {lo:+5d}..{hi:+5d} N {len(s):4d}  R_V {s.rv.median():.2f}  dA_V {s.dA.median():+.3f}")


def stage_closure(out: Path, min_per_node: int = 30) -> None:
    """The R_V closure table from the plx10 free runs (A_V 1, R_V 3.05 injected):
    per grid node and log g class, k = A_V * median(1/R_V_fit - 1/R_V_inj), the
    additive residual of the corrected template in the R_V direction (it scales
    as 1/A_V: 0.030 / 0.029 / 0.036 at A_V 0.5 / 1 / 2 for the cool dwarfs).
    Nodes with fewer stars take the class's pooled cool (< 5000 K dwarfs,
    < 4500 K giants) or warm value.  Written into the corrections table
    (rv_closure_k, rv_closure_n) for ensemble.measure_law."""
    corr = dict(np.load(out, allow_pickle=False))
    edges = corr["teff_edges"]
    nT, nG = len(edges) - 1, corr["ratio"].shape[2]
    k_tab = np.zeros((nT, nG))
    n_tab = np.zeros((nT, nG), int)
    pooled = {}
    for gcls, sample in enumerate(("dwarfs", "giants")):
        ws = workspace(sample)
        p = ws.path("inject_%s_av1.00_rv3.05_free_plx10.csv" % sample)
        if not p.exists():
            print(f"closure: no {p.name} for {sample} - skipped")
            continue
        d = pd.read_csv(p)
        d = d[np.isfinite(d.rv) & (d.rv > 0) & (d.teff_true < 6000)]
        dk = d.av_inj * (1.0 / d.rv - 1.0 / d.rv_inj)
        cool = d.teff_true < (5000.0 if sample == "dwarfs" else 4500.0)
        pooled[gcls] = (float(dk[cool].median()), float(dk[~cool].median()))
        iT = np.clip(np.digitize(d.teff_true.values, edges) - 1, 0, nT - 1)
        for t in range(nT):
            m = iT == t
            n_tab[t, gcls] = int(m.sum())
            if m.sum() >= min_per_node:
                k_tab[t, gcls] = float(dk[m].median())
            else:
                mid = 0.5 * (edges[t] + edges[t + 1])
                k_tab[t, gcls] = pooled[gcls][0 if mid < (5000.0 if sample == "dwarfs" else 4500.0) else 1]
        print(f"closure {sample}: N {len(d)}, pooled k cool {pooled[gcls][0]:+.4f} warm {pooled[gcls][1]:+.4f}")
    corr["rv_closure_k"] = k_tab
    corr["rv_closure_n"] = n_tab
    np.savez_compressed(out, **corr)
    nodes = 0.5 * (edges[:-1] + edges[1:])
    for gcls, name in enumerate(("dwarfs", "giants")):
        print(f"  {name}: " + " ".join(f"{nodes[t]:.0f}:{k_tab[t, gcls]:+.3f}({n_tab[t, gcls]})"
                                       for t in range(nT) if n_tab[t, gcls] >= min_per_node))
    print(f"wrote {out}")


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "closure":
        from dustline import assets
        stage_closure(Path(sys.argv[2]) if len(sys.argv) > 2
                      else assets.cache_dir() / "template_corrections.npz")
    elif len(sys.argv) == 2 and sys.argv[1].endswith(".csv"):
        summarize(pd.read_csv(sys.argv[1]))
    else:
        main()

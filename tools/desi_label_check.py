"""Validate the DESI DR1 MWS log g and [Fe/H] labels against dustline's own NewEra fits.

The fits normally take DESI log g and [Fe/H] as Gaussian priors (catalogs.desi.shape). Here they
are removed, so log g and [M/H] come only from the Gaia XP spectrum, the photometry and the
parallax radius prior (MIST radius at the model's T_eff and log g), and are compared with DESI.

Variants (per field, one fit pass at the field's converged zero points, R_V 3.1, column mode):
  tlock  T_eff locked to the scale-corrected DESI label (calib.desi_teff_offset), log g / [M/H] free
  free   T_eff free too (a 3000 K prior keeps the spectroscopic [M/H] grid)
[M/H] runs over the whole NewEra grid (-2..+0.5); the clamp that holds unlabelled stars to
-0.5..+0.5 (fit.MH_FREE_RANGE) is switched off. NewEra [M/H] is scaled-solar: DESI is also
compared as [M/H]_DESI = [Fe/H] + log10(0.638 x 10^[a/Fe] + 0.362) (Salaris+93).

  python tools/desi_label_check.py run --fields F.csv --out DIR [--jobs 40]
  python tools/desi_label_check.py summary --out DIR
"""
import argparse
import json
import os
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
RADIUS = 30.0


def one(args):
    name, ra, dec, variant, threads, out = args
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(threads)
    from dustline import bands as bands_mod
    from dustline import calib
    from dustline import fit as F
    from dustline.api import Sightline
    from dustline.catalogs import gaia
    t0 = time.time()
    try:
        F.MH_FREE_RANGE = (-9.0, 9.0)                      # no [M/H] clamp for unlabelled stars
        sl = Sightline(ra, dec, RADIUS)
        if not sl.ws.has("xp_sampled.npz"):
            sl.ws.seed_from_sibling()
        g = gaia.cone(sl.ws); xp = gaia.calibrate_xp(sl.ws)
        stars, bands = bands_mod.assemble(sl.ws, g, sl.plan)
        off = json.loads(sl.ws.path("phot_offsets.json").read_text()) if sl.ws.has("phot_offsets.json") else {}
        offsets = {b: float(v) for b, v in off.items() if isinstance(v, (int, float))}
        s = stars.copy()
        lab = s[["source_id", "teff_spec", "logg_spec", "feh_spec", "spec_snr"]].rename(
            columns=dict(teff_spec="teff_desi", logg_spec="logg_desi", feh_spec="feh_desi", spec_snr="snr_desi"))
        if "alphafe_spec" in s:
            lab["alphafe_desi"] = s["alphafe_spec"].values
        snr = s["spec_snr"].values.astype(float)
        t = s["teff_spec"].values.astype(float)
        t_scaled = t - calib.desi_teff_offset(t, snr)
        s["teff_spec"] = t_scaled
        s["teff_spec_err"] = np.where(snr >= calib.DESI_LOWSNR, calib.TEFF_PRIOR_SIGMA, calib.DESI_LOWSNR_SIGMA) \
            if variant == "tlock" else 3000.0
        s["logg_spec"] = np.nan; s["feh_spec"] = np.nan
        f = F.fit_stars(sl.ws, s, xp, bands, offsets=offsets, spectro_priors=True, rv_fixed=3.1, plx_inflate=1.0)
        f = f.merge(lab, on="source_id", how="left")
        f["teff_desi_scaled"] = f.teff_desi - calib.desi_teff_offset(f.teff_desi.values, f.snr_desi.values)
        f["field"] = name; f["variant"] = variant
        f.to_csv(os.path.join(out, f"{name}.{variant}.csv"), index=False)
        return dict(name=name, variant=variant, ok=True, n=int(len(f)), seconds=round(time.time() - t0))
    except Exception as e:  # noqa: BLE001
        return dict(name=name, variant=variant, ok=False, error=f"{type(e).__name__}: {e}", tb=traceback.format_exc()[-1500:])


def cmd_run(a):
    os.makedirs(a.out, exist_ok=True)
    fl = pd.read_csv(a.fields)
    todo = [(r.name, r.ra, r.dec, v, a.threads, a.out) for r in fl.itertuples() for v in ("tlock", "free")
            if not os.path.exists(os.path.join(a.out, f"{r.name}.{v}.csv"))]
    print(f"{len(todo)} fits", flush=True)
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for fut in as_completed([ex.submit(one, t) for t in todo]):
            r = fut.result()
            print(f"  {r['name']:8s} {r['variant']:6s} " + (f"ok {r['n']} stars {r['seconds']} s" if r["ok"] else "FAIL " + r["error"]), flush=True)


def rs(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return 1.4826 * np.median(np.abs(v - np.median(v))) if len(v) else np.nan


def cmd_summary(a):
    import glob
    d = pd.concat([pd.read_csv(p) for p in glob.glob(os.path.join(a.out, "*.csv")) if p.count(".") >= 2 and "summary" not in p],
                  ignore_index=True)
    plx_snr = d.parallax / d.parallax_error
    q = ((d.chi2_best * 3 / d.n_xp < 2.5) & (d.ruwe < 1.4) & (plx_snr > 5) & d.logg_desi.notna() & d.feh_desi.notna()
         & (d.snr_desi >= 10) & (d.n_phot >= 4))
    d = d[q].copy()
    d["mh_desi"] = d.feh_desi + np.log10(0.638 * 10 ** d.get("alphafe_desi", pd.Series(0.0, index=d.index)).fillna(0.0) + 0.362)
    d["dlogg"] = d.logg - d.logg_desi
    d["dfeh"] = d.mh - d.feh_desi
    d["dmh"] = d.mh - d.mh_desi
    d.round(4).to_csv(os.path.join(a.out, "label_check_summary.csv"), index=False)
    L = []
    for var, s in d.groupby("variant"):
        L.append(f"=== variant {var}: {len(s)} stars in {s.field.nunique()} fields (DESI S/N >= 10, plx S/N > 5)")
        L.append(f"  log g: fit - DESI median {np.median(s.dlogg):+.3f} (robust sd {rs(s.dlogg):.3f}); "
                 f"[M/H]_fit - [Fe/H]_DESI {np.median(s.dfeh):+.3f} ({rs(s.dfeh):.3f}); "
                 f"- [M/H]_DESI (alpha-corrected) {np.median(s.dmh):+.3f} ({rs(s.dmh):.3f})")
        L.append(f"  median fit errors: log g {s.logg_err.median():.3f}, [M/H] {s.mh_err.median():.3f}; T_eff fit - DESI(scaled) "
                 f"{np.median(s.teff - s.teff_desi_scaled):+.0f} K")
        for lab, col, edges in (("T_eff (scaled DESI)", "teff_desi_scaled", [4000, 4750, 5250, 5750, 6250, 7000]),
                                ("log g (DESI)", "logg_desi", [0, 3.5, 4.0, 4.4, 6]),
                                ("[Fe/H] (DESI)", "feh_desi", [-3, -1.5, -1.0, -0.5, -0.2, 0.1, 0.6]),
                                ("DESI S/N", "snr_desi", [10, 20, 40, 80, 1e5]),
                                ("Gaia G", "phot_g_mean_mag", [8, 14, 15.5, 16.5, 18])):
            L.append(f"  by {lab}:  bin: N  dlogg (sd)  d[Fe/H] (sd)  d[M/H]_alpha")
            b = np.digitize(s[col], edges) - 1
            for i in range(len(edges) - 1):
                m = b == i
                if m.sum() >= 15:
                    L.append(f"    {edges[i]:g}-{edges[i + 1]:g}: {m.sum():5d}  {np.median(s.dlogg[m]):+.3f} ({rs(s.dlogg[m]):.3f})  "
                             f"{np.median(s.dfeh[m]):+.3f} ({rs(s.dfeh[m]):.3f})  {np.median(s.dmh[m]):+.3f}")
        c = np.polyfit(s.feh_desi, s.mh, 1)
        L.append(f"  [M/H]_fit = {c[1]:+.3f} + {c[0]:.3f} x [Fe/H]_DESI (slope 1 = same scale)")
    txt = "\n".join(L)
    open(os.path.join(a.out, "label_check_summary.txt"), "w").write(txt + "\n")
    print(txt)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("--fields", required=True); r.add_argument("--out", required=True)
    r.add_argument("--jobs", type=int, default=40); r.add_argument("--threads", type=int, default=4)
    s = sub.add_parser("summary"); s.add_argument("--out", required=True)
    a = ap.parse_args()
    {"run": cmd_run, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    main()

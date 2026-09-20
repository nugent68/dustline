"""Build the empirical NewEra dwarf-template corrections (template_corrections.npz)
from unreddened DESI x Gaia-XP dwarfs.  See dustline.calib for the rationale.

  python tools/build_template_corrections.py pull    # calibrators, Gaia, XP (~1.5 h), PS1, WISE, GALEX
  python tools/build_template_corrections.py fit     # A_V = 0 fits with DESI priors
  python tools/build_template_corrections.py build   # aggregate -> template_corrections.npz
  python tools/build_template_corrections.py all

Everything lives in the calibration workspace ~/.cache/dustline/sightlines/ra+000..._calib.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from dustline import calib, fit as fit_mod
from dustline.cache import Workspace
from dustline.catalogs import datalab, gaia, galex, ps1, vizier, wise
from dustline.catalogs.xmatch import match_to_gaia

BANDS_PS1 = {"g": "PS1_g", "r": "PS1_r", "i": "PS1_i", "z": "PS1_z", "y": "PS1_y"}


def workspace() -> Workspace:
    return Workspace(0.0, 0.0, 0.0, dict(calib="template_corrections"))


def stage_pull(ws: Workspace, max_per_bin: int = 250) -> None:
    if not ws.has("calibrators.csv"):
        d = calib.select_calibrators(max_per_bin=max_per_bin)
        d.to_csv(ws.path("calibrators.csv"), index=False)
        print(f"{len(d)} calibrators")
    cal = pd.read_csv(ws.path("calibrators.csv"))
    if not ws.has("gaia.csv"):
        g = gaia.by_ids(cal.source_id.values)
        g.to_csv(ws.path("gaia.csv"), index=False)
    g = pd.read_csv(ws.path("gaia.csv"), low_memory=False)
    if not ws.has("xp_sampled.npz"):
        gaia.fetch_xp(ws, g, chunk=200)
    gaia.calibrate_xp(ws)
    if not ws.has("ps1_gaia.csv"):
        p = ps1.by_positions(g.ra.values, g.dec.values)
        res = pd.DataFrame(dict(source_id=g.source_id.values))
        for b, name in BANDS_PS1.items():
            good = ((p[f"{b}QfPerfect"].values > 0.85) & (p[f"{b}MeanPSFMagNpt"].values >= 2)
                    & (p[f"{b}MeanPSFMag"].values > ps1.SATURATION[b]))
            res[f"mag_{name}"] = np.where(good, p[f"{b}MeanPSFMag"].values, np.nan)
            res[f"magerr_{name}"] = np.where(good, np.sqrt(p[f"{b}MeanPSFMagErr"].values ** 2
                                                          + ps1.SYS_FLOOR ** 2), np.nan)
        res.round(4).to_csv(ws.path("ps1_gaia.csv"), index=False)
        print(f"PS1: {np.isfinite(res.mag_PS1_g).sum()} with g")
    if not ws.has("wise_gaia.csv"):
        # positions are scattered: chunked OR-box query, then the standard shaping per chunk centre
        d = datalab.positions_query("allwise.source", wise._COLS, g.ra.values, g.dec.values, 3.0)
        res = _shape_scattered(wise.shape, d, g, "ra", "dec")
        res.round(4).to_csv(ws.path("wise_gaia.csv"), index=False)
        print(f"AllWISE: {np.isfinite(res.mag_WISE_W1).sum()} with W1")
    if not ws.has("galex_gaia.csv"):
        d = vizier.positions_query("II/335/galex_ais", galex._COLS, g.ra.values, g.dec.values, 3.0)
        res = _shape_scattered(galex.shape, d, g, "RAJ2000", "DEJ2000")
        res.round(4).to_csv(ws.path("galex_gaia.csv"), index=False)
        print(f"GALEX: {np.isfinite(res.mag_GALEX_NUV).sum()} with NUV")


def _shape_scattered(shape, d: pd.DataFrame, g: pd.DataFrame, ra_col: str, dec_col: str) -> pd.DataFrame:
    """The catalog modules' shape() assume one field centre; run them per star."""
    out = []
    for _, row in g.iterrows():
        sep = np.hypot((d[ra_col] - row.ra) * np.cos(np.radians(row.dec)), d[dec_col] - row.dec) * 3600
        near = d[sep < 5.0]
        if len(near):
            out.append(shape(near.reset_index(drop=True), g[g.source_id == row.source_id], row.ra, row.dec))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(dict(source_id=[]))


def stars_table(ws: Workspace) -> tuple[pd.DataFrame, list[str]]:
    g = pd.read_csv(ws.path("gaia.csv"), low_memory=False)
    for b in ("J", "H", "Ks"):
        g = g.rename(columns={f"mag_{b}": f"mag_2MASS_{b}", f"magerr_{b}": f"magerr_2MASS_{b}"})
    for f in ("ps1_gaia.csv", "wise_gaia.csv", "galex_gaia.csv"):
        x = pd.read_csv(ws.path(f))
        x = x.drop(columns=[c for c in x.columns if c.endswith("_sep")])
        if f.startswith("ps1"):
            for b, name in BANDS_PS1.items():      # saturation cut (older caches)
                sat = x[f"mag_{name}"] <= ps1.SATURATION[b]
                x.loc[sat, [f"mag_{name}", f"magerr_{name}"]] = np.nan
        g = g.merge(x, on="source_id", how="left")
    cal = pd.read_csv(ws.path("calibrators.csv"))
    sp = pd.DataFrame(dict(source_id=cal.source_id, teff_spec=cal.teff,
                           teff_spec_err=np.full(len(cal), calib.TEFF_PRIOR_SIGMA),
                           logg_spec=cal.logg, logg_spec_err=np.sqrt(cal.logg_err ** 2 + 0.1 ** 2),
                           feh_spec=cal.feh, feh_spec_err=np.sqrt(cal.feh_err ** 2 + 0.1 ** 2),
                           spec_snr=cal.snr_med))
    g = g.merge(sp, on="source_id", how="left")
    bands = sorted({c[4:] for c in g.columns if c.startswith("mag_") and f"magerr_{c[4:]}" in g.columns})
    return g, bands


def stage_fit(ws: Workspace) -> None:
    stars, bands = stars_table(ws)
    xp = dict(np.load(ws.path("xp_sampled.npz")))
    # the calibrators are nearby but not dust-free (z ~ 80-190 pc at |b| > 40 is behind
    # most of the dust layer): remove each star's 3D-map A_V before the A_V = 0 fit
    av_map = calib.map_extinction(stars.ra.values, stars.dec.values, 1000.0 / stars.parallax.values)
    stars["av_map"] = av_map
    print(f"3D-map A_V of the calibrators: median {np.median(av_map):.3f}, "
          f"16-84 % {np.percentile(av_map, 16):.3f}-{np.percentile(av_map, 84):.3f}")
    xp, stars = calib.deredden(xp, stars, bands, av_map)
    # the T_eff coordinate: the empirical dwarf T_eff of the dereddened BP-RP, locked
    stars["teff_desi"] = stars["teff_spec"]
    stars["teff_spec"] = calib.teff_from_bprp(stars.bp_rp.values, av_map, stars.feh_spec.values)
    stars["teff_spec_err"] = calib.TEFF_PRIOR_SIGMA
    ok = np.isfinite(stars.teff_spec)
    print(f"colour T_eff for {int(ok.sum())} of {len(stars)} calibrators; "
          f"median colour-DESI = {np.nanmedian(stars.teff_spec - stars.teff_desi):+.0f} K")
    stars = stars[ok].reset_index(drop=True)
    keep = np.isin(xp["source_id"], stars.source_id.values)
    xp = {k: (v[keep] if k != "wave_nm" else v) for k, v in xp.items()}
    fit = fit_mod.fit_stars(ws, stars, xp, bands, offsets=None, av_fixed=0.0, rv_fixed=3.1,
                            corrections=None)
    fit = fit.merge(stars[["source_id", "av_map"]], on="source_id", how="left")
    fit.round(5).to_csv(ws.path("xp_stars_av0.csv"), index=False)
    grid = fit_mod.build_grid(ws, bands, mh=None, corrections=None)
    ids = list(xp["source_id"])
    ratios = np.full((len(fit), len(xp["wave_nm"])), np.nan)
    for k, row in enumerate(fit.itertuples()):
        j = ids.index(row.source_id)
        ratios[k] = calib.ratio_spectrum(xp["flux"][j].astype(float), xp["flux_err"][j].astype(float),
                                         xp["wave_nm"], grid["m_xp"][int(row.imodel)])
    np.savez_compressed(ws.path("xp_ratios.npz"), source_id=fit.source_id.values, ratio=ratios,
                        wave_nm=xp["wave_nm"])
    print(f"fitted {len(fit)} calibrators at A_V = 0")


def stage_build(ws: Workspace, out: Path) -> None:
    fit = pd.read_csv(ws.path("xp_stars_av0.csv"))
    r = np.load(ws.path("xp_ratios.npz"))
    assert (r["source_id"] == fit.source_id.values).all()
    bands = [c[3:] for c in fit.columns if c.startswith("dm_")]
    # the locked-T_eff fits are poor by construction before the correction exists:
    # keep everything but gross failures
    ok = ((fit.chi2_best * 3 / fit.n_xp < 25.0) & (fit.n_phot >= 5)).values
    print(f"{ok.sum()} of {len(fit)} calibrators used")
    corr = calib.aggregate(fit[ok].reset_index(drop=True), r["ratio"][ok], r["wave_nm"], bands,
                           teff_edges=calib.node_edges(calib.teff_nodes()))
    np.savez_compressed(out, **corr)
    print(f"wrote {out}: bins with data\n{corr['n']}")
    nodes = calib.teff_nodes()
    for t in range(len(nodes)):
        z = 1
        if corr["n"][t, z] and nodes[t] % 500 == 0:
            R = corr["ratio"][t, z]
            print(f"  node {nodes[t]:.0f} K (N {corr['n'][t, z]:3d}): "
                  f"ratio at 400/500/600/900 nm = "
                  + "/".join(f"{np.interp(l, r['wave_nm'], R):.3f}" for l in (400, 500, 600, 900))
                  + "  dm " + " ".join(f"{b.split('_')[-1]}{corr['band_dm'][t, z, j]:+.3f}"
                                       for j, b in enumerate(bands)))


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    ws = workspace()
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else ws.path("template_corrections.npz")
    if stage in ("pull", "all"):
        stage_pull(ws)
    if stage in ("fit", "all"):
        stage_fit(ws)
    if stage in ("build", "all"):
        stage_build(ws, out)


if __name__ == "__main__":
    main()

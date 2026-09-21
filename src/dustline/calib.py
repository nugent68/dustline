"""Empirical template corrections: how a NewEra dwarf template of a given DESI
T_eff / [Fe/H] label differs from real, unreddened stars.

Why: on the COSMOS test the K dwarfs returned A_V ~ 0.2 against a 0.06
foreground.  Against the empirical (IRFM) Mamajek BP-RP locus, DESI T_eff is
50-85 K too hot for K dwarfs and NewEra's K-dwarf colours are too blue at
4500-5000 K (BP-RP 1.347 vs 1.396); pinned at the DESI T_eff, the fit pays for
the colour mismatch with reddening.  A correction indexed by the DESI label
absorbs both at once.

Calibrators: DESI DR1 MWS dwarfs (log g > 4) within 250 pc at |b| > 40 deg
(A_V <~ 0.02) with Gaia XP, PS1, 2MASS and AllWISE.  Each is fitted with the
pipeline at A_V = 0 with its DESI priors; per (T_eff, [Fe/H]) bin we keep the
median XP obs/model ratio spectrum and the median per-band magnitude offsets.
``apply`` multiplies the model grid by them (dwarf models only, log g >= LOGG_MIN).

Product: ``template_corrections.npz`` (release asset; tools/build_template_corrections.py).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

TEFF_EDGES = np.array([3500., 3750., 4000., 4250., 4500., 4750., 5000., 5250., 5500.,
                       5750., 6000., 6500., 7000.])   # (legacy; the bins are the grid nodes)
FEH_EDGES = np.array([-3.0, -0.4, 1.0])      # metal-poor / solar-ish
LOGG_MIN = 3.5                                # dwarf / giant boundary of the log g classes
LOGG_CLASSES = ((3.5, 6.5), (0.0, 3.5))       # 0: dwarfs (BP-RP-locus T_eff), 1: giants (ASPCAP T_eff)
MIN_PER_BIN = 15
TEFF_PRIOR_SIGMA = 1.0    # K: the calibration fits (A_V = 0) AND the science fits LOCK T_eff
                          # to the grid point nearest the DESI label (grid step 100 K), so the
                          # corrections are measured where they are applied (25-100 K priors let
                          # the A_V = 0 fits drift 150-200 K cooler: the XP chi2 gain wins)
_cache: dict = {}


def select_calibrators(max_per_bin: int = 300, dist_pc: float = 250.0, bmin: float = 40.0,
                       snr_min: float = 20.0) -> pd.DataFrame:
    """DESI DR1 MWS dwarfs with Gaia XP within dist_pc at |b| > bmin, capped per
    T_eff bin (highest DESI S/N first).  Data Lab join on Gaia DR3."""
    from .catalogs.datalab import TAP, UA  # noqa: F401
    import io
    import urllib.parse
    import urllib.request

    cols = ("m.source_id, m.teff, m.teff_err, m.logg, m.logg_err, m.feh, m.feh_err, m.alphafe, "
            "m.snr_med, m.program, g.ra, g.dec, g.b, g.parallax, g.phot_g_mean_mag, g.bp_rp")
    adql = (f"SELECT {cols} FROM desi_dr1.mws m JOIN gaia_dr3.gaia_source g ON g.source_id = m.source_id "
            f"WHERE m.rr_spectype = 'STAR' AND m.rvs_warn = 0 AND m.zcat_primary = 't' "
            f"AND m.logg > 4.0 AND m.snr_med > {snr_min} AND g.has_xp_continuous = 1 "
            f"AND g.parallax > {1000.0 / dist_pc} AND g.parallax_over_error > 20 AND g.ruwe < 1.4 "
            f"AND g.ipd_frac_multi_peak <= 10 AND ABS(g.b) > {bmin} AND m.teff BETWEEN 3500 AND 7000")
    q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql))
    txt = urllib.request.urlopen(urllib.request.Request(TAP + "?" + q, headers=UA),
                                 timeout=900).read().decode()
    d = pd.read_csv(io.StringIO(txt)).drop_duplicates("source_id")
    d = d.sort_values("snr_med", ascending=False)
    ib = np.digitize(d.teff, TEFF_EDGES) - 1
    keep = pd.concat([d[ib == k].head(max_per_bin) for k in range(len(TEFF_EDGES) - 1)])
    return keep.sort_values("teff").reset_index(drop=True)


E_BPRP_PER_AV = 0.44      # G23 at R_V 3.1: A_BP/A_V 1.08, A_RP/A_V 0.64
# d(BP-RP)/d[M/H] at fixed T_eff (mag/dex; NewEra dwarfs, log g 4.5, [M/H] -1..0): a
# metal-poor star is bluer, so the solar locus reads it as hotter (COSMOS F/G stars at
# [Fe/H] -0.6 came out +146 K and paid with 0.15 A_V)
_FEH_SLOPE_T = np.array([3500., 4000., 4500., 5000., 5500., 6000., 6500., 7000.])
_FEH_SLOPE = np.array([0.331, 0.119, 0.062, 0.051, 0.035, 0.031, 0.019, 0.008])
_locus = {}


def teff_from_bprp(bprp, av=0.0, feh=0.0) -> np.ndarray:
    """Empirical dwarf T_eff from the dereddened Gaia BP-RP (Mamajek 2022 main-
    sequence locus, packaged in data/mamajek_dwarf_locus.dat), with a
    metallicity term from the NewEra colour sensitivity (_FEH_SLOPE).  This -
    not the DESI label - is the T_eff coordinate of the template corrections:
    the DESI T_eff varies with S/N and brightness at the 100-150 K level
    (calibrators of one DESI label differ by 0.1 mag in BP-RP between G 12 and
    G 14), while BP-RP is robust at every magnitude.  NaN where BP-RP is
    missing or outside the locus."""
    if "locus" not in _locus:
        from . import assets

        c, t = np.loadtxt(assets.filter_dir().parent / "mamajek_dwarf_locus.dat", usecols=(0, 1)).T
        o = np.argsort(c)
        _locus["locus"] = (c[o], t[o])
    c, t = _locus["locus"]
    x = np.asarray(bprp, float) - E_BPRP_PER_AV * np.asarray(av, float)
    z = np.nan_to_num(np.asarray(feh, float), nan=0.0) * np.ones_like(x)
    t0 = np.interp(x, c, t)                                   # solar-locus estimate
    k = np.interp(t0, _FEH_SLOPE_T, _FEH_SLOPE)
    x1 = x - k * z                                            # the colour a solar star of this T_eff has
    out = np.interp(x1, c, t)
    out[~np.isfinite(x1) | (x1 < c.min()) | (x1 > c.max())] = np.nan
    return out


def map_extinction(ra, dec, dist_pc, rv: float = 3.1) -> np.ndarray:
    """A_V of each calibrator from the Edenhofer+2023 3D dust map (dustmaps;
    data under $DUSTLINE_CACHE_DIR/dustmaps).  The map's E is in ZGR23 units:
    A_V = 2.8 E (Edenhofer+2023 Sect. 2).  Stars beyond the map (1.25 kpc)
    take the last distance bin."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from dustmaps.config import config
    from dustmaps.edenhofer2023 import Edenhofer2023Query

    from . import assets

    config["data_dir"] = str(assets.cache_dir() / "dustmaps")
    q = Edenhofer2023Query(integrated=True)
    c = SkyCoord(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg,
                 distance=np.clip(np.asarray(dist_pc, float), 70.0, 1240.0) * u.pc, frame="icrs")
    e = np.asarray(q(c), float)
    return 2.8 * np.nan_to_num(e, nan=0.0)


def deredden(xp: dict, stars: pd.DataFrame, bands: list[str], av: np.ndarray,
             rv: float = 3.1) -> tuple[dict, pd.DataFrame]:
    """Remove a known A_V from each star's XP spectrum (G23 at rv) and photometry
    (band ratios of the reference SED), so the A_V = 0 calibration fit sees the
    intrinsic star."""
    from . import ensemble, extinction

    ext = extinction.curve(np.asarray(xp["wave_nm"]) * 10.0, rv)
    out = dict(xp)
    ids = list(xp["source_id"])
    av_xp = pd.Series(av, index=stars.source_id.values).reindex(ids).fillna(0.0).values
    scale = 10.0 ** (0.4 * av_xp[:, None] * ext[None, :])
    out["flux"] = xp["flux"] * scale
    out["flux_err"] = xp["flux_err"] * scale
    st = stars.copy()
    for b in bands:
        st[f"mag_{b}"] = st[f"mag_{b}"] - av * ensemble.band_ratio_at_rv(b, rv)
    return out, st


def select_giant_calibrators(max_per_node: int = 120, dist_pc: float = 1200.0,
                             snr_min: float = 50.0) -> pd.DataFrame:
    """APOGEE DR17 giants (log g < 3.5, T_eff 3500-5500) with Gaia XP within dist_pc
    at |b| > 30 deg (> 20 deg below 4500 K, where nearby giants are rare), capped
    per 100 K node (highest S/N first).  ASPCAP T_eff (IRFM-calibrated for giants)
    is the T_eff anchor of the giant corrections."""
    import io
    import urllib.parse
    import urllib.request

    from .catalogs.datalab import TAP, UA

    cols = ("a.gaiaedr3_source_id AS source_id, a.teff, a.logg, a.fe_h AS feh, a.fe_h_err AS feh_err, "
            "a.alpha_m AS alphafe, a.snr AS snr_med, a.telescope AS program, "
            "g.ra, g.dec, g.b, g.parallax, g.phot_g_mean_mag, g.bp_rp")
    adql = (f"SELECT {cols} FROM sdss_dr17.apogee2_allstar a JOIN gaia_dr3.gaia_source g "
            f"ON g.source_id = a.gaiaedr3_source_id WHERE a.teff BETWEEN 3500 AND 5500 "
            f"AND a.logg > 0 AND a.logg < {LOGG_MIN} AND a.snr > {snr_min} AND a.fe_h > -3 "
            f"AND g.has_xp_continuous = 1 AND g.parallax > {1000.0 / dist_pc} "
            f"AND g.parallax_over_error > 10 AND g.ruwe < 1.4 AND g.ipd_frac_multi_peak <= 10 "
            f"AND ((ABS(g.b) > 30) OR (ABS(g.b) > 20 AND a.teff < 4500))")
    q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql))
    txt = urllib.request.urlopen(urllib.request.Request(TAP, data=q.encode(), headers=UA),
                                 timeout=900).read().decode()
    d = pd.read_csv(io.StringIO(txt)).drop_duplicates("source_id")
    d["teff_err"] = 0.0
    d["logg_err"] = 0.0
    d = d.sort_values("snr_med", ascending=False)
    nodes = teff_nodes()
    inode = np.abs(nodes[None, :] - d.teff.values[:, None]).argmin(axis=1)
    keep = pd.concat([d[inode == k].head(max_per_node) for k in range(len(nodes))])
    return keep.sort_values("teff").reset_index(drop=True)


def ratio_spectrum(flux, err, wave, model, edge=(340.0, 1015.0)):
    """obs / (C * model) with the scale C profiled over the XP range; NaN outside."""
    ok = np.isfinite(flux) & np.isfinite(err) & (wave >= edge[0]) & (wave <= edge[1]) & (model > 0)
    C = (flux[ok] * model[ok] / err[ok] ** 2).sum() / ((model[ok] ** 2) / err[ok] ** 2).sum()
    out = np.full(len(wave), np.nan)
    out[ok] = flux[ok] / (C * model[ok])
    return out


def teff_nodes() -> np.ndarray:
    """The model cache's T_eff grid points in the fit range: each correction bin is one
    node, i.e. 'stars locked to this template'."""
    from . import models

    meta = models.load_cache()[2]
    return np.array(sorted(set(meta[(meta[:, 0] >= 3200) & (meta[:, 0] <= 12000), 0].tolist())))


def node_edges(nodes: np.ndarray) -> np.ndarray:
    """Bin edges halfway between consecutive nodes (open at the ends)."""
    mid = 0.5 * (nodes[1:] + nodes[:-1])
    return np.concatenate([[nodes[0] - 50.0], mid, [nodes[-1] + 50.0]])


def aggregate(fit: pd.DataFrame, ratios: np.ndarray, xp_wave: np.ndarray, bands: list[str],
              smooth_nm: float = 0.0, teff_edges: np.ndarray | None = None,
              stat: str = "median", clip: float = 3.0) -> dict:
    """Per (T_eff, [Fe/H]) bin medians of the XP ratio spectra (unsmoothed by
    default: the cool-template residuals have 10 nm structure at the TiO band
    heads, common to every star of a node, and the median of >= 15 stars is
    already at the 1 % level; a 20 nm sigma left a -24 % dip at 400 nm in the
    3500 K node which a free A_V then read as extra reddening. smooth_nm > 0
    applies a Gaussian of that sigma. NOT renormalised, so that the XP and band
    corrections share the joint scale of the calibration fit) and of the
    per-band offsets dm_<band>
    (observed - synthetic, mag).  Bins with < MIN_PER_BIN stars fall back to the
    T_eff bin over all [Fe/H], then to no correction."""
    edges = TEFF_EDGES if teff_edges is None else np.asarray(teff_edges, float)
    nT, nZ, nG = len(edges) - 1, len(FEH_EDGES) - 1, len(LOGG_CLASSES)
    R = np.ones((nT, nZ, nG, len(xp_wave)))
    dm = np.zeros((nT, nZ, nG, len(bands)))
    n = np.zeros((nT, nZ, nG), int)
    # bin by the calibrator's BEST-FIT MODEL (T_eff, [M/H], log g), which is how
    # apply() looks the correction up for each model: binning by the star's own
    # [Fe/H] mixed stars fitted with -0.5 and 0.0 templates (whose blue residuals
    # differ by 30 % at 3500 K) into one median that fitted neither
    t = fit.teff_best.values if "teff_best" in fit else fit.teff_spec.values
    z = fit.mh_best.values if "mh_best" in fit else fit.feh_spec.values
    lg = (fit.logg_best.values if "logg_best" in fit
          else fit.logg_spec.values if "logg_spec" in fit else np.full(len(fit), 4.5))
    iT = np.clip(np.digitize(t, edges) - 1, 0, nT - 1)
    iZ = np.clip(np.digitize(z, FEH_EDGES) - 1, 0, nZ - 1)
    iG = np.where(np.nan_to_num(lg.astype(float), nan=4.5) >= LOGG_MIN, 0, 1)
    kern = np.exp(-0.5 * ((np.arange(-30, 31) * 2.0) / smooth_nm) ** 2) if smooth_nm > 0 else None
    if kern is not None:
        kern /= kern.sum()

    def med_ratio(m):
        r = np.nanmedian(ratios[m], axis=0)
        if stat == "mean":
            # sigma-clipped mean about the median: the fit is least squares, so the
            # star-to-star residual it sees averages, and a skewed residual
            # distribution biases the median-template fit (R_V -0.13 for cool dwarfs)
            x = ratios[m]
            mad = 1.4826 * np.nanmedian(np.abs(x - r[None, :]), axis=0)
            x = np.where(np.abs(x - r[None, :]) <= clip * np.maximum(mad, 1e-3)[None, :], x, np.nan)
            r = np.nanmean(x, axis=0)
        good = np.isfinite(r)
        r = np.interp(xp_wave, xp_wave[good], r[good])
        if kern is None:
            return r
        return np.convolve(np.pad(r, 30, mode="edge"), kern, mode="valid")

    for gcls in range(nG):
        for t in range(nT):
            for z in range(nZ):
                m = (iT == t) & (iZ == z) & (iG == gcls)
                if m.sum() < MIN_PER_BIN:
                    m = (iT == t) & (iG == gcls)
                if m.sum() < MIN_PER_BIN:
                    continue
                n[t, z, gcls] = int(m.sum())
                R[t, z, gcls] = med_ratio(m)
                for j, b in enumerate(bands):
                    d = fit.loc[m, f"dm_{b}"].dropna()
                    dm[t, z, gcls, j] = float(d.median()) if len(d) >= MIN_PER_BIN else 0.0
    return dict(teff_edges=edges, feh_edges=FEH_EDGES, xp_wave=xp_wave, ratio=R,
                bands=np.array(bands), band_dm=dm, n=n, logg_min=LOGG_MIN,
                logg_classes=np.array(LOGG_CLASSES))


def load(path=None) -> dict | None:
    """The corrections asset (None when absent or disabled: DUSTLINE_TEMPLATE_CORR=none)."""
    import os

    from . import assets

    key = str(path or os.environ.get("DUSTLINE_TEMPLATE_CORR") or "template_corrections.npz")
    if key.lower() == "none":
        return None
    if key not in _cache:
        from pathlib import Path

        p = Path(key) if Path(key).is_file() else None
        if p is None:
            if key not in assets.ASSETS or assets.ASSETS[key].get("sha256") is None \
                    and not (assets.cache_dir() / key).exists():
                _cache[key] = None
                return None
            p = assets.fetch(key)
        d = np.load(p, allow_pickle=False)
        out = {k: d[k] for k in d.files}
        out["bands"] = [str(b) for b in out["bands"]]
        _cache[key] = out
    return _cache[key]


def tag(corr: dict | None) -> str:
    """Short identity of a corrections table for cache keys."""
    if corr is None:
        return ""
    import hashlib

    return "tc" + hashlib.sha256(corr["ratio"].tobytes() + corr["band_dm"].tobytes()).hexdigest()[:8]


def rv_closure(corr: dict | None, teff, logg) -> np.ndarray:
    """k = A_V * (1/R_V_fit - 1/R_V_true) of the reddening-injection closure test
    (tools/inject_reddening.py closure) for stars of the given fitted T_eff and
    log g; zeros when the table has none.  ensemble.measure_law subtracts k/A_V
    from 1/R_V per star."""
    teff = np.asarray(teff, float)
    if corr is None or "rv_closure_k" not in corr:
        return np.zeros(len(teff))
    k = corr["rv_closure_k"]
    nT = len(corr["teff_edges"]) - 1
    iT = np.clip(np.digitize(teff, corr["teff_edges"]) - 1, 0, nT - 1)
    iG = np.minimum(np.where(np.asarray(logg, float) >= float(corr["logg_min"]), 0, 1), k.shape[1] - 1)
    out = k[iT, iG]
    out[~np.isfinite(teff)] = 0.0
    return out


def apply(corr: dict, meta: np.ndarray, m_xp: np.ndarray, xp_wave: np.ndarray,
          fnu0: np.ndarray, bands: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Multiply the grid products by the corrections of each model's (T_eff, [M/H])
    bin (dwarf models only).  Returns (m_xp, fnu0, corrected_flag)."""
    nT, nZ = len(corr["teff_edges"]) - 1, len(corr["feh_edges"]) - 1
    iT = np.clip(np.digitize(meta[:, 0], corr["teff_edges"]) - 1, 0, nT - 1)
    iZ = np.clip(np.digitize(meta[:, 2], corr["feh_edges"]) - 1, 0, nZ - 1)
    ratio, band_dm, n = corr["ratio"], corr["band_dm"], corr["n"]
    if ratio.ndim == 3:                      # legacy dwarf-only table
        ratio, band_dm, n = ratio[:, :, None], band_dm[:, :, None], n[:, :, None]
    iG = np.where(meta[:, 1] >= float(corr["logg_min"]), 0, 1)
    iG = np.minimum(iG, n.shape[2] - 1)
    has = n[iT, iZ, iG] > 0
    if ratio.shape[2] == 1:                  # legacy: dwarfs only
        has &= meta[:, 1] >= float(corr["logg_min"])
    use = has
    R = np.ones((len(meta), len(xp_wave)), np.float32)
    R[use] = np.stack([np.interp(xp_wave, corr["xp_wave"], ratio[t, z, gc])
                       for t, z, gc in zip(iT[use], iZ[use], iG[use])]).astype(np.float32)
    F = np.ones((len(meta), len(bands)))
    cb = list(corr["bands"])
    for j, b in enumerate(bands):
        if b in cb:
            F[use, j] = 10.0 ** (-0.4 * band_dm[iT[use], iZ[use], iG[use], cb.index(b)])
    return m_xp * R, fnu0 * F, use

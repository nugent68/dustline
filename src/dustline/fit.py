"""Per-star fits of Gaia XP spectra + broadband photometry with PHOENIX NewEra
models x a G23(R_V) extinction curve, the flux scale (R/D)^2 profiled and tied
to the Gaia parallax through a MIST radius prior.

Model:  F_obs(lam) = C * F_NewEra(Teff, logg, [M/H]=0) * 10^(-0.4 A_V ext(lam; R_V)),
        C = (R/D)^2 profiled analytically per grid point; the radius prior
        (log R_impl - log R_MIST(Teff, logg))^2 / sigma^2, R_impl = D sqrt(C),
        D = 1/parallax, breaks the Teff-A_V degeneracy of A/F stars at R ~ 50.
Grid:   NewEra ([M/H] 0, Teff >= 3200 K, all log g) x A_V 0-8 (0.1) x
        R_V 2.3-5.55 (0.25, the G23 range).  With spectroscopic priors
        (DESI MWS T_eff/log g/[Fe/H] on the star table) the [M/H] axis opens
        to -0.5/0/+0.5 and a Gaussian prior on (Teff, logg, [M/H]) is added
        to the per-star chi2 (spec_prior_chi2).

Hard-won pipeline rules baked in (see the dustline method notes):
- never fit XP alone (T_eff biased cool, R_V ~ 4.4): photometry is required;
- [M/H] fixed at 0 unless a per-star spectroscopic prior constrains it (free
  [M/H] absorbs XP/NewEra residual systematics);
- Lindegren+2021 parallax zero point -0.017 mas;
- photometric zero-point offsets iterated, but NIR offsets FROZEN after the
  first pass (iterating them drifts a common JHK shift, R_V +0.04/pass).
"""

from __future__ import annotations

import hashlib
import json
import time

import numpy as np
import pandas as pd

from . import calib, extinction, filters, models
from .cache import Workspace

AV_GRID = np.round(np.concatenate([np.arange(-0.1, 0.5, 0.01), np.arange(0.5, 8.01, 0.1)]), 3)
# 0.01 steps below 0.5 (column mode: a 0.1 grid quantises the posterior once T_eff is locked);
# the grid reaches -0.1 so that the posterior mean of an unreddened star is unbiased (a
# grid starting at 0 pushes every near-zero star to +0.01..0.02 and a field median can
# never fall below 0)
RV_GRID = np.arange(2.3, 5.56, 0.25)
XP_FLOOR = 0.02            # fractional flux systematic per XP sample
XP_PEAK_FLOOR = 0.005      # fraction of the spectrum peak
XP_EDGE = (340.0, 1015.0)  # nm; outside is calibration edge, masked
XP_SCALE = 3.0             # chi2_XP divisor (343 samples ~ 110 resolution elements)
PLX_ZP = -0.017            # Lindegren+2021 global zero point (mas)
MIN_PLX_SNR_PRIOR = 3.0    # below this the radius prior is off
PHOT_SYS_DEFAULT = 0.03    # per-band systematic added in quadrature
NIR_PREFIXES = ("2MASS", "VISTA", "WISE")   # zero points frozen after pass 1
UV_PREFIXES = ("GALEX",)   # zero points iterated PER T_EFF BIN (the model UV bias is
                           # T_eff dependent: ~0.35 mag too bright in NUV for G stars)
UV_SYS = 0.15              # model UV flux systematic (GALEX bands): the NUV residual MAD about
                           # the per-T_eff offsets is 0.15-0.23 mag on COSMOS (chromospheres etc.)
UV_TEFF_EDGES = np.array([5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 8000.0, 12000.0])
UV_MIN_PER_BIN = 12
MH_DEFAULT = (0.0,)
MH_SPECTRO = None          # with spectroscopic priors: every [M/H] the model cache holds
MH_FREE_RANGE = (-0.5, 0.5)   # stars WITHOUT a prior stay on this range (free [M/H] absorbs
                              # the XP/NewEra residuals: COSMOS no-prior stars drifted to -2)
MH_PENALTY = 1e4
SPEC_COLS = ("teff_spec", "teff_spec_err", "logg_spec", "logg_spec_err",
             "feh_spec", "feh_spec_err")
BATCH_MODELS = 64 * 735            # stars x models per batch (memory: ~2.5 GB of cubes)


def _phot_sys(band: str) -> float:
    if band.startswith(UV_PREFIXES):
        return UV_SYS
    return 0.04 if band.endswith("_Y") or band.endswith("_y") else PHOT_SYS_DEFAULT


def build_grid(ws: Workspace, bands: list[str], law: str = "g23",
               mh: tuple[float, ...] = MH_DEFAULT, corrections="default") -> dict:
    """Model products on the fit sub-grid: XP-resolution spectra, band f_nu,
    per-(R_V, A_V) band extinctions, the radius prior.  Cached in the shared
    asset cache keyed by the band list + law (+ the [M/H] axis when it is not
    the solar default, the model cache and the template corrections), reused
    across sightlines.

    corrections: "default" applies the empirical dwarf-template corrections
    (dustline.calib; $DUSTLINE_TEMPLATE_CORR overrides, "none" disables), None
    applies none (used when building the corrections themselves)."""
    from . import assets

    corr = calib.load() if corrections == "default" else corrections
    wave, flux, meta = models.load_cache()
    if mh is None:
        mh = tuple(sorted(set(meta[:, 2].tolist())))
    mh = tuple(float(m) for m in mh)
    tag = "" if mh == MH_DEFAULT else "mh" + ",".join(f"{m:+.1f}" for m in mh)
    ctag = "" if models.cache_name() == "newera_full_cache.npz" else models.cache_tag()
    atag = "" if AV_GRID[0] == 0.0 and len(AV_GRID) == 81 else f"av{AV_GRID[0]:+.2f}:{len(AV_GRID)}"
    key = hashlib.sha256(("|".join(bands) + law + tag + ctag + calib.tag(corr) + atag).encode()).hexdigest()[:10]
    path = assets.cache_dir() / f"xp_grid_{law}_{key}.npz"
    if path.exists():
        d = np.load(path, allow_pickle=False)
        out = {k: d[k] for k in d.files}
        out["bands"] = [str(b) for b in out["bands"]]
        return out
    t0 = time.time()
    sel = models.sub_grid(meta, teff=(3200, 12000), logg=(0.0, 6.0), mh=mh)
    meta, flux = meta[sel], flux[sel]
    xp_wave = np.arange(336.0, 1021.0, 2.0)
    K = models.xp_lsf_matrix(wave, xp_wave)
    m_xp = (flux @ K.T) * models.RD2                       # W m-2 nm-1 at C=(R/D)^2=1
    ext_xp = np.stack([extinction.curve(xp_wave * 10.0, rv, law) for rv in RV_GRID])
    ph = models.Photometry(wave, bands)
    ext_full = np.stack([extinction.curve(wave, rv, law) for rv in RV_GRID])
    fl = flux * 100.0 * models.RD2                         # erg/s/cm2/A at C=1
    fnu0 = np.stack([ph.fnu(fl)[b] for b in bands], axis=1)
    corrected = np.zeros(len(meta), bool)
    if corr is not None:
        m_xp, fnu0, corrected = calib.apply(corr, meta, m_xp, xp_wave, fnu0, bands)
        print(f"  grid: empirical template corrections applied to {corrected.sum()} dwarf models",
              flush=True)
    A_band = np.zeros((len(RV_GRID), len(AV_GRID), len(meta), len(bands)), np.float32)
    for r, rv in enumerate(RV_GRID):
        for a, av in enumerate(AV_GRID):
            if av == 0:
                continue
            fe = ph.fnu(fl * (10.0 ** (-0.4 * av * ext_full[r]))[None, :])
            A_band[r, a] = np.stack(
                [-2.5 * np.log10(fe[b] / fnu0[:, k]) for k, b in enumerate(bands)], axis=1)
        print(f"  grid: band extinction R_V {rv:.2f} ({time.time() - t0:.0f} s)", flush=True)
    # the MIST radius prior at each model's own metallicity (nearest MIST table:
    # -0.5/0/+0.5, so the metal-poor models use the -0.5 isochrones)
    rprior = np.zeros((len(meta), 2))
    for m in np.unique(meta[:, 2]):
        k = meta[:, 2] == m
        rprior[k] = models.radius_prior(meta[k], feh=float(m))
    out = dict(meta=meta, xp_wave=xp_wave, m_xp=m_xp.astype(np.float32),
               ext_xp=ext_xp.astype(np.float32), fnu0=fnu0.astype(np.float64),
               A_band=A_band, rprior=rprior, av=AV_GRID, rv=RV_GRID,
               bands=np.array(bands), corrected=corrected)
    np.savez_compressed(path, **out)
    out["bands"] = bands
    print(f"grid: {len(meta)} models x {len(AV_GRID)} A_V x {len(RV_GRID)} R_V "
          f"-> {path} ({time.time() - t0:.0f} s)", flush=True)
    return out


def _fit_batch(grid, flux, err, plx, plx_err, phot_f, phot_w):
    """chi2 over (R_V, A_V, model) for a batch: (chi2_data, chi2_prior, C)."""
    w = 1.0 / err ** 2
    wd = (w * flux).astype(np.float32).T
    wT = w.astype(np.float32).T
    sdd = (w * flux ** 2).sum(axis=1) / XP_SCALE
    sdd = sdd + (phot_w * phot_f ** 2).sum(axis=1)
    NR, NA, NM = len(grid["rv"]), len(grid["av"]), len(grid["meta"])
    Ns = flux.shape[0]
    A = np.empty((NR, NA, NM, Ns), np.float64)
    B = np.empty_like(A)
    for r in range(NR):
        for a in range(NA):
            m = grid["m_xp"] * (10.0 ** (-0.4 * grid["av"][a] * grid["ext_xp"][r]))[None, :]
            A[r, a] = (m @ wd) / XP_SCALE
            B[r, a] = ((m * m) @ wT) / XP_SCALE
            fnu = grid["fnu0"] * 10.0 ** (-0.4 * grid["A_band"][r, a].astype(np.float64))
            A[r, a] += fnu @ (phot_w * phot_f).T
            B[r, a] += (fnu * fnu) @ phot_w.T
    C = A / B
    chi2 = sdd[None, None, None, :] - A * A / B
    D = 1.0 / plx
    sig_d = 0.4343 * plx_err / plx
    logR = 0.5 * np.log10(np.maximum(C, 1e-30)) + np.log10(D)[None, None, None, :]
    mu = grid["rprior"][:, 0][None, None, :, None]
    sig = np.sqrt(grid["rprior"][:, 1][None, None, :, None] ** 2 + sig_d[None, None, None, :] ** 2)
    chi2_prior = ((logR - mu) / sig) ** 2
    return chi2, chi2_prior, C


def spec_prior_chi2(meta: np.ndarray, spec: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Gaussian spectroscopic prior on (Teff, logg, [M/H]) per (model, star):
    (chi2 (NM, Ns), has_prior (Ns,)).  Stars without a finite prior get zeros.
    The [Fe/H] prior is clamped to the grid's [M/H] range, so a star outside it
    (e.g. [Fe/H] -1.5 on a -0.5..+0.5 grid) sits at the edge without a penalty;
    the caller flags those (mh_clamped)."""
    ns = len(spec)
    out = np.zeros((len(meta), ns))
    if ns == 0 or "teff_spec" not in spec.columns:
        return out, np.zeros(ns, bool)
    t, te = spec.teff_spec.values.astype(float), spec.teff_spec_err.values.astype(float)
    g, ge = spec.logg_spec.values.astype(float), spec.logg_spec_err.values.astype(float)
    z, ze = spec.feh_spec.values.astype(float), spec.feh_spec_err.values.astype(float)
    has = np.isfinite(t) | np.isfinite(g) | np.isfinite(z)     # any of the three terms
    if not has.any():
        return out, has
    te = np.maximum(np.nan_to_num(te, nan=100.0), 1.0)
    ge = np.maximum(np.nan_to_num(ge, nan=0.3), 0.05)
    ze = np.maximum(np.nan_to_num(ze, nan=0.3), 0.05)
    zc = np.clip(np.nan_to_num(z, nan=0.0), meta[:, 2].min(), meta[:, 2].max())
    chi = np.where(np.isfinite(t)[None, :], ((meta[:, 0][:, None] - np.nan_to_num(t)[None, :]) / te[None, :]) ** 2, 0.0)
    chi += np.where(np.isfinite(g)[None, :], ((meta[:, 1][:, None] - np.nan_to_num(g)[None, :]) / ge[None, :]) ** 2, 0.0)
    chi += np.where(np.isfinite(z)[None, :], ((meta[:, 2][:, None] - zc[None, :]) / ze[None, :]) ** 2, 0.0)
    out[:, has] = chi[:, has]
    return out, has


def _summarize(grid, chi2, chi2_prior, C, bands):
    """Per-star best-fit + posterior summaries from the chi2 cubes."""
    tot = chi2 + chi2_prior
    meta = grid["meta"]
    teff, logg, mh = meta[:, 0], meta[:, 1], meta[:, 2]
    rows = []
    for s in range(tot.shape[-1]):
        t = tot[..., s]
        k = np.unravel_index(np.argmin(t), t.shape)
        L = np.exp(-0.5 * (t - t[k]))
        L /= L.sum()
        pr = L.sum(axis=(0, 1))
        pav = L.sum(axis=(0, 2))
        prv = L.sum(axis=(1, 2))

        def mstd(x, p):
            m = (x * p).sum()
            return m, np.sqrt(max((x * x * p).sum() - m * m, 0.0))

        tm, ts = mstd(teff, pr)
        gm, gs = mstd(logg, pr)
        zm, zs = mstd(mh, pr)
        am, as_ = mstd(grid["av"], pav)
        rm, rs = mstd(grid["rv"], prv)
        Ab = np.einsum("ram,ramb->b", L, grid["A_band"])
        Ab2 = np.einsum("ram,ramb->b", L, grid["A_band"] ** 2)
        Abs = np.sqrt(np.maximum(Ab2 - Ab ** 2, 0.0))
        row = dict(teff_best=teff[k[2]], logg_best=logg[k[2]], mh_best=mh[k[2]], imodel=k[2],
                   av_best=grid["av"][k[1]], rv_best=grid["rv"][k[0]],
                   chi2_best=chi2[k[0], k[1], k[2], s],
                   chi2_prior_best=chi2_prior[k[0], k[1], k[2], s],
                   C_best=C[k[0], k[1], k[2], s],
                   teff=tm, teff_err=ts, logg=gm, logg_err=gs, mh=zm, mh_err=zs,
                   av=am, av_err=as_, rv=rm, rv_err=rs)
        for j, b in enumerate(bands):
            row[f"A_{b}"] = Ab[j]
            row[f"A_{b}_err"] = Abs[j]
        rows.append(row)
    return pd.DataFrame(rows)


def fit_stars(ws: Workspace, stars: pd.DataFrame, xp, bands: list[str],
              offsets: dict[str, float] | None = None, law: str = "g23",
              batch: int = 0, max_stars: int = 0, spectro_priors: bool = True,
              rv_fixed: float | None = None,
              star_offsets: pd.DataFrame | None = None, av_fixed: float | None = None,
              corrections="default") -> pd.DataFrame:
    """Fit every XP star; returns the per-star results table.

    With spectro_priors and ``teff_spec``/``logg_spec``/``feh_spec`` columns on
    ``stars`` (DESI MWS, see catalogs.desi) the grid opens to [M/H] -0.5..+0.5
    and each star with a prior gets spec_prior_chi2 added to its chi2.
    rv_fixed: restrict the R_V axis to the nearest grid value (column mode: at
    A_V ~ 0.1 R_V is unconstrained and only adds noise).
    star_offsets: per-star zero-point offsets (source_id + off_<band> columns,
    the T_eff-binned UV offsets) applied on top of the global ``offsets``.
    av_fixed: hold A_V at the nearest grid value (calibration fits at A_V = 0).
    corrections: see build_grid.
    """
    use_spec = spectro_priors and all(c in stars.columns for c in SPEC_COLS) \
        and (np.isfinite(stars.teff_spec).any() or np.isfinite(stars.feh_spec).any())
    grid = build_grid(ws, bands, law, mh=MH_SPECTRO if use_spec else MH_DEFAULT,
                      corrections=corrections)
    if rv_fixed is not None:
        r = int(np.argmin(np.abs(grid["rv"] - rv_fixed)))
        grid = dict(grid, rv=grid["rv"][r:r + 1], ext_xp=grid["ext_xp"][r:r + 1],
                    A_band=grid["A_band"][r:r + 1])
    if av_fixed is not None:
        a = int(np.argmin(np.abs(grid["av"] - av_fixed)))
        grid = dict(grid, av=grid["av"][a:a + 1], A_band=grid["A_band"][:, a:a + 1])
    batch = batch or max(8, BATCH_MODELS // len(grid["meta"]))
    offsets = {} if offsets is None else dict(offsets)
    offsets = {b: offsets.get(b, 0.0) for b in bands}
    zp = np.array([filters.get(b).zp_jy for b in bands]) * 1e-23
    ids = xp["source_id"]
    g = stars.set_index("source_id").loc[ids].reset_index()
    if star_offsets is not None and len(star_offsets):
        so = star_offsets.set_index("source_id").reindex(g.source_id).fillna(0.0)
        for b in bands:
            if f"off_{b}" in so:
                g[f"mag_{b}"] = g[f"mag_{b}"].values - so[f"off_{b}"].values
    flux_all = xp["flux"].astype(np.float64)
    err_all = xp["flux_err"].astype(np.float64)
    wave = xp["wave_nm"]
    n = len(g) if not max_stars else min(max_stars, len(g))
    n_spec = int((np.isfinite(g.teff_spec.iloc[:n]) | np.isfinite(g.logg_spec.iloc[:n])
                  | np.isfinite(g.feh_spec.iloc[:n])).sum()) if use_spec else 0
    print(f"fitting {n} stars, {len(grid['meta'])} models, bands {bands}"
          + (f", spectroscopic priors for {n_spec}" if use_spec else "")
          + (f", R_V fixed at {grid['rv'][0]:.2f}" if rv_fixed is not None else "")
          + (f", A_V fixed at {grid['av'][0]:.2f}" if av_fixed is not None else "")
          + (f", {int(grid['corrected'].sum())} corrected templates" if "corrected" in grid
             and grid["corrected"].any() else ""), flush=True)
    mh_lo, mh_hi = grid["meta"][:, 2].min(), grid["meta"][:, 2].max()
    out = []
    t0 = time.time()
    for k0 in range(0, n, batch):
        sl = slice(k0, min(k0 + batch, n))
        flux, err = flux_all[sl].copy(), err_all[sl].copy()
        peak = np.nanmax(np.abs(flux), axis=1, keepdims=True)
        err = np.sqrt(err ** 2 + (XP_FLOOR * np.abs(flux)) ** 2 + (XP_PEAK_FLOOR * peak) ** 2)
        bad = (~np.isfinite(flux) | ~np.isfinite(err)
               | (wave < XP_EDGE[0])[None, :] | (wave > XP_EDGE[1])[None, :])
        flux[bad] = 0.0
        err[bad] = np.inf
        gs = g.iloc[sl]
        plx = gs.parallax.values - PLX_ZP
        plx_err = gs.parallax_error.values
        good = np.isfinite(plx) & (plx / plx_err > MIN_PLX_SNR_PRIOR)
        plx_fit = np.where(good, plx, 1.0)
        plxe_fit = np.where(good, plx_err, 1.0)
        mags = np.stack([gs[f"mag_{b}"].values - offsets[b] for b in bands], axis=1)
        merr = np.stack([np.sqrt(gs[f"magerr_{b}"].values ** 2 + _phot_sys(b) ** 2)
                         for b in bands], axis=1)
        ok = np.isfinite(mags) & np.isfinite(merr) & (merr < 0.5)
        phot_f = np.where(ok, zp * 10.0 ** (-0.4 * np.nan_to_num(mags)), 0.0)
        phot_w = np.where(ok, 1.0 / (0.9210 * merr * phot_f) ** 2, 0.0)
        phot_w[~ok] = 0.0
        chi2, chi2_prior, C = _fit_batch(grid, flux, err, plx_fit, plxe_fit, phot_f, phot_w)
        chi2_prior[..., ~good] = 0.0
        if use_spec:
            chi2_spec, has_spec = spec_prior_chi2(grid["meta"], gs)
            # no [M/H] prior: keep the star on the solar-ish [M/H] range
            has_feh = np.isfinite(gs.feh_spec.values.astype(float))
            outside = (grid["meta"][:, 2] < MH_FREE_RANGE[0]) | (grid["meta"][:, 2] > MH_FREE_RANGE[1])
            chi2_spec[np.ix_(outside, ~has_feh)] = MH_PENALTY
            chi2_prior += chi2_spec[None, None, :, :]
        else:
            has_spec = np.zeros(len(gs), bool)
        df = _summarize(grid, chi2, chi2_prior, C, bands)
        df.insert(0, "source_id", gs.source_id.values)
        df["n_xp"] = np.isfinite(err).sum(axis=1)
        df["plx_used"] = good
        df["n_phot"] = (phot_w > 0).sum(axis=1)
        df["spec_prior"] = has_spec
        if use_spec:
            z = gs.feh_spec.values.astype(float)
            df["mh_clamped"] = has_spec & ((z < mh_lo - 1e-6) | (z > mh_hi + 1e-6))
        # observed - synthetic magnitude at the best fit, for the zero-point iteration
        for j, b in enumerate(bands):
            im = df.imodel.values
            ia = np.abs(grid["av"][None, :] - df.av_best.values[:, None]).argmin(axis=1)
            ir = np.abs(grid["rv"][None, :] - df.rv_best.values[:, None]).argmin(axis=1)
            syn = -2.5 * np.log10(grid["fnu0"][im, j]
                                  * 10 ** (-0.4 * grid["A_band"][ir, ia, im, j])
                                  * df.C_best.values / zp[j])
            df[f"dm_{b}"] = mags[:, j] - syn
        out.append(df)
        print(f"  {sl.stop}/{n} ({time.time() - t0:.0f} s)", flush=True)
    res = pd.concat(out, ignore_index=True)
    keep = ["source_id", "ra", "dec", "parallax", "parallax_error", "parallax_over_error",
            "ruwe", "ipd_frac_multi_peak", "phot_g_mean_mag", "bp_rp"]
    if use_spec:
        keep += [c for c in SPEC_COLS + ("spec_snr", "teff_desi", "teff_ap", "logg_ap", "feh_ap", "ap_snr")
                 if c in g.columns and c not in keep]
    res = g[keep].iloc[:n].merge(res, on="source_id")
    res["D_kpc"] = 1.0 / (res.parallax - PLX_ZP)
    res["D_err"] = res.parallax_error / (res.parallax - PLX_ZP) ** 2
    return res


def measure_offsets(stars_fit: pd.DataFrame, bands: list[str], prev: dict[str, float],
                    min_bands: int, freeze_nir: bool,
                    freeze_all: bool = False) -> tuple[dict[str, float], float]:
    """Cumulative per-band zero-point offsets (observed - synthetic at the fit) from
    well-fit stars; returns (offsets, largest optical move).  NIR bands are frozen
    after the first pass (freeze_nir=True) - iterating them drifts R_V.  With
    freeze_all the offsets are measured (the largest one is returned as the
    "move") but not applied: prev is returned unchanged."""
    r = stars_fit
    ok = ((r.chi2_best * 3 / r.n_xp < 2.5) & (r.ruwe < 1.4)
          & (r.ipd_frac_multi_peak <= 10) & (r.n_phot >= min_bands))
    r = r[ok]
    off = dict(prev)
    max_move = 0.0
    for b in bands:
        d = r[f"dm_{b}"].dropna()
        if len(d) < 20:
            off.setdefault(b, 0.0)
            continue
        med = float(d.median())
        if b.startswith(UV_PREFIXES):
            off.setdefault(b, 0.0)       # UV: handled per T_eff bin (uv_offsets_by_teff)
            continue
        is_nir = any(b.startswith(p) for p in NIR_PREFIXES)
        if freeze_all:
            off.setdefault(b, 0.0)
            max_move = max(max_move, abs(med))
            continue
        if is_nir and freeze_nir and b in prev:
            continue    # NIR offsets frozen at their first-pass values
        off[b] = round(prev.get(b, 0.0) + med, 4)
        if not is_nir:
            max_move = max(max_move, abs(med))
    return off, max_move


def uv_offsets_by_teff(stars_fit: pd.DataFrame, bands: list[str], prev: dict,
                       min_bands: int) -> tuple[dict, pd.DataFrame, float]:
    """Cumulative UV zero-point offsets per T_eff bin (UV_TEFF_EDGES) from well-fit
    stars, applied per star through its fitted T_eff (or its spectroscopic T_eff
    when it has one).  Returns (table {band: {edges, offsets}}, per-star DataFrame
    (source_id, off_<band>), largest move).  Bins with < UV_MIN_PER_BIN stars take
    the global median.  The NewEra UV flux is too bright for G stars (~0.35 mag
    in NUV at 5800 K) and the bias depends on T_eff, so a global offset - let
    alone a frozen one - would leak into A_V."""
    uv = [b for b in bands if b.startswith(UV_PREFIXES)]
    r = stars_fit
    teff = np.where(np.isfinite(r.get("teff_spec", pd.Series(np.nan, index=r.index))),
                    r.get("teff_spec", pd.Series(np.nan, index=r.index)), r.teff)
    ok = ((r.chi2_best * 3 / r.n_xp < 2.5) & (r.ruwe < 1.4)
          & (r.ipd_frac_multi_peak <= 10) & (r.n_phot >= min_bands)).values
    table, per_star, max_move = dict(prev), pd.DataFrame(dict(source_id=r.source_id.values)), 0.0
    ibin = np.clip(np.digitize(teff, UV_TEFF_EDGES) - 1, 0, len(UV_TEFF_EDGES) - 2)
    for b in uv:
        d = r[f"dm_{b}"].values
        prev_off = np.array(prev.get(b, {}).get("offsets", [0.0] * (len(UV_TEFF_EDGES) - 1)))
        good = ok & np.isfinite(d)
        glob = float(np.median(d[good])) if good.sum() >= UV_MIN_PER_BIN else 0.0
        moves = np.zeros(len(prev_off))
        for k in range(len(prev_off)):
            m = good & (ibin == k)
            moves[k] = float(np.median(d[m])) if m.sum() >= UV_MIN_PER_BIN else glob
        new = np.round(prev_off + moves, 4)
        table[b] = dict(edges=UV_TEFF_EDGES.tolist(), offsets=new.tolist(),
                        n=[int((good & (ibin == k)).sum()) for k in range(len(new))])
        per_star[f"off_{b}"] = new[ibin]
        max_move = max(max_move, float(np.abs(moves[np.array(table[b]["n"]) >= UV_MIN_PER_BIN]).max()
                                       if any(np.array(table[b]["n"]) >= UV_MIN_PER_BIN) else 0.0))
    return table, per_star, max_move


def run_fit_with_offsets(ws: Workspace, stars: pd.DataFrame, xp, bands: list[str],
                         law: str = "g23", min_bands_offsets: int = 4,
                         max_passes: int = 3, converge: float = 0.02,
                         force: bool = False, spectro_priors: bool = True,
                         freeze_offsets: bool = False,
                         rv_fixed: float | None = None,
                         teff_from_colour: bool = False,
                         desi_teff: bool = False) -> pd.DataFrame:
    """The fit + zero-point iteration: fit, measure offsets, refit until the
    optical offsets move < converge mag (NIR frozen after pass 1). Cached.

    freeze_offsets: single pass with all zero points held at 0 (the offsets are
    still measured and written to phot_offsets.json for the record) - the A/B
    control for low-extinction fields, where a uniform A_V screen and the
    optical zero points are partly degenerate.
    teff_from_colour (column mode): every star's T_eff prior is the empirical
    dwarf T_eff of its Gaia BP-RP (calib.teff_from_bprp), dereddened by the
    previous pass's A_V (0 in pass 1; converges in two passes at A_V ~ 0.05),
    locked like the DESI label would be; DESI still supplies log g / [Fe/H].
    """
    out_path = ws.path("xp_stars.csv")
    off_path = ws.path("phot_offsets.json")
    if out_path.exists() and not force:
        return pd.read_csv(out_path)
    offsets: dict[str, float] = {}
    uv_table: dict = {}
    star_off = None
    has_uv = any(b.startswith(UV_PREFIXES) for b in bands)
    fit = None
    stars = stars.copy()
    for c in SPEC_COLS:
        if c not in stars:
            stars[c] = np.nan
    if "teff_desi" not in stars:
        stars["teff_desi"] = stars["teff_spec"]
    if not desi_teff and not teff_from_colour:
        # the DESI T_eff label is not a 100 K anchor (S/N-dependent by +/-200 K against
        # the colour scale): log g / [Fe/H] priors only; T_eff from the data + radius prior
        stars["teff_spec"] = np.nan
    n_pass = 1 if freeze_offsets else max_passes + (1 if has_uv else 0)
    for p in range(max(n_pass, 2 if teff_from_colour else 1)):
        if teff_from_colour:
            av_prev = (fit.set_index("source_id").av.reindex(stars.source_id).fillna(0.0).values
                       if fit is not None else np.zeros(len(stars)))
            feh = stars.feh_spec.values.astype(float)
            feh = np.where(np.isfinite(feh), feh, np.nanmedian(feh) if np.isfinite(feh).any() else 0.0)
            t_col = calib.teff_from_bprp(stars.bp_rp.values, av_prev, feh)
            stars["teff_spec"] = t_col
            stars["teff_spec_err"] = calib.TEFF_PRIOR_SIGMA
            print(f"    T_eff prior from BP-RP (dereddened by pass-{p} A_V) for "
                  f"{int(np.isfinite(t_col).sum())} stars", flush=True)
        print(f"--- fit pass {p + 1} (offsets: { {k: v for k, v in offsets.items()} }"
              + (f"; UV bins {uv_table}" if uv_table else "") + ")", flush=True)
        fit = fit_stars(ws, stars, xp, bands, offsets=offsets, law=law,
                        spectro_priors=spectro_priors, rv_fixed=rv_fixed, star_offsets=star_off)
        offsets, moved = measure_offsets(fit, bands, offsets, min_bands_offsets,
                                         freeze_nir=(p >= 1), freeze_all=freeze_offsets)
        if has_uv and not freeze_offsets:
            uv_table, star_off, uv_moved = uv_offsets_by_teff(fit, bands, uv_table, min_bands_offsets)
            moved = max(moved, uv_moved)
        json.dump(dict(offsets, **{b: v for b, v in uv_table.items()}), open(off_path, "w"), indent=1)
        print(f"    offsets now {offsets} (largest move {moved:.3f})"
              + (f"\n    UV offsets by T_eff: {uv_table}" if uv_table else ""), flush=True)
        if p >= 1 and moved < converge:
            break
    if star_off is not None:
        fit = fit.merge(star_off, on="source_id", how="left")
    if "teff_desi" not in fit:
        fit = fit.merge(stars[["source_id", "teff_desi"]], on="source_id", how="left")
    fit.round(5).to_csv(out_path, index=False)
    print(f"wrote {out_path}: {len(fit)} stars", flush=True)
    return fit

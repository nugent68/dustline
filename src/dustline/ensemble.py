"""Ensemble products from the per-star fits: the sightline extinction law
(R_V, band ratios) and the extinction-distance run A_V(D) from the Gaia
parallax stars, with an optional red-clump bulge bridge; at high latitude
(column mode) the assumed-law dict and the foreground column instead.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import extinction, models

D_FINE = np.arange(0.1, 12.01, 0.1)
MIN_AV_LAW = 2.0        # A_V above which a star enters the law average
RV_ASSUMED = 3.1        # column mode: the law is not measurable, G23 at this R_V is assumed
D_BEHIND_KPC = 0.5      # column mode: stars beyond this are behind all the local dust
TEFF_CLEAN = 5500.0     # column mode: K/M dwarfs return spurious A_V ~ 0.2 (template/XP
                        # systematics at low A_V; COSMOS test), F/G stars do not
CHI2N_MAX = 2.5
RUWE_MAX = 1.4
IPD_MAX = 10


def good_sample(r: pd.DataFrame, min_snr: float = 5.0, min_bands: int = 4) -> pd.Series:
    """Quality cuts for the dust-run sample (needs a parallax distance)."""
    return ((r.parallax_over_error > min_snr) & r.plx_used & (r.ruwe < RUWE_MAX)
            & (r.ipd_frac_multi_peak <= IPD_MAX) & (r.chi2_best * 3 / r.n_xp < CHI2N_MAX)
            & (r.n_phot >= min_bands) & (r.av < 7.9))


def law_sample(r: pd.DataFrame, min_snr: float = 2.0, min_bands: int = 4,
               min_av: float = MIN_AV_LAW) -> pd.Series:
    """Quality cuts for the LAW sample (R_V needs reddening, not distance).
    Below plx S/N 3 the radius prior is off, so those stars must be cool."""
    return ((r.parallax_over_error > min_snr) & (r.ruwe < RUWE_MAX)
            & (r.ipd_frac_multi_peak <= IPD_MAX) & (r.chi2_best * 3 / r.n_xp < CHI2N_MAX)
            & (r.n_phot >= min_bands) & (r.av < 7.9) & (r.av >= min_av)
            & (r.plx_used | (r.teff < 5500)))


def running_profile(D: np.ndarray, A: np.ndarray, frac: float = 0.25,
                    nmin: int = 12) -> np.ndarray:
    """Running median/16/84 % of A in a +/-frac*D window: (len(D_FINE), 4)."""
    out = np.full((len(D_FINE), 4), np.nan)
    for k, d in enumerate(D_FINE):
        m = np.abs(D - d) <= frac * d
        if m.sum() >= nmin:
            out[k] = np.percentile(A[m], [50, 16, 84]).tolist() + [m.sum()]
    return out


def measure_law(fit: pd.DataFrame, bands: list[str], min_snr: float = 2.0,
                min_bands: int = 4, min_av: float = MIN_AV_LAW, closure: bool = True) -> dict:
    """The sightline law: R_V statistics + per-star band ratios A_b/A_V.

    closure: subtract the reddening-injection closure offset of each star's
    template node (calib.rv_closure: k/A_V in 1/R_V; cool dwarfs and giants read
    R_V ~0.25 low at A_V = 1) - rv_raw keeps the uncorrected median."""
    from . import calib

    w = fit[law_sample(fit, min_snr, min_bands, min_av)]
    if len(w) < 10:
        raise RuntimeError(
            f"only {len(w)} stars with A_V >= {min_av} pass the law cuts - "
            "not enough reddening on this sightline to measure R_V; "
            "the A(D) run is still produced (at the G23 default shape)")
    rv_raw = w.rv.values
    k = calib.rv_closure(calib.load(), w.teff.values, w.logg.values) if closure else np.zeros(len(w))
    rv = np.where(k != 0, 1.0 / np.clip(1.0 / rv_raw - k / np.maximum(w.av.values, 0.3), 1e-3, None), rv_raw)
    rve = np.maximum(w.rv_err.values, 0.05)
    wm = 1.0 / rve ** 2
    rv_med = float(np.median(rv))
    rv_mad = float(1.4826 * np.median(np.abs(rv - rv_med)))
    ratios, ratios_mad = {}, {}
    for b in bands:
        q = (w[f"A_{b}"] / w["av"]).values     # A_band / A_V per star
        if closure and (k != 0).any():
            # move each corrected star's ratio along the law's R_V dependence
            rg, qg = _ratio_curve(b)
            q = q * np.interp(rv, rg, qg) / np.interp(rv_raw, rg, qg)
        med = float(np.median(q))
        ratios[b] = med
        ratios_mad[b] = float(1.4826 * np.median(np.abs(q - med)))
    return dict(rv=rv_med, rv_mad=rv_mad,
                rv_mean=float((wm * rv).sum() / wm.sum()),
                rv_mean_err=float(1.0 / np.sqrt(wm.sum())),
                rv_raw=float(np.median(rv_raw)), n_closure=int((k != 0).sum()),
                n_stars=int(len(w)), min_av=float(min_av),
                ratios_av=ratios, ratios_av_mad=ratios_mad,
                rv_16=float(np.percentile(rv, 16)), rv_84=float(np.percentile(rv, 84)),
                systematics="NIR zero-point freeze +/-0.05 on R_V; XP blue-end "
                            "systematics for faint heavily-reddened stars (see docs)")


def default_law(bands: list[str], rv: float = RV_ASSUMED, reason: str = "") -> dict:
    """The law dict when R_V cannot be measured: G23 at the assumed R_V, band
    ratios from the reference red-giant SED (band_ratio_at_rv). n_stars = 0 and
    rv_assumed = True mark it; rv_mad / rv_16 / rv_84 are 0."""
    ratios = {b: band_ratio_at_rv(b, rv) for b in bands}
    return dict(rv=float(rv), rv_mad=0.0, rv_mean=float(rv), rv_mean_err=0.0, n_stars=0,
                min_av=float(MIN_AV_LAW), ratios_av=ratios,
                ratios_av_mad={b: 0.0 for b in bands}, rv_16=float(rv), rv_84=float(rv),
                rv_assumed=True, mode="column",
                systematics="R_V not measured (no A_V >= 2 stars): G23 at R_V = "
                            f"{rv:g} assumed" + (f"; {reason}" if reason else ""))


def foreground_column(fit: pd.DataFrame, d_min_kpc: float = D_BEHIND_KPC,
                      min_snr: float = 5.0, min_bands: int = 4, ncell: int = 4,
                      n_boot: int = 300, seed: int = 0) -> dict:
    """The total foreground A_V of a high-latitude field from the parallax
    stars behind the local dust (D > d_min_kpc): median, MAD, bootstrap error
    of the median, inverse-variance mean, N; the same for the subsets with and
    without a spectroscopic prior, for T_eff above/below TEFF_CLEAN and by G
    magnitude; and an ncell x ncell map of the median A_V across the field
    (tangent plane, cells ordered by increasing x then y).

    ``clean`` (T_eff >= TEFF_CLEAN) is the recommended column: on the COSMOS
    test field the K/M dwarfs (mostly G > 15.5) return A_V ~ 0.2-0.3 against
    an SFD column of 0.05 - a cool-template / faint-XP systematic that the
    F/G stars (0.07-0.09) do not show. ``all`` keeps everything for reference."""
    s = fit[good_sample(fit, min_snr, min_bands) & (fit.D_kpc > d_min_kpc)]
    has_spec = s["spec_prior"].astype(bool) if "spec_prior" in s else pd.Series(False, index=s.index)

    def stats(a: pd.DataFrame) -> dict:
        if len(a) == 0:
            return dict(n=0, av=np.nan, av_mad=np.nan, av_err=np.nan, av_wmean=np.nan)
        av = a.av.values
        med = float(np.median(av))
        mad = float(1.4826 * np.median(np.abs(av - med)))
        rng = np.random.default_rng(seed)
        boot = [np.median(rng.choice(av, len(av))) for _ in range(n_boot)] if len(av) > 3 else [med]
        w = 1.0 / np.maximum(a.av_err.values, 0.03) ** 2
        return dict(n=int(len(a)), av=med, av_mad=mad, av_err=float(np.std(boot)),
                    av_wmean=float((w * av).sum() / w.sum()))

    out = dict(d_min_kpc=float(d_min_kpc), teff_clean=float(TEFF_CLEAN), **stats(s))
    hot = s.teff >= TEFF_CLEAN
    out["clean"] = stats(s[hot])
    out["clean_with_spec_prior"] = stats(s[hot & has_spec])
    if "dm_GALEX_NUV" in s:
        has_nuv = np.isfinite(s["dm_GALEX_NUV"])
        out["clean_with_nuv"] = stats(s[hot & has_nuv])
        out["clean_without_nuv"] = stats(s[hot & ~has_nuv])
    out["with_spec_prior"] = stats(s[has_spec])
    out["without_spec_prior"] = stats(s[~has_spec])
    out["teff_hot"] = out["clean"]
    out["teff_cool"] = stats(s[~hot])
    gmag = s["phot_g_mean_mag"] if "phot_g_mean_mag" in s else pd.Series(np.nan, index=s.index)
    out["by_gmag"] = {f"{lo:g}-{hi:g}": stats(s[(gmag >= lo) & (gmag < hi)])
                      for lo, hi in ((8, 14), (14, 15.5), (15.5, 16.5), (16.5, 18))}
    if has_spec.any():
        a = s[has_spec]
        out["teff_spec_minus_fit"] = float(np.median(a.teff_spec - a.teff))
        out["n_mh_clamped"] = int(a["mh_clamped"].sum()) if "mh_clamped" in a else 0
    # spatial map: tangent-plane cells over the field (clean stars)
    s = s[hot]
    if len(s) and ncell > 1:
        ra0, dec0 = float(fit.ra.median()), float(fit.dec.median())
        x = (s.ra - ra0) * np.cos(np.radians(dec0))
        y = s.dec - dec0
        r = float(np.hypot(x, y).max())
        edges = np.linspace(-r, r, ncell + 1)
        ix = np.clip(np.digitize(x, edges) - 1, 0, ncell - 1)
        iy = np.clip(np.digitize(y, edges) - 1, 0, ncell - 1)
        cells = []
        for j in range(ncell):
            for i in range(ncell):
                m = (ix == i) & (iy == j)
                cells.append(dict(x=float(0.5 * (edges[i] + edges[i + 1])),
                                  y=float(0.5 * (edges[j] + edges[j + 1])),
                                  n=int(m.sum()),
                                  av=float(np.median(s.av[m])) if m.sum() >= 5 else None))
        out["map"] = dict(ra0=ra0, dec0=dec0, half_width_deg=r, ncell=ncell, cells=cells)
    return out


def dust_run(fit: pd.DataFrame, min_snr: float = 5.0, min_bands: int = 4) -> pd.DataFrame:
    """A_V(D) running profile of the parallax stars: D_kpc, A_V median/16/84, N."""
    s = fit[good_sample(fit, min_snr, min_bands)]
    if len(s) < 20:
        raise RuntimeError(f"only {len(s)} stars pass the run cuts (need >= 20); "
                           "try a larger radius")
    prof = running_profile(s.D_kpc.values, s.av.values)
    df = pd.DataFrame(dict(D_kpc=np.round(D_FINE, 2), AV_med=prof[:, 0],
                           AV_16=prof[:, 1], AV_84=prof[:, 2], n_stars=prof[:, 3]))
    df = df[np.isfinite(df.AV_med)].reset_index(drop=True)
    df["bridged"] = False
    return df


def bridge_to_clump(run: pd.DataFrame, anchor: dict) -> pd.DataFrame:
    """Extend the A_V(D) run beyond the last parallax bin: linear bridge to the
    red-clump bulge column at D_RC, held flat behind it (bridged rows flagged,
    n_stars = 0). When the parallax run already reaches D_RC, the run is
    extended flat at the clump column (short 0.5 kpc ramp from the last bin).
    anchor: the dict from clump.find_clump()."""
    D1 = float(run.D_kpc.iloc[-1])
    A1, lo1, hi1 = (float(run[c].iloc[-1]) for c in ("AV_med", "AV_16", "AV_84"))
    Drc, Arc = anchor["D_RC"], anchor["AV_column"]
    sig = anchor["AV_column_err"]
    ext = D_FINE[D_FINE > D1]
    if not len(ext):
        return run
    Dj = max(Drc, D1 + 0.5)      # junction: the clump, or a short ramp past the data
    f = np.clip((ext - D1) / max(Dj - D1, 0.1), 0.0, 1.0)
    mu = A1 + f * (Arc - A1)
    lo = lo1 + f * (Arc - sig - lo1)
    hi = hi1 + f * (Arc + sig - hi1)
    tail = pd.DataFrame(dict(D_kpc=np.round(ext, 2), AV_med=mu, AV_16=lo, AV_84=hi,
                             n_stars=0, bridged=True))
    return pd.concat([run, tail], ignore_index=True)


_ratio_curves: dict = {}


def _ratio_curve(band: str) -> tuple[np.ndarray, np.ndarray]:
    """band_ratio_at_rv tabulated on R_V 2.0-6.0 (memoised)."""
    if band not in _ratio_curves:
        rg = np.linspace(2.3, 5.6, 34)           # the G23 validity range
        _ratio_curves[band] = (rg, np.array([band_ratio_at_rv(band, r) for r in rg]))
    return _ratio_curves[band]


def band_ratio_at_rv(band: str, rv: float, teff: float = 4500.0, logg: float = 2.5,
                     law: str = "g23", av: float = 2.0) -> float:
    """Band-integrated A_band/A_V for a NewEra reference SED under G23(R_V).

    This converts the measured R_V into A_X for ANY filter X with a curve -
    the output side of the pipeline (Cousins I by default in the API).
    """
    wave, flux, meta = models.load_cache()
    d2 = (meta[:, 0] - teff) ** 2 / 500.0 ** 2 + (meta[:, 1] - logg) ** 2 + meta[:, 2] ** 2 * 100
    k = int(np.argmin(d2))
    ph = models.Photometry(wave, [band])
    A = extinction.band_extinction(ph, flux[k] * 100.0, av, rv, law)
    return float(A[band] / av)

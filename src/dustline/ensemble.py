"""Ensemble products from the per-star fits: the sightline extinction law
(R_V, band ratios) and the extinction-distance run A_V(D) from the Gaia
parallax stars, with an optional red-clump bulge bridge.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import extinction, models

D_FINE = np.arange(0.1, 12.01, 0.1)
MIN_AV_LAW = 2.0        # A_V above which a star enters the law average
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
                min_bands: int = 4, min_av: float = MIN_AV_LAW) -> dict:
    """The sightline law: R_V statistics + per-star band ratios A_b/A_V."""
    w = fit[law_sample(fit, min_snr, min_bands, min_av)]
    if len(w) < 10:
        raise RuntimeError(
            f"only {len(w)} stars with A_V >= {min_av} pass the law cuts - "
            "not enough reddening on this sightline to measure R_V; "
            "the A(D) run is still produced (at the G23 default shape)")
    rv = w.rv.values
    rve = np.maximum(w.rv_err.values, 0.05)
    wm = 1.0 / rve ** 2
    rv_med = float(np.median(rv))
    rv_mad = float(1.4826 * np.median(np.abs(rv - rv_med)))
    ratios, ratios_mad = {}, {}
    for b in bands:
        q = (w[f"A_{b}"] / w["av"]).values     # A_band / A_V per star
        med = float(np.median(q))
        ratios[b] = med
        ratios_mad[b] = float(1.4826 * np.median(np.abs(q - med)))
    return dict(rv=rv_med, rv_mad=rv_mad,
                rv_mean=float((wm * rv).sum() / wm.sum()),
                rv_mean_err=float(1.0 / np.sqrt(wm.sum())),
                n_stars=int(len(w)), min_av=float(min_av),
                ratios_av=ratios, ratios_av_mad=ratios_mad,
                rv_16=float(np.percentile(rv, 16)), rv_84=float(np.percentile(rv, 84)),
                systematics="NIR zero-point freeze +/-0.05 on R_V; XP blue-end "
                            "systematics for faint heavily-reddened stars (see docs)")


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

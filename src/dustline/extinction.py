"""Extinction curves and band-integrated extinctions.

Wraps dust_extinction's G23 (default, R_V 2.3-5.6) and F99 parameter-averaged
curves, plus a simple power law; provides band-integrated A_band for a source
spectrum through any registered filter (the correct way to convert the
measured R_V into A_X for an arbitrary filter X).
"""

from __future__ import annotations

import numpy as np

RV_RANGE_G23 = (2.3, 5.6)


def curve(wave_A: np.ndarray, rv: float, law: str = "g23") -> np.ndarray:
    """A_lambda / A_V at wave_A for the dust_extinction G23 (default) or F99 curve."""
    import astropy.units as u
    from dust_extinction.parameter_averages import F99, G23

    x = 1.0 / (np.asarray(wave_A, float) * 1e-4)   # 1/micron
    mod = (G23 if law == "g23" else F99)(Rv=rv)
    lo, hi = mod.x_range
    xx = np.clip(x, lo + 1e-6, hi - 1e-6)
    return np.asarray(mod(xx / u.micron))


def powerlaw_curve(wave_A: np.ndarray, alpha: float, lam0_A: float = 7774.0) -> np.ndarray:
    """A_lambda / A(lam0) for a power law (lam/lam0)^-alpha."""
    return (np.asarray(wave_A, float) / lam0_A) ** (-alpha)


def band_extinction(photometry, flux_lam: np.ndarray, av: float, rv: float,
                    law: str = "g23") -> dict[str, float]:
    """Band-integrated A_band (mag) of a spectrum under A_V * curve(R_V).

    photometry: a models.Photometry instance (defines the bands and wavelength grid);
    flux_lam:   the unreddened spectrum on that grid (erg/s/cm2/A, any scale).
    """
    ext = curve(photometry.wave, rv, law)
    f0 = photometry.fnu(flux_lam)
    fe = photometry.fnu(flux_lam * 10.0 ** (-0.4 * av * ext))
    return {b: -2.5 * np.log10(fe[b] / f0[b]) for b in photometry.bands}

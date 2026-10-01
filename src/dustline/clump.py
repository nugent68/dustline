"""Red-clump bulge anchor: the extinction column and distance of the red-clump
overdensity in the field's near-infrared CMD, used to bridge the A(D) run
beyond the last Gaia-parallax bin on bulge sightlines.

Method (port of the research clump_anchor/vvv_clump work):
- clump peak in the (J-Ks, Ks) CMD (2D histogram mode + window medians);
- intrinsic clump colour/absolute magnitude from the NewEra clump model
  (4700 K, log g 2.5, [M/H] 0, log L/L_sun 1.72 -> M_Ks ~ -1.61) computed
  through the ACTUAL J/Ks curves in use (2MASS or VISTA);
- E(J-Ks) -> A_Ks with the MEASURED law ratios (A_J/A_V, A_Ks/A_V);
- D_RC from the dereddened Ks.

Calibration: the column is multiplied by ANCHOR_CAL (0.96; $DUSTLINE_ANCHOR_CAL or the
``cal`` argument override). Against the 85 microlensed bulge dwarfs of Bensby et al. (2017)
(benchmarks/bensby2017, Stage A/B) the J-Ks anchor converted to Cousins I is 1.050 x the
spectroscopic source A_I and 1.017 x the OGLE-III clump A_I of Nataf et al. (2013); the gap
between the two is the sources' median E(V-I) offset from their clump (-0.035 mag), which may be
dust or the microlensing colour calibration, so 0.96 sits between 1/1.050 and 1/1.017. The
uncalibrated values are kept (A_Ks_raw, AV_column_raw) and a 2 % calibration error is added.

Sanity gates: the anchor is only accepted if the window is well populated and
D_RC lands at a bulge-like distance (5-12 kpc); otherwise no bridge is applied
(the run then simply stops at the last parallax bin).
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

from . import models

CLUMP_MODEL = (4700.0, 2.5, 0.0)   # Teff, log g, [M/H]
LOGL_RC = 1.72                     # log L/L_sun that gives M_Ks ~ -1.61 (Nataf+13-like)
D_RC_RANGE = (5.0, 12.0)           # accepted bulge distances (kpc)
MIN_WINDOW = 60                    # stars in the clump window
ANCHOR_CAL = float(os.environ.get("DUSTLINE_ANCHOR_CAL", 0.96))   # see the module docstring
ANCHOR_CAL_ERR = 0.02              # relative calibration uncertainty, added in quadrature


def clump_intrinsic(j_band: str, ks_band: str) -> dict:
    """NewEra clump intrinsic (J-Ks)_0 and M_J/M_Ks through the given curves."""
    wave, flux, meta = models.load_cache()
    teff, logg, mh = CLUMP_MODEL
    k = np.where((meta[:, 0] == teff) & (meta[:, 1] == logg) & (meta[:, 2] == mh))[0]
    if not len(k):
        d2 = (meta[:, 0] - teff) ** 2 + 1e6 * (meta[:, 1] - logg) ** 2 + 1e8 * (meta[:, 2] - mh) ** 2
        k = [int(np.argmin(d2))]
    fl = flux[k[0]] * 100.0
    R = np.sqrt(10 ** LOGL_RC * 3.828e33
                / (4 * np.pi * 5.6704e-5 * teff ** 4)) / models.RSUN_CM
    ph = models.Photometry(wave, [j_band, ks_band])
    m0 = ph.mags(fl * models.RD2 * (R / 0.01) ** 2)      # absolute magnitudes (10 pc)
    return dict(JK0=m0[j_band] - m0[ks_band], MJ=m0[j_band], MKs=m0[ks_band])


def nir_bands(stars: pd.DataFrame, plan_nir: str, min_pairs: int = 500) -> tuple[str, str]:
    """(J, Ks) band names for the clump: the plan's survey when its photometry is actually there,
    else 2MASS (VVV comes from the DECaPS stellar-inference table and is missing where DECaPS
    does not reach, e.g. l > ~5.5 deg in the bulge)."""
    for pre in (("VISTA", "2MASS") if plan_nir == "vvv" else ("2MASS",)):
        J, K = stars.get(f"mag_{pre}_J"), stars.get(f"mag_{pre}_Ks")
        if J is not None and K is not None and int((np.isfinite(J) & np.isfinite(K)).sum()) >= min_pairs:
            return f"{pre}_J", f"{pre}_Ks"
    return "2MASS_J", "2MASS_Ks"


def find_clump(stars: pd.DataFrame, law: dict, j_band: str, ks_band: str,
               cal: float | None = None) -> dict | None:
    """Locate the red clump in (J-Ks, Ks) and derive the bulge column.

    stars: the field table with mag_<j_band>, mag_<ks_band> columns (ALL stars,
    not just the XP subset). law: the measured-law dict (ratios_av per A_V).
    cal: column calibration factor (default ANCHOR_CAL; 1.0 = the uncalibrated J-Ks column).
    Returns the anchor dict, or None when no credible clump is found.
    """
    from scipy.ndimage import gaussian_filter

    J = stars.get(f"mag_{j_band}")
    K = stars.get(f"mag_{ks_band}")
    if J is None or K is None:
        return None
    m = np.isfinite(J) & np.isfinite(K)
    jk, ks = (J - K)[m].values, K[m].values
    if m.sum() < 500:
        return None

    rJ = law["ratios_av"].get(j_band)
    rK = law["ratios_av"].get(ks_band)
    if not rJ or not rK:
        return None

    intr = clump_intrinsic(j_band, ks_band)
    # search the reddened-clump region: redward of the intrinsic colour, Ks 11-15.5, and on the
    # clump's reddening track at a bulge distance (each star dereddened as if it were a clump
    # giant); without the track cut, low-extinction fields (b ~ -5) peak on the bulge turnoff
    # at Ks ~ 15.3 instead of the clump at Ks ~ 13
    xr = (intr["JK0"] + 0.1, intr["JK0"] + 3.0)
    yr = (11.0, 15.5)
    mu_track = ks - (jk - intr["JK0"]) * rK / (rJ - rK) - intr["MKs"]
    mu_lo, mu_hi = (5 * np.log10(d * 100.0) for d in D_RC_RANGE)
    inwin = ((jk > xr[0]) & (jk < xr[1]) & (ks > yr[0]) & (ks < yr[1])
             & (mu_track > mu_lo) & (mu_track < mu_hi))
    if inwin.sum() < 200:
        return None
    H, xe, ye = np.histogram2d(jk[inwin], ks[inwin],
                               bins=[np.arange(*xr, 0.05), np.arange(*yr, 0.1)])
    Hs = gaussian_filter(H, 1.2)
    i0, j0 = np.unravel_index(np.argmax(Hs), Hs.shape)
    jk0 = 0.5 * (xe[i0] + xe[i0 + 1])
    ks0 = 0.5 * (ye[j0] + ye[j0 + 1])
    w = (np.abs(jk - jk0) < 0.2) & (np.abs(ks - ks0) < 0.5)
    if w.sum() < MIN_WINDOW:
        return None
    jk_rc = float(np.median(jk[w]))
    ks_rc = float(np.median(ks[w]))
    mad_jk = float(1.4826 * np.median(np.abs(jk[w] - jk_rc)))

    E_JK = jk_rc - intr["JK0"]
    if E_JK < 0.15:                     # essentially unreddened peak: a disk feature
        return None
    cal = ANCHOR_CAL if cal is None else float(cal)
    A_Ks_raw = E_JK * rK / (rJ - rK)    # measured-law A_Ks from the colour excess
    A_Ks = cal * A_Ks_raw               # calibrated (Bensby+17 microlensed dwarfs; module docstring)
    A_V_col = A_Ks / rK
    D_RC = 10 ** ((ks_rc - A_Ks - intr["MKs"]) / 5.0 - 2.0)
    if not (D_RC_RANGE[0] < D_RC < D_RC_RANGE[1]):
        return None
    # differential-extinction spread across the clump window, in A_V, (+) the calibration error
    sig_AV = cal * mad_jk * rK / (rJ - rK) / rK
    sig_AV = float(np.hypot(max(sig_AV, 0.1), ANCHOR_CAL_ERR * A_V_col if cal != 1.0 else 0.0))

    return dict(jk_rc=jk_rc, ks_rc=ks_rc, mad_jk=mad_jk, n_window=int(w.sum()),
                JK0=float(intr["JK0"]), MKs=float(intr["MKs"]),
                E_JK=float(E_JK), A_Ks=float(A_Ks), A_Ks_raw=float(A_Ks_raw),
                AV_column=float(A_V_col), AV_column_raw=float(A_Ks_raw / rK), cal=cal,
                AV_column_err=sig_AV,
                D_RC=float(D_RC), j_band=j_band, ks_band=ks_band,
                clump_model=dict(teff=CLUMP_MODEL[0], logg=CLUMP_MODEL[1],
                                 mh=CLUMP_MODEL[2], logL=LOGL_RC))

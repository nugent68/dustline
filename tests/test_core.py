"""Unit tests that need no model assets or network: filters, extinction curves,
footprint logic, synthetic photometry on toy spectra, ensemble statistics."""

import numpy as np
import pandas as pd
import pytest

from dustline import extinction, filters
from dustline.catalogs import footprint
from dustline.models import C_AA, Photometry


def test_builtin_filters_load():
    for name in ("DECam_i", "PS1_g", "2MASS_Ks", "I", "Cousins_I", "Gaia_G"):
        b = filters.get(name)
        assert b.wave_A.size > 5
        assert b.zp_jy > 100
        assert 3000 < b.lam_eff_A < 30000


def test_cousins_i_effective_wavelength():
    b = filters.get("I")
    assert 7600 < b.lam_eff_A < 8200   # Cousins I ~ 7900 A


def test_unknown_filter_raises():
    with pytest.raises(KeyError):
        filters.get("not_a_band")


def test_extinction_curve_monotonic_optical():
    wave = np.linspace(4000, 9000, 50)
    ext = extinction.curve(wave, rv=3.1)
    assert np.all(np.diff(ext) < 0)          # optical extinction falls with wavelength
    assert 1.0 < ext[0] < 2.0                # A(4000A)/A_V
    assert 0.3 < ext[-1] < 0.8               # A(9000A)/A_V


def test_extinction_rv_dependence():
    wave = np.array([4400.0])                # B band: lower R_V = steeper = more A_B/A_V
    a_lo = extinction.curve(wave, rv=2.5)[0]
    a_hi = extinction.curve(wave, rv=4.5)[0]
    assert a_lo > a_hi


def test_synthetic_photometry_flat_spectrum_ab():
    """A flat-f_nu spectrum must have AB colour 0 in every AB band."""
    wave = np.arange(3000.0, 11001.0, 2.0)
    fnu0 = 3631e-23 * 10 ** (-0.4 * 20.0)             # AB mag 20
    flam = fnu0 * C_AA / wave ** 2                    # erg/s/cm2/A
    ph = Photometry(wave, ["DECam_g", "DECam_r", "DECam_i", "PS1_z"])
    mags = ph.mags(flam)
    for b, m in mags.items():
        assert m == pytest.approx(20.0, abs=0.01), b


def test_band_extinction_ratios_sensible():
    """Band-integrated extinction of a flat spectrum: A_g > A_r > A_i > A_z."""
    wave = np.arange(3000.0, 11001.0, 2.0)
    flam = 1e-8 * (wave / 5000.0) ** -2
    ph = Photometry(wave, ["DECam_g", "DECam_r", "DECam_i", "DECam_z"])
    A = extinction.band_extinction(ph, flam, av=2.0, rv=3.1)
    assert A["DECam_g"] > A["DECam_r"] > A["DECam_i"] > A["DECam_z"]
    assert A["DECam_g"] / A["DECam_i"] == pytest.approx(1.9, abs=0.25)


def test_footprint_bulge():
    p = footprint.plan(267.86642, -33.13517)       # ob170095: l 357, b -3.2
    assert p.optical == "decaps"
    assert p.nir == "vvv"
    assert p.in_bulge_window


def test_footprint_ps1_default():
    p = footprint.plan(275.089125, -18.269389)     # ob240669: l 13.2 (no DECaPS/VVV)
    assert p.optical == "ps1"
    assert p.nir == "2mass"


def test_footprint_far_south_high_lat():
    p = footprint.plan(80.0, -70.0)                # LMC-ish: Dec < -30, high |b|
    assert p.optical == "decals"


def test_running_profile():
    from dustline.ensemble import D_FINE, running_profile
    rng = np.random.default_rng(0)
    D = rng.uniform(0.5, 5.0, 2000)
    A = 0.5 * D + rng.normal(0, 0.05, 2000)        # linear dust run
    prof = running_profile(D, A)
    k = np.argmin(np.abs(D_FINE - 2.0))
    assert prof[k, 0] == pytest.approx(1.0, abs=0.05)
    assert prof[k, 1] < prof[k, 0] < prof[k, 2]


def test_law_sample_cuts():
    from dustline.ensemble import law_sample
    r = pd.DataFrame(dict(
        parallax_over_error=[10, 10, 1, 10, 10],
        ruwe=[1.0, 2.0, 1.0, 1.0, 1.0],
        ipd_frac_multi_peak=[0, 0, 0, 50, 0],
        chi2_best=[100, 100, 100, 100, 100],
        n_xp=[343] * 5,
        n_phot=[6, 6, 6, 6, 2],
        av=[3.0, 3.0, 3.0, 3.0, 3.0],
        plx_used=[True, True, False, True, True],
        teff=[4500] * 5,
    ))
    sel = law_sample(r)
    assert sel.tolist() == [True, False, False, False, False]

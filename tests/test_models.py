"""Tests that need the model assets in the local cache (skipped when absent).

Validation anchors from the research work (dustline HANDOFF):
- the NewEra Sun model gives M_g 5.05, M_i 4.52, M_Ks 3.27 (absolute scale);
- G23 band ratios at R_V 3.15 for a 6000 K dwarf: A_g/A_i 1.885, A_r/A_i 1.335,
  A_z/A_i 0.767, A_J/A_i 0.442, A_Ks/A_i 0.179.
"""

import numpy as np
import pytest

from dustline import assets

pytestmark = pytest.mark.skipif(
    not (assets.cache_dir() / "newera_full_cache.npz").exists()
    or not (assets.cache_dir() / "mist_v1.2_basic.npz").exists(),
    reason="model assets not in the local cache",
)


def _sun_model():
    from dustline import models
    wave, flux, meta = models.load_cache()
    k = np.where((meta[:, 0] == 5800) & (meta[:, 1] == 4.5) & (meta[:, 2] == 0))[0]
    if not len(k):
        d2 = (meta[:, 0] - 5772) ** 2 + 1e6 * (meta[:, 1] - 4.44) ** 2 + 1e8 * meta[:, 2] ** 2
        k = [int(np.argmin(d2))]
    return wave, flux[k[0]], meta[k[0]]


def test_sun_absolute_magnitudes():
    from dustline import models
    wave, flux, meta = _sun_model()
    # absolute magnitude: R = 1 R_sun at D = 10 pc = 0.01 kpc
    scale = 100.0 * models.RD2 / 0.01 ** 2      # to erg/s/cm2/A at 10 pc
    ph = models.Photometry(wave, ["DECam_g", "DECam_i", "VISTA_Ks"])
    mags = ph.mags(flux * scale)
    assert mags["DECam_g"] == pytest.approx(5.05, abs=0.15)
    assert mags["DECam_i"] == pytest.approx(4.52, abs=0.15)
    assert mags["VISTA_Ks"] == pytest.approx(3.27, abs=0.20)


def test_band_ratios_at_rv315():
    """Reproduce the published ob170095 law ratios (6000 K source, G23 R_V 3.15)."""
    from dustline.ensemble import band_ratio_at_rv
    ref = {"DECam_g": 1.885, "DECam_r": 1.335, "DECam_z": 0.767,
           "VISTA_J": 0.442, "VISTA_Ks": 0.179}
    ri = band_ratio_at_rv("DECam_i", 3.15, teff=6000, logg=4.0)
    for band, expected in ref.items():
        r = band_ratio_at_rv(band, 3.15, teff=6000, logg=4.0)
        assert r / ri == pytest.approx(expected, rel=0.03), band


def test_cousins_i_ratio_close_to_decam_i():
    from dustline.ensemble import band_ratio_at_rv
    r_c = band_ratio_at_rv("I", 3.1)
    r_d = band_ratio_at_rv("DECam_i", 3.1)
    assert r_c == pytest.approx(r_d, rel=0.10)   # Cousins I (7900A) vs DECam i (7770A)


def test_radius_prior_dwarf_and_giant():
    from dustline import models
    meta = np.array([[5800.0, 4.5, 0.0],       # solar dwarf: R ~ 1 R_sun
                     [4800.0, 2.5, 0.0]])      # red giant: R ~ 10 R_sun
    rp = models.radius_prior(meta)
    assert rp[0, 0] == pytest.approx(0.0, abs=0.1)     # log R ~ 0
    assert 0.7 < rp[1, 0] < 1.4                        # log R ~ 1
    assert np.all(rp[:, 1] >= 0.08)


def test_xp_lsf_matrix_normalised():
    from dustline import models
    wave, flux, meta = models.load_cache()
    xp_wave = np.arange(336.0, 1021.0, 2.0)
    K = models.xp_lsf_matrix(wave, xp_wave)
    assert K.shape == (len(xp_wave), len(wave))
    np.testing.assert_allclose(K.sum(axis=1), 1.0, rtol=1e-4)

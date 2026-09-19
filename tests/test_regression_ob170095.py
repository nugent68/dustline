"""Regression of the ensemble statistics against the published ob170095 numbers,
using the research repo's per-star fit table as a fixture (skipped when absent).

Published (dustline research repo, 2026-09-17):
- law: R_V = 3.15 (median), MAD 0.23 (195 stars with A_V >= 2, plx S/N > 5);
- run: A_i 0.59 / 1.19 / 1.33 / 1.39 / 1.81 at 0.5-1 / 1-1.5 / 1.5-2 / 2-3 / 3-5 kpc.
"""

import os

import numpy as np
import pandas as pd
import pytest

FIXTURE = os.path.expanduser("~/claude/dustline/results/ob170095/xp_stars_phot.csv")
pytestmark = pytest.mark.skipif(not os.path.exists(FIXTURE),
                                reason="ob170095 research fixture not available")

# the research table uses short band names; A_i is DECam i
BANDS = ["g", "r", "i", "z", "Y", "J", "H", "Ks"]


@pytest.fixture(scope="module")
def stars():
    return pd.read_csv(FIXTURE)


def test_law_rv_regression(stars):
    from dustline.ensemble import law_sample, measure_law
    sel = law_sample(stars, min_snr=5.0, min_bands=5, min_av=2.0)
    law = measure_law(stars, BANDS, min_snr=5.0, min_bands=5, min_av=2.0)
    assert law["n_stars"] == pytest.approx(195, abs=20)
    assert law["rv"] == pytest.approx(3.15, abs=0.05)
    assert law["rv_mad"] == pytest.approx(0.23, abs=0.05)
    # band ratios vs A_V for the same stars: A_g/A_i via the two A_x/A_V ratios
    g_over_i = law["ratios_av"]["g"] / law["ratios_av"]["i"]
    assert g_over_i == pytest.approx(1.87, abs=0.08)
    assert sel.sum() == law["n_stars"]


def test_dust_run_regression(stars):
    from dustline.ensemble import good_sample
    s = stars[good_sample(stars, min_snr=5.0, min_bands=5)]
    ref = {(0.5, 1.0): 0.59, (1.0, 1.5): 1.19, (1.5, 2.0): 1.33,
           (2.0, 3.0): 1.39, (3.0, 5.0): 1.81}
    for (lo, hi), expected in ref.items():
        m = (s.D_kpc >= lo) & (s.D_kpc < hi)
        assert m.sum() >= 3, (lo, hi)
        med = float(np.median(s.A_i[m]))
        assert med == pytest.approx(expected, abs=0.05), (lo, hi)


def test_run_profile_shape(stars):
    from dustline.ensemble import good_sample, running_profile
    s = stars[good_sample(stars, min_snr=5.0, min_bands=5)]
    prof = running_profile(s.D_kpc.values, s.A_i.values)
    valid = np.isfinite(prof[:, 0])
    assert valid.sum() > 20
    # extinction increases outward overall
    first = prof[valid, 0][0]
    last = prof[valid, 0][-1]
    assert last > first

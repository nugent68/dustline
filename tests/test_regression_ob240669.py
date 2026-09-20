"""Regression of the ensemble products against the packaged OGLE-2024-BLG-0669
example (examples/ob240669), using the per-star fit table and the Gaia field
table (with its 2MASS join) of the cached sightline workspace as fixtures
(skipped when that workspace is not in the local cache; tools/run_ob240669.py
builds it).

Pinned (dustline v0.5.0: UV-IR model cache with the empirical dwarf-template
corrections, tiled PS1 fetch with saturation cut, 2MASS, 9'):
- law: R_V = 2.80 (median), MAD 0.46 (895 stars with A_V >= 2);
- clump: E(J-Ks) 0.52, A_V column 3.27 +/- 0.61 at D_RC 6.2 kpc;
- run: A_I 0.35 / 0.82 / 0.99 / 1.16 / 1.26 at 1 / 2 / 3 / 4 / 5 kpc, bridged to
  1.90 beyond the clump (A_I/A_V = 0.579).

The corresponding figures are in tests/figures/ob240669/ (law_rv.png,
extinction_run_I.png) for visual comparison when the numbers drift.
"""

import numpy as np
import pandas as pd
import pytest

from dustline.cache import Workspace

RA, DEC, RADIUS = 275.089125, -18.269389, 9.0
BANDS = ["2MASS_H", "2MASS_J", "2MASS_Ks", "PS1_g", "PS1_i", "PS1_r", "PS1_y", "PS1_z"]


def _workspace():
    from dustline import calib
    tag = calib.tag(calib.load())
    if not tag:
        return None
    ws = Workspace(RA, DEC, RADIUS, dict(optical="ps1", nir="2mass", user_phot=None,
                                         model_cache="uvir_cache", template_corr=tag))
    return ws if ws.has("xp_stars.csv") and ws.has("gaia.csv") else None


pytestmark = pytest.mark.skipif(_workspace() is None,
                                reason="ob240669 sightline workspace not in the local cache")


@pytest.fixture(scope="module")
def fit():
    return pd.read_csv(_workspace().path("xp_stars.csv"))


@pytest.fixture(scope="module")
def law(fit):
    from dustline.ensemble import measure_law
    return measure_law(fit, BANDS)


@pytest.fixture(scope="module")
def field():
    """The full Gaia field with the 2MASS join, named as bands.assemble does."""
    g = pd.read_csv(_workspace().path("gaia.csv"), usecols=["mag_J", "mag_Ks"])
    return g.rename(columns={"mag_J": "mag_2MASS_J", "mag_Ks": "mag_2MASS_Ks"})


@pytest.fixture(scope="module")
def anchor(field, law):
    from dustline.clump import find_clump
    return find_clump(field, law, "2MASS_J", "2MASS_Ks")


def test_law_rv_regression(law):
    assert law["n_stars"] == pytest.approx(895, abs=20)
    assert law["rv"] == pytest.approx(2.80, abs=0.05)
    assert law["rv_mad"] == pytest.approx(0.46, abs=0.05)
    assert law["ratios_av"]["PS1_i"] == pytest.approx(0.622, abs=0.02)
    assert law["ratios_av"]["2MASS_Ks"] == pytest.approx(0.096, abs=0.01)


def test_clump_anchor_regression(anchor):
    assert anchor is not None
    assert anchor["E_JK"] == pytest.approx(0.52, abs=0.03)
    assert anchor["AV_column"] == pytest.approx(3.27, abs=0.15)
    assert anchor["D_RC"] == pytest.approx(6.1, abs=0.3)


def test_dust_run_regression(fit, law, anchor):
    from dustline.ensemble import band_ratio_at_rv, bridge_to_clump, dust_run
    run = bridge_to_clump(dust_run(fit), anchor)
    ratio = band_ratio_at_rv("I", law["rv"])
    assert ratio == pytest.approx(0.579, abs=0.01)
    ref = {1.0: 0.35, 2.0: 0.82, 3.0: 0.99, 4.0: 1.16, 5.0: 1.26, 8.0: 1.90, 12.0: 1.90}
    for D, expected in ref.items():
        r = run[np.isclose(run.D_kpc, D)].iloc[0]
        assert r.AV_med * ratio == pytest.approx(expected, abs=0.03), D
        assert bool(r.bridged) == (D > 6.3), D
    assert run.D_kpc.iloc[-1] == pytest.approx(12.0, abs=0.01)

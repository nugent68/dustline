"""Regression of the ensemble products against the packaged OGLE-2017-BLG-0095
example (examples/ob170095), using the per-star fit table and the band-assembled
field of the cached sightline workspace as fixtures (skipped when that workspace
is not in the local cache; `dustline run 267.86642 -33.13517 --radius 5 --plx-inflate 1.7`
builds it in ~2 h, the XP fetch dominating).

Pinned (dustline v0.7.1: v0.7.0 stack + parallax errors x1.7 and parallax x
photometric distance posteriors; DECaPS grizY + VVV JHKs, 5'):
- law: R_V = 3.16 (median), MAD 0.23 (326 stars with A_V >= 2); raw 3.07;
- clump: E(J-Ks) 0.44, A_V column 2.71 +/- 0.44 at D_RC 8.4 kpc;
- in-field parallax-error underestimate x1.65 (robust) from the clump window;
- run: A_I 0.67 / 1.23 / 1.32 / 1.40 / 1.55 at 1 / 2 / 3 / 4 / 5 kpc, bridged to
  1.66 beyond the clump (A_I/A_V = 0.609).
"""

import numpy as np
import pandas as pd
import pytest

from dustline.cache import Workspace

RA, DEC, RADIUS = 267.86642, -33.13517, 5.0
BANDS = ["DECam_Y", "DECam_g", "DECam_i", "DECam_r", "DECam_z", "VISTA_H", "VISTA_J", "VISTA_Ks"]


def _workspace():
    from dustline import calib
    tag = calib.tag(calib.load())
    if not tag:
        return None
    ws = Workspace(RA, DEC, RADIUS, dict(optical="decaps", nir="vvv", user_phot=None,
                                         plx_inflate=1.7, model_cache="uvir_cache", template_corr=tag))
    return ws if ws.has("xp_stars.csv") and ws.has("gaia.csv") and ws.has("vvv_gaia.csv") else None


pytestmark = pytest.mark.skipif(_workspace() is None,
                                reason="ob170095 sightline workspace not in the local cache")


@pytest.fixture(scope="module")
def fit():
    return pd.read_csv(_workspace().path("xp_stars.csv"))


@pytest.fixture(scope="module")
def law(fit):
    from dustline.ensemble import measure_law
    return measure_law(fit, BANDS)


@pytest.fixture(scope="module")
def field():
    """The Gaia field with the VVV JHKs, named as bands.assemble does."""
    ws = _workspace()
    g = pd.read_csv(ws.path("gaia.csv"), low_memory=False)
    v = pd.read_csv(ws.path("vvv_gaia.csv"))
    return g.merge(v.drop(columns=["vvv_sep"]), on="source_id", how="left")


@pytest.fixture(scope="module")
def anchor(field, law):
    from dustline.clump import find_clump
    return find_clump(field, law, "VISTA_J", "VISTA_Ks")


def test_law_rv_regression(law):
    assert law["n_stars"] == pytest.approx(326, abs=15)
    assert law["rv"] == pytest.approx(3.16, abs=0.05)
    assert law["rv_raw"] == pytest.approx(3.07, abs=0.05)
    assert law["rv_mad"] == pytest.approx(0.23, abs=0.04)
    assert law["ratios_av"]["DECam_i"] == pytest.approx(0.615, abs=0.02)
    assert law["ratios_av"]["VISTA_Ks"] == pytest.approx(0.111, abs=0.01)


def test_clump_anchor_regression(anchor):
    assert anchor is not None
    assert anchor["E_JK"] == pytest.approx(0.44, abs=0.03)
    assert anchor["AV_column"] == pytest.approx(2.71, abs=0.15)
    assert anchor["D_RC"] == pytest.approx(8.4, abs=0.3)


def test_parallax_inflation_regression(field, anchor):
    from dustline.ensemble import parallax_inflation
    infl = parallax_inflation(field, anchor)
    assert infl is not None and infl["n"] > 1000
    assert infl["r_robust"] == pytest.approx(1.65, abs=0.1)


def test_dust_run_regression(fit, law, anchor):
    from dustline.ensemble import band_ratio_at_rv, bridge_to_clump, dust_run
    assert (fit.plx_inflate == 1.7).all()
    run = bridge_to_clump(dust_run(fit), anchor)
    ratio = band_ratio_at_rv("I", law["rv"])
    assert ratio == pytest.approx(0.609, abs=0.01)
    ref = {1.0: 0.67, 2.0: 1.23, 3.0: 1.32, 4.0: 1.40, 5.0: 1.55, 8.0: 1.66, 12.0: 1.66}
    for D, expected in ref.items():
        r = run[np.isclose(run.D_kpc, D)].iloc[0]
        assert r.AV_med * ratio == pytest.approx(expected, abs=0.04), D
        assert bool(r.bridged) == (D > 5.2), D

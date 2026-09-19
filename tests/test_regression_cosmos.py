"""Regression of column mode + DESI priors against the packaged COSMOS example
(examples/cosmos), using the cached sightline workspace as the fixture (skipped
when absent; `dustline run 150.12 2.21 --radius 30` builds it in ~25 min).

Pinned (dustline v0.4.0: UV-IR cache, PS1 + 2MASS + WISE + GALEX, 30', DESI DR1
priors, R_V fixed 3.05):
- 592 XP stars, 422 with DESI priors; F/G column A_V 0.072 +/- 0.009 (N 149);
  K/M stars 0.19; T_eff(DESI) - T_eff(fit) ~ 33 K; WISE W1 offset ~ 0.
"""

import numpy as np
import pandas as pd
import pytest

from dustline.cache import Workspace

RA, DEC, RADIUS = 150.12, 2.21, 30.0


def _workspace():
    ws = Workspace(RA, DEC, RADIUS, dict(optical="ps1", nir="2mass", user_phot=None,
                                         spectro="desi", mode="column", uv="galex", mir="wise",
                                         model_cache="uvir_cache"))
    return ws if ws.has("xp_stars.csv") else None


pytestmark = pytest.mark.skipif(_workspace() is None,
                                reason="COSMOS sightline workspace not in the local cache")


@pytest.fixture(scope="module")
def fit():
    return pd.read_csv(_workspace().path("xp_stars.csv"))


def test_priors_present_and_consistent(fit):
    assert len(fit) == 592
    assert fit.spec_prior.sum() == 422
    sp = fit[fit.spec_prior.astype(bool)]
    assert abs(np.median(sp.teff_spec - sp.teff)) < 60
    assert abs(np.median(sp.logg_spec - sp.logg)) < 0.15
    assert set(np.unique(fit.mh_best)) <= {-2.0, -1.5, -1.0, -0.5, 0.0, 0.5}
    free = fit[~fit.spec_prior.astype(bool)]
    assert free.mh_best.between(-0.5, 0.5).all()              # no prior -> solar range only
    assert np.isfinite(fit.dm_WISE_W1).sum() > 400 and np.isfinite(fit.dm_GALEX_NUV).sum() > 100
    assert np.allclose(fit.rv, 3.05) and (fit.rv_err == 0).all()   # R_V fixed in column mode


def test_column_regression(fit):
    from dustline.ensemble import foreground_column
    c = foreground_column(fit)
    assert c["clean"]["n"] == pytest.approx(149, abs=10)
    assert c["clean"]["av"] == pytest.approx(0.072, abs=0.02)
    assert c["clean"]["av_mad"] < 0.08
    assert c["teff_cool"]["av"] > 0.15                        # the cool-star systematic
    assert c["clean_with_spec_prior"]["av_mad"] < 0.05
    assert c["clean_with_nuv"]["n"] > 80
    assert c["map"]["ncell"] == 4

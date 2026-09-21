"""Regression of column mode + DESI priors against the packaged COSMOS example
(examples/cosmos), using the cached sightline workspace as the fixture (skipped
when absent; `dustline run 150.12 2.21 --radius 30` builds it in ~25 min).

Pinned (dustline v0.5.0/v0.6.0: UV-IR cache, empirical template corrections, T_eff locked
to the BP-RP locus, PS1 + 2MASS + WISE + GALEX, 30', DESI log g/[Fe/H], R_V fixed 3.05):
- 592 XP stars; F/G column A_V 0.029 +/- 0.004 (MAD 0.052, N 132); K/M stars 0.039 (v0.7.0)
  (the two agree - the K-dwarf systematic is gone); G-trend 0.01/0.03/0.03/0.05.
"""

import numpy as np
import pandas as pd
import pytest

from dustline.cache import Workspace

RA, DEC, RADIUS = 150.12, 2.21, 30.0


def _workspace():
    from dustline import calib
    tag = calib.tag(calib.load())
    if not tag:
        return None
    ws = Workspace(RA, DEC, RADIUS, dict(optical="ps1", nir="2mass", user_phot=None,
                                         spectro="desi", mode="column", teff_lock="colour",
                                         uv="galex", mir="wise", model_cache="uvir_cache",
                                         template_corr=tag))
    return ws if ws.has("xp_stars.csv") else None


pytestmark = pytest.mark.skipif(_workspace() is None,
                                reason="COSMOS sightline workspace not in the local cache")


@pytest.fixture(scope="module")
def fit():
    return pd.read_csv(_workspace().path("xp_stars.csv"))


def test_lock_and_priors(fit):
    assert len(fit) == 592
    assert fit.spec_prior.sum() >= 585                        # every star with BP-RP is locked
    sp = fit[fit.spec_prior.astype(bool)]
    assert abs(np.median(sp.teff_spec - sp.teff)) < 15          # locked to the nearest node
    assert set(np.unique(fit.mh_best)) <= {-2.0, -1.5, -1.0, -0.5, 0.0, 0.5}
    nofeh = fit[~np.isfinite(fit.feh_spec)]
    assert nofeh.mh_best.between(-0.5, 0.5).all()             # no [Fe/H] prior -> solar range
    assert np.isfinite(fit.dm_WISE_W1).sum() > 400 and np.isfinite(fit.dm_GALEX_NUV).sum() > 100
    assert np.allclose(fit.rv, 3.05) and (fit.rv_err == 0).all()   # R_V fixed in column mode


def test_column_regression(fit):
    from dustline.ensemble import foreground_column
    c = foreground_column(fit)
    assert c["clean"]["n"] == pytest.approx(132, abs=10)
    assert c["clean"]["av"] == pytest.approx(0.029, abs=0.02)
    assert c["clean"]["av_mad"] < 0.07
    assert abs(c["teff_cool"]["av"] - c["clean"]["av"]) < 0.02   # cool and hot stars agree
    assert c["by_gmag"]["16.5-18"]["av"] < 0.08                  # faint-end excess < 0.05
    assert c["clean_with_nuv"]["n"] > 80
    assert c["map"]["ncell"] == 4

"""Tests of the user-photometry adapter and band assembly (no network)."""

import numpy as np
import pandas as pd
import pytest

from dustline.userphot import UserPhotometry


@pytest.fixture
def phot_csv(tmp_path):
    df = pd.DataFrame(dict(
        ra=[10.0001, 10.0002, 10.5],
        dec=[-30.0001, -30.0002, -30.5],
        mag_g=[15.0, 16.0, 17.0],
        magerr_g=[0.01, 0.02, 0.03],
        mag_K=[13.0, 14.0, 15.0],
        magerr_K=[0.02, 0.02, 0.02],
    ))
    p = tmp_path / "phot.csv"
    df.to_csv(p, index=False)
    return p


def test_userphot_validates_bands(phot_csv):
    with pytest.raises(KeyError):
        UserPhotometry(phot_csv, bands={"g": "NoSuchBand"})


def test_userphot_missing_column(tmp_path):
    p = tmp_path / "bad.csv"
    pd.DataFrame(dict(ra=[1.0], dec=[2.0], mag_g=[15.0])).to_csv(p, index=False)
    with pytest.raises(ValueError, match="magerr_g"):
        UserPhotometry(p, bands={"g": "DECam_g"})


def test_userphot_match(phot_csv):
    up = UserPhotometry(phot_csv, bands={"g": "DECam_g", "K": "VISTA_Ks"})
    gaia = pd.DataFrame(dict(
        source_id=[1, 2, 3],
        ra=[10.0001, 10.0002, 11.0],
        dec=[-30.0001, -30.0002, -31.0],
    ))
    out = up.matched_to_gaia(gaia, 10.0, -30.0)
    assert set(out.columns) == {"source_id", "mag_DECam_g", "magerr_DECam_g",
                                "mag_VISTA_Ks", "magerr_VISTA_Ks"}
    assert len(out) == 2                      # the third star is 1 deg away
    row = out[out.source_id == 1].iloc[0]
    assert row.mag_DECam_g == pytest.approx(15.0)
    assert row.mag_VISTA_Ks == pytest.approx(13.0)


def test_workspace_caching(tmp_path, monkeypatch):
    monkeypatch.setenv("DUSTLINE_CACHE_DIR", str(tmp_path))
    from dustline.cache import Workspace
    ws1 = Workspace(10.0, -30.0, 5.0, config=dict(optical="ps1"))
    ws2 = Workspace(10.0, -30.0, 5.0, config=dict(optical="ps1"))
    ws3 = Workspace(10.0, -30.0, 5.0, config=dict(optical="decaps"))
    assert ws1.dir == ws2.dir                 # same position + config = same cache
    assert ws1.dir != ws3.dir                 # different config = different cache
    np.savetxt(ws1.path("x.txt"), [1.0])
    assert ws2.has("x.txt")

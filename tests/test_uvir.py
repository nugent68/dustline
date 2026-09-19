"""GALEX / AllWISE catalog shaping, the UV-IR plan slots and the model-cache
coverage guard (synthetic inputs; the live queries are 'network' tests)."""

import numpy as np
import pandas as pd
import pytest

from dustline import assets, filters
from dustline.catalogs import footprint, galex, wise

HAS_ASSETS = (assets.cache_dir() / "newera_full_cache.npz").exists()


def _gaia(n=4):
    return pd.DataFrame(dict(source_id=np.arange(1, n + 1), ra=150.0 + np.arange(n) * 0.01,
                             dec=np.full(n, 2.0), bp_rp=[0.5, 0.8, 1.3, 0.4]))


def test_footprint_uv_mir_slots():
    p = footprint.plan(150.12, 2.21)
    assert p.uv == "galex" and p.mir == "wise"
    q = footprint.plan(275.089125, -18.269389)          # b -1.6: no GALEX AIS, WISE confused
    assert q.uv == "none" and q.mir == "none"


def test_galex_filters_registered():
    for b in ("GALEX_FUV", "GALEX_NUV"):
        band = filters.get(b)
        assert band.system == "AB" and band.zp_jy == 3631.0
    assert filters.get("GALEX_NUV").wave_A.min() > 1600


def test_galex_shape_colour_and_flag_cuts():
    g = _gaia()
    d = pd.DataFrame({
        "RAJ2000": g.ra.values, "DEJ2000": g.dec.values,
        "NUVmag": [19.0, 19.5, 18.0, 20.0], "e_NUVmag": [0.05, 0.05, 0.05, 0.5],
        "FUVmag": [20.0, 21.0, 19.0, 21.0], "e_FUVmag": [0.1, 0.1, 0.1, 0.1],
        "Nafl": [0, 0, 0, 0], "Nexf": [0, 0, 0, 0], "Fafl": [0, 1, 0, 0], "Fexf": [0, 0, 0, 0],
        "E(B-V)": [0.018] * 4})
    r = galex.shape(d, g, 150.015, 2.0).set_index("source_id")
    assert np.isfinite(r.mag_GALEX_NUV[1]) and np.isfinite(r.mag_GALEX_NUV[2])
    assert np.isnan(r.mag_GALEX_NUV[3])                    # BP-RP 1.3: too cool for NUV
    assert np.isnan(r.mag_GALEX_NUV[4])                    # error 0.5 > MAX_ERR
    assert np.isfinite(r.mag_GALEX_FUV[1]) and np.isnan(r.mag_GALEX_FUV[2])   # artefact flag
    assert np.isnan(r.mag_GALEX_FUV[3])                    # BP-RP 1.3 > FUV colour cut
    assert r.magerr_GALEX_NUV[1] == pytest.approx(np.hypot(0.05, galex.SYS_FLOOR))


def test_wise_shape_flags_and_saturation():
    g = _gaia()
    d = pd.DataFrame({
        "ra": g.ra.values, "dec": g.dec.values,
        "w1mpro": [12.0, 7.5, 13.0, 14.0], "w1sigmpro": [0.03, 0.03, 0.03, 0.03],
        "w2mpro": [12.1, 7.6, 13.1, 14.1], "w2sigmpro": [0.03, 0.03, 0.03, 0.4],
        "cc_flags": ["0000", "0000", "d000", "0000"], "ext_flg": [0, 0, 0, 0]})
    r = wise.shape(d, g, 150.015, 2.0).set_index("source_id")
    assert np.isfinite(r.mag_WISE_W1[1]) and np.isfinite(r.mag_WISE_W2[1])
    assert np.isnan(r.mag_WISE_W1[2])                      # saturated
    assert np.isnan(r.mag_WISE_W1[3]) and np.isfinite(r.mag_WISE_W2[3])   # W1 cc flag only
    assert np.isnan(r.mag_WISE_W2[4])                      # error cut


@pytest.mark.skipif(not HAS_ASSETS, reason="model assets not in the local cache")
def test_cache_coverage_guard(monkeypatch):
    from dustline import models
    monkeypatch.delenv("DUSTLINE_MODEL_CACHE", raising=False)
    wave = models.load_cache()[0]
    assert models.cache_covers("PS1_g", wave) and models.cache_covers("2MASS_Ks", wave)
    if wave.min() > 2000:                                   # the 2500-25000 A cache
        assert not models.cache_covers("GALEX_NUV", wave)
        assert not models.cache_covers("WISE_W1", wave)
        with pytest.raises(ValueError, match="not covered"):
            models.Photometry(wave, ["GALEX_NUV"])
    else:                                                   # the UV-IR cache
        assert models.cache_covers("GALEX_FUV", wave) and models.cache_covers("WISE_W2", wave)


@pytest.mark.network
def test_galex_and_wise_fetch_cosmos(tmp_path, monkeypatch):
    from dustline.cache import Workspace
    monkeypatch.setenv("DUSTLINE_CACHE_DIR", str(tmp_path))
    ws = Workspace(150.12, 2.21, 6.0, dict(uv="galex", mir="wise"))
    g = pd.DataFrame(dict(source_id=[], ra=[], dec=[], bp_rp=[]))
    raw = galex.box_query("II/335/galex_ais", galex._COLS, 150.12, 2.21, 0.1)
    assert len(raw) > 20 and np.isfinite(raw.NUVmag).sum() > 20
    assert galex.fetch(ws, g).shape[0] == 0 and ws.has("galex_sfd.json")
    raw = wise.box_query("allwise.source", wise._COLS, 150.12, 2.21, 0.1)
    assert len(raw) > 100 and np.isfinite(raw.w1mpro).sum() > 100
    assert wise.fetch(ws, g).shape[0] == 0 and ws.has("wise_gaia.csv")

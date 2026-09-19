"""Tests of the red-clump bulge bridge (synthetic CMD; needs model assets only
for clump_intrinsic, skipped when absent)."""

import numpy as np
import pandas as pd
import pytest

from dustline import assets

HAS_ASSETS = (assets.cache_dir() / "newera_full_cache.npz").exists()


def _fake_law():
    # G23-ish NIR ratios per A_V
    return dict(rv=2.7, rv_mad=0.3, ratios_av={"2MASS_J": 0.25, "2MASS_Ks": 0.09})


def _fake_field(n=4000, e_jk=0.75, jk0=0.62, mks=-1.61, d_rc=8.2, seed=1):
    """A bulge CMD: a reddened clump blob + disk foreground + noise."""
    rng = np.random.default_rng(seed)
    n_rc = n // 2
    jk_rc = jk0 + e_jk + rng.normal(0, 0.06, n_rc)
    ks_rc = (mks + 5 * np.log10(d_rc * 100)
             + 0.09 / (0.25 - 0.09) * e_jk         # A_Ks from the fake law
             + rng.normal(0, 0.22, n_rc))
    n_fg = n - n_rc
    jk_fg = rng.uniform(0.2, 2.4, n_fg)
    ks_fg = rng.uniform(11.0, 16.5, n_fg)
    return pd.DataFrame({
        "mag_2MASS_J": np.concatenate([jk_rc + ks_rc, jk_fg + ks_fg]),
        "mag_2MASS_Ks": np.concatenate([ks_rc, ks_fg]),
    })


@pytest.mark.skipif(not HAS_ASSETS, reason="model assets not in the local cache")
def test_find_clump_recovers_column():
    from dustline.clump import find_clump
    law = _fake_law()
    stars = _fake_field()
    anchor = find_clump(stars, law, "2MASS_J", "2MASS_Ks")
    assert anchor is not None
    assert anchor["E_JK"] == pytest.approx(0.75, abs=0.08)
    # A_Ks = E_JK * rK/(rJ-rK) = 0.75 * 0.5625; A_V = A_Ks / 0.09
    assert anchor["AV_column"] == pytest.approx(0.75 * 0.09 / 0.16 / 0.09, rel=0.12)
    assert anchor["D_RC"] == pytest.approx(8.2, abs=0.8)


@pytest.mark.skipif(not HAS_ASSETS, reason="model assets not in the local cache")
def test_find_clump_rejects_unreddened():
    from dustline.clump import find_clump
    stars = _fake_field(e_jk=0.05)      # nearly no reddening -> E(J-Ks) gate
    assert find_clump(stars, _fake_law(), "2MASS_J", "2MASS_Ks") is None


def test_find_clump_needs_bands():
    from dustline.clump import find_clump
    assert find_clump(pd.DataFrame({"x": [1.0]}), _fake_law(),
                      "2MASS_J", "2MASS_Ks") is None


def test_bridge_to_clump():
    from dustline.ensemble import bridge_to_clump
    run = pd.DataFrame(dict(D_kpc=[1.0, 2.0, 3.0], AV_med=[1.0, 1.5, 2.0],
                            AV_16=[0.8, 1.3, 1.8], AV_84=[1.2, 1.7, 2.2],
                            n_stars=[50, 40, 20], bridged=[False] * 3))
    anchor = dict(D_RC=8.0, AV_column=3.0, AV_column_err=0.3)
    out = bridge_to_clump(run, anchor)
    assert out.bridged.sum() > 0
    assert (out[out.bridged].n_stars == 0).all()
    # at D_RC the bridge hits the column; beyond it stays flat
    at_rc = out[np.isclose(out.D_kpc, 8.0)].iloc[0]
    assert at_rc.AV_med == pytest.approx(3.0, abs=0.02)
    end = out.iloc[-1]
    assert end.AV_med == pytest.approx(3.0, abs=0.02)
    assert end.AV_16 == pytest.approx(2.7, abs=0.05)
    assert end.AV_84 == pytest.approx(3.3, abs=0.05)
    # continuous at the junction
    j = out[np.isclose(out.D_kpc, 3.1)].iloc[0]
    assert abs(j.AV_med - 2.0) < 0.05


def test_bridge_anchor_inside_run_extends_flat():
    """When D_RC falls inside the parallax range, the run is still extended
    (short ramp then flat at the clump column) rather than left truncated."""
    from dustline.ensemble import bridge_to_clump
    run = pd.DataFrame(dict(D_kpc=[8.0, 9.0], AV_med=[2.0, 2.1], AV_16=[1.8, 1.9],
                            AV_84=[2.2, 2.3], n_stars=[20, 15], bridged=[False] * 2))
    out = bridge_to_clump(run, dict(D_RC=7.0, AV_column=3.0, AV_column_err=0.3))
    assert out.bridged.sum() > 0
    end = out.iloc[-1]
    assert end.D_kpc == pytest.approx(12.0, abs=0.11)
    assert end.AV_med == pytest.approx(3.0, abs=0.02)
    # junction ramp is short (0.5 kpc past the last bin)
    at_10 = out[np.isclose(out.D_kpc, 10.0)].iloc[0]
    assert at_10.AV_med == pytest.approx(3.0, abs=0.02)

"""Column mode and the DESI spectroscopic template priors (synthetic inputs;
the live Data Lab query is a 'network' test)."""

import numpy as np
import pandas as pd
import pytest

from dustline import ensemble
from dustline.catalogs import desi, footprint


def test_footprint_cosmos_column_mode():
    p = footprint.plan(150.12, 2.21)                   # COSMOS: l 237, b +42
    assert p.spectro == "desi" and p.mode == "column"
    assert p.optical == "ps1" and p.nir == "2mass" and not p.in_bulge_window


def test_footprint_bulge_unchanged():
    p = footprint.plan(275.089125, -18.269389)        # ob240669: b -1.6
    assert p.spectro == "none" and p.mode == "law"


def test_desi_shape_dedups_and_floors():
    raw = pd.DataFrame(dict(
        source_id=[1, 1, 2, 3, np.nan], target_ra=[150.0] * 5, target_dec=[2.0] * 5,
        teff=[5000.0, 5100.0, 6000.0, np.nan, 4500.0], teff_err=[10.0, 10.0, 20.0, 10.0, 10.0],
        logg=[4.5, 4.4, 4.0, 4.0, 4.0], logg_err=[0.01, 0.01, 0.05, 0.05, 0.05],
        feh=[-0.3, -0.2, -1.2, 0.0, 0.0], feh_err=[0.01, 0.01, 0.05, 0.05, 0.05],
        alphafe=[0.1] * 5, snr_med=[20.0, 50.0, 30.0, 10.0, 10.0],
        survey=["main"] * 5, program=["bright"] * 5))
    out = desi.shape(raw)
    assert list(out.source_id) == [1, 2]               # dup -> highest S/N; NaN teff / id dropped
    assert out.teff_spec.iloc[0] == 5100.0
    assert out.teff_spec_err.iloc[0] == pytest.approx(np.hypot(10.0, desi.TEFF_SYS))
    assert out.feh_spec_err.iloc[1] == pytest.approx(np.hypot(0.05, desi.FEH_SYS))


def test_spec_prior_chi2_pulls_and_clamps():
    from dustline.fit import spec_prior_chi2
    meta = np.array([[5000.0, 4.5, -0.5], [5000.0, 4.5, 0.0], [6000.0, 4.0, 0.0],
                     [6000.0, 4.0, 0.5]])
    spec = pd.DataFrame(dict(
        teff_spec=[6000.0, np.nan, 5000.0], teff_spec_err=[60.0, np.nan, 60.0],
        logg_spec=[4.0, np.nan, 4.5], logg_spec_err=[0.1, np.nan, 0.1],
        feh_spec=[0.0, np.nan, -1.5], feh_spec_err=[0.1, np.nan, 0.1]))
    chi, has = spec_prior_chi2(meta, spec)
    assert chi.shape == (4, 3) and list(has) == [True, False, True]
    assert np.argmin(chi[:, 0]) == 2 and chi[2, 0] == 0.0      # star 0 -> its model
    assert np.all(chi[:, 1] == 0.0)                             # no prior -> no penalty
    # star 2 has [Fe/H] -1.5: clamped to the grid edge -0.5, no penalty there
    assert np.argmin(chi[:, 2]) == 0 and chi[0, 2] == 0.0
    assert chi[1, 2] == pytest.approx(25.0)                    # ([M/H] 0 vs -0.5)/0.1


def test_spec_prior_pulls_posterior():
    """A star whose data prefer model 0 but whose prior points to model 1."""
    from dustline.fit import _summarize
    grid = dict(meta=np.array([[5000.0, 4.5, 0.0], [6000.0, 4.0, 0.0]]),
                av=np.array([0.0, 1.0]), rv=np.array([3.1]),
                A_band=np.zeros((1, 2, 2, 1), np.float32))
    chi2 = np.zeros((1, 2, 2, 1))
    chi2[0, 0, 1, 0] = 4.0                                      # data mildly prefer model 0
    C = np.ones_like(chi2)
    df0 = _summarize(grid, chi2, np.zeros_like(chi2), C, ["X"])
    prior = np.zeros_like(chi2)
    prior[0, :, 0, 0] = 100.0                                   # strong prior against model 0
    df1 = _summarize(grid, chi2, prior, C, ["X"])
    assert df0.teff_best.iloc[0] == 5000.0 and df1.teff_best.iloc[0] == 6000.0
    assert df1.teff.iloc[0] > df0.teff.iloc[0]


def _fake_fit(n=400, av0=0.06, seed=3):
    rng = np.random.default_rng(seed)
    D = rng.uniform(0.2, 3.0, n)
    av = np.where(D > 0.3, av0, 0.0) + rng.normal(0, 0.05, n)
    return pd.DataFrame(dict(
        ra=150.12 + rng.uniform(-0.5, 0.5, n), dec=2.21 + rng.uniform(-0.5, 0.5, n),
        D_kpc=D, D_err=0.05 * D, av=av, av_err=np.full(n, 0.05),
        teff=rng.uniform(4000, 7000, n), parallax_over_error=np.full(n, 20.0),
        phot_g_mean_mag=rng.uniform(12, 17.5, n),
        plx_used=np.ones(n, bool), ruwe=np.ones(n), ipd_frac_multi_peak=np.zeros(n),
        chi2_best=np.full(n, 100.0), n_xp=np.full(n, 343), n_phot=np.full(n, 8),
        spec_prior=rng.uniform(size=n) < 0.7, mh_clamped=np.zeros(n, bool),
        teff_spec=rng.uniform(4000, 7000, n)))


def test_foreground_column_recovers_screen():
    fit = _fake_fit()
    c = ensemble.foreground_column(fit)
    assert c["n"] > 250
    assert c["av"] == pytest.approx(0.06, abs=0.01)
    assert 0 < c["av_err"] < 0.01
    assert c["with_spec_prior"]["n"] + c["without_spec_prior"]["n"] == c["n"]
    assert c["clean"]["n"] + c["teff_cool"]["n"] == c["n"]
    assert c["clean"]["av"] == pytest.approx(0.06, abs=0.012)
    assert c["map"]["ncell"] == 4 and len(c["map"]["cells"]) == 16
    assert sum(cell["n"] for cell in c["map"]["cells"]) == c["clean"]["n"]


def test_foreground_column_empty():
    fit = _fake_fit(n=20)
    c = ensemble.foreground_column(fit[fit.D_kpc < 0.3])
    assert c["n"] == 0 and np.isnan(c["av"])


@pytest.mark.network
def test_desi_fetch_cosmos(tmp_path, monkeypatch):
    from dustline.cache import Workspace
    monkeypatch.setenv("DUSTLINE_CACHE_DIR", str(tmp_path))
    ws = Workspace(150.12, 2.21, 6.0, dict(spectro="desi"))
    gaia = pd.DataFrame(dict(source_id=[]))
    raw = desi.box_query("desi_dr1.mws", desi._COLS, 150.12, 2.21, 0.15,
                         ra_col="target_ra", dec_col="target_dec", extra=desi._EXTRA)
    assert len(raw) > 40 and np.isfinite(raw.teff).sum() > 40
    out = desi.fetch(ws, gaia)                     # no Gaia sources -> empty but no crash
    assert len(out) == 0 and ws.has("desi_gaia.csv")

"""Empirical template corrections: aggregation and application on synthetic input."""

import numpy as np
import pandas as pd

from dustline import calib


def _corr():
    wave = np.arange(336.0, 1021.0, 2.0)
    n = 60
    rng = np.random.default_rng(0)
    fit = pd.DataFrame(dict(teff_spec=np.r_[np.full(n // 2, 4600.0), np.full(n // 2, 5600.0)],
                            feh_spec=np.full(n, 0.0),
                            dm_PS1_g=np.r_[np.full(n // 2, 0.10), np.full(n // 2, 0.0)] + rng.normal(0, 0.01, n),
                            dm_2MASS_Ks=np.zeros(n)))
    slope = np.where(fit.teff_spec.values[:, None] < 5000, 1 - 0.1 * (700 - wave[None, :]) / 350, 1.0)
    ratios = slope * (1 + rng.normal(0, 0.01, (n, len(wave))))
    return calib.aggregate(fit, ratios, wave, ["PS1_g", "2MASS_Ks"]), wave


def test_aggregate_bins_and_normalisation():
    corr, wave = _corr()
    t_cool = np.digitize(4600.0, calib.TEFF_EDGES) - 1
    t_hot = np.digitize(5600.0, calib.TEFF_EDGES) - 1
    z = np.digitize(0.0, calib.FEH_EDGES) - 1
    assert corr["n"][t_cool, z] == 30 and corr["n"][t_hot, z] == 30
    R = corr["ratio"][t_cool, z]
    assert abs(np.interp(750, wave, R) - 1.0) < 0.02           # injected ratio is 1 at 700 nm
    assert np.interp(400, wave, R) < 0.95                       # the injected blue deficit
    assert abs(np.interp(400, wave, corr["ratio"][t_hot, z]) - 1.0) < 0.01
    assert corr["band_dm"][t_cool, z, 0] == pytest_approx(0.10, 0.01)
    assert corr["n"][0, z] == 0                                  # empty bin -> no correction


def pytest_approx(v, tol):
    class A:
        def __eq__(self, other):
            return abs(other - v) < tol
    return A()


def test_apply_only_dwarfs_in_populated_bins():
    corr, wave = _corr()
    meta = np.array([[4600.0, 4.5, 0.0], [4600.0, 2.0, 0.0], [3600.0, 4.5, 0.0], [5600.0, 4.5, -1.0]])
    m_xp = np.ones((4, len(wave)))
    fnu0 = np.ones((4, 2))
    mx, fn, used = calib.apply(corr, meta, m_xp, wave, fnu0, ["PS1_g", "2MASS_Ks"])
    assert list(used) == [True, False, False, True]            # giant and empty-bin model untouched
    assert np.interp(400, wave, mx[0]) < 0.95 and mx[1].min() == 1.0 and mx[2].min() == 1.0
    assert fn[0, 0] == pytest_approx(10 ** (-0.4 * 0.10), 0.01) and fn[0, 1] == 1.0
    # metal-poor 5600 K model: its own [Fe/H] bin is empty -> falls back to the T_eff bin
    assert used[3] and abs(mx[3].min() - 1.0) < 0.02


def test_load_disabled(monkeypatch):
    monkeypatch.setenv("DUSTLINE_TEMPLATE_CORR", "none")
    assert calib.load() is None and calib.tag(None) == ""

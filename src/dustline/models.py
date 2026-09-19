"""PHOENIX NewEra model toolkit: spectral cache, XP-resolution resampling,
synthetic photometry on the flexible band registry, and the MIST radius prior
that breaks the T_eff-A_V degeneracy for hot stars.

Flux convention: the NewEra cache carries the surface flux in W m^-2 nm^-1
(x100 = erg s^-1 cm^-2 A^-1; checked against sigma T^4 to 1-4 %).  The
observed flux is (R/D)^2 * F_surface; with R in R_sun and D in kpc the scale
is RD2 = (R_sun/kpc)^2 = 5.086e-22.

The cache and the MIST tables are release assets fetched on first use
(dustline.assets); nothing here reads the network after that.
"""

from __future__ import annotations

import numpy as np

from . import assets, filters

RSUN_CM = 6.957e10
KPC_CM = 3.0857e21
RD2 = (RSUN_CM / KPC_CM) ** 2          # (R/D)^2 for R = 1 R_sun, D = 1 kpc
GMSUN = 1.32712e26                     # cgs G M_sun
C_AA = 2.99792458e18                   # c in Angstrom/s

_cache = {}


def load_cache():
    """(wave_A, flux (Nmod, Nw) float32 W m^-2 nm^-1 surface, meta (Nmod, 3) [Teff, logg, MH])."""
    if "newera" not in _cache:
        d = np.load(assets.fetch("newera_full_cache.npz"))
        _cache["newera"] = (d["wave"].astype(float), d["flux"].astype(np.float32),
                            d["meta"].astype(float))
    return _cache["newera"]


def sub_grid(meta, teff=(3500, 12000), logg=(0.0, 6.0), mh=(0.0,)):
    return np.where((meta[:, 0] >= teff[0]) & (meta[:, 0] <= teff[1]) & (meta[:, 1] >= logg[0])
                    & (meta[:, 1] <= logg[1]) & np.isin(meta[:, 2], mh))[0]


class Photometry:
    """Synthetic photon-counting magnitudes of spectra on a common wavelength grid.

    bands: list of names resolved through the dustline.filters registry.
    """

    def __init__(self, wave_A: np.ndarray, bands: list[str]):
        self.wave = np.asarray(wave_A, float)
        self.bands = list(bands)
        self.T, self.zp, self.norm = {}, {}, {}
        for name in self.bands:
            band = filters.get(name)
            T = np.interp(self.wave, band.wave_A, band.T, left=0.0, right=0.0)
            self.T[name] = T * self.wave               # photon counting: weight lambda T
            self.zp[name] = band.zp_jy
            # f_nu-weighted normalisation: int (c/lam^2) * zp * lam T dlam in f_lam units
            self.norm[name] = np.trapezoid(self.T[name] / self.wave ** 2, self.wave)
        self.dl = np.gradient(self.wave)

    def fnu(self, flux_lam):
        """Band-averaged f_nu (erg/s/cm2/Hz for flux_lam in erg/s/cm2/A), spectra (..., Nw)."""
        out = {}
        for b in self.bands:
            num = (flux_lam * (self.T[b] * self.dl)).sum(axis=-1)
            out[b] = num / (C_AA * self.norm[b])
        return out

    def mags(self, flux_lam):
        """Magnitudes (AB or Vega per band) for observer-frame flux_lam in erg/s/cm2/A."""
        fn = self.fnu(flux_lam)
        return {b: -2.5 * np.log10(fn[b] / (self.zp[b] * 1e-23)) for b in self.bands}


def xp_lsf_matrix(wave_A, xp_wave_nm, fwhm_nm=None):
    """(Nxp, Nw) matrix: Gaussian LSF at each XP sample, integrating the model flux
    density over the cache grid.  Default FWHM: BP/RP resolving power ~30-100
    (Montegriffo+2023): 12 nm at 336, 8 nm at 640, 10 nm at 1020, linear in between."""
    wnm = np.asarray(wave_A, float) / 10.0
    xp_wave_nm = np.asarray(xp_wave_nm, float)
    if fwhm_nm is None:
        fwhm = np.interp(xp_wave_nm, [336.0, 640.0, 1020.0], [12.0, 8.0, 10.0])
    else:
        fwhm = np.full(len(xp_wave_nm), float(fwhm_nm))
    sig = fwhm / 2.3548
    dl = np.gradient(wnm)
    K = np.exp(-0.5 * ((wnm[None, :] - xp_wave_nm[:, None]) / sig[:, None]) ** 2) * dl[None, :]
    K /= K.sum(axis=1, keepdims=True)
    return K.astype(np.float32)


def load_mist():
    """MIST v1.2 basic isochrone table (npz release asset): dict feh -> dict of column arrays."""
    if "mist" not in _cache:
        d = np.load(assets.fetch("mist_v1.2_basic.npz"), allow_pickle=False)
        # layout: for each feh tag (m0.50, p0.00, p0.50): <tag>_cols (names), <tag>_data (N, ncol)
        out = {}
        for tag in ("m0.50", "p0.00", "p0.50"):
            cols = [c for c in d[f"{tag}_cols"]]
            a = d[f"{tag}_data"]
            feh = (-1 if tag[0] == "m" else 1) * float(tag[1:])
            out[feh] = {str(c): a[:, k] for k, c in enumerate(cols)}
        _cache["mist"] = out
    return _cache["mist"]


def radius_prior(meta, feh=0.0, ages=(8.0, 10.1), dlogT=0.02, dlogg=0.2):
    """Expected log10 R/R_sun and its 1-sigma for every model (Teff, logg) from MIST
    isochrones between the given log ages: median and MAD-based sigma (floor 0.08 dex)
    of the isochrone points within dlogT/dlogg; fallback: R from log g with M = 1 M_sun,
    sigma 0.15."""
    tables = load_mist()
    iso = tables[min(tables, key=lambda f: abs(f - feh))]
    sel = ((iso["log10_isochrone_age_yr"] >= ages[0]) & (iso["log10_isochrone_age_yr"] <= ages[1])
           & (iso["phase"] <= 4))
    lT, lg, lR = iso["log_Teff"][sel], iso["log_g"][sel], iso["log_R"][sel]
    out = np.zeros((len(meta), 2))
    for k, (teff, logg, _) in enumerate(meta):
        m = (np.abs(lT - np.log10(teff)) < dlogT) & (np.abs(lg - logg) < dlogg)
        if m.sum() >= 3:
            med = np.median(lR[m])
            sig = max(1.4826 * np.median(np.abs(lR[m] - med)), 0.08)
        else:
            med = 0.5 * (np.log10(GMSUN) - logg) - np.log10(RSUN_CM)
            sig = 0.15
        out[k] = med, sig
    return out

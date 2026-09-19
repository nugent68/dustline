"""Filter/band registry.

A *band* is a name bound to a transmission curve file and a magnitude system.
Built-in curves ship as package data (``dustline/data/filters/*.dat``, SVO
photon-counting transmissions with a ``ZeroPoint`` header in Jy for Vega
systems). Users may register their own curves at runtime.

    from dustline import filters
    filters.get("DECam_i")          # -> Band(wave_A, T, system, zp_jy)
    filters.register("myband", "/path/to/curve.dat", system="AB")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import assets


@dataclass(frozen=True)
class Band:
    name: str
    wave_A: np.ndarray = field(repr=False)
    T: np.ndarray = field(repr=False)
    system: str            # "AB" | "Vega"
    zp_jy: float           # 3631 for AB, the file header value for Vega

    @property
    def lam_eff_A(self) -> float:
        """Photon-weighted effective wavelength (A) of the curve itself."""
        w = self.T * self.wave_A
        return float((self.wave_A * w).sum() / w.sum())


# name -> (curve file stem, system); the packaged curves
BUILTIN: dict[str, tuple[str, str]] = {
    # DECam (AB)
    "DECam_g": ("DECam_g", "AB"), "DECam_r": ("DECam_r", "AB"), "DECam_i": ("DECam_i", "AB"),
    "DECam_z": ("DECam_z", "AB"), "DECam_Y": ("DECam_Y", "AB"),
    # Pan-STARRS1 (AB)
    "PS1_g": ("PS1_g", "AB"), "PS1_r": ("PS1_r", "AB"), "PS1_i": ("PS1_i", "AB"),
    "PS1_z": ("PS1_z", "AB"), "PS1_y": ("PS1_y", "AB"),
    # 2MASS / VISTA (Vega)
    "2MASS_J": ("2MASS_J", "Vega"), "2MASS_H": ("2MASS_H", "Vega"), "2MASS_Ks": ("2MASS_Ks", "Vega"),
    "VISTA_J": ("VISTA_J", "Vega"), "VISTA_H": ("VISTA_H", "Vega"), "VISTA_Ks": ("VISTA_Ks", "Vega"),
    # Gaia (Vega)
    "Gaia_G": ("Gaia_G", "Vega"), "Gaia_BP": ("Gaia_BP", "Vega"), "Gaia_RP": ("Gaia_RP", "Vega"),
    # Johnson-Cousins (Vega) -- the default output filter I
    "I": ("Cousins_I", "Vega"), "Cousins_I": ("Cousins_I", "Vega"),
    "V": ("Johnson_V", "Vega"), "Johnson_V": ("Johnson_V", "Vega"),
    # Keck
    "Keck_Kp": ("Keck_Kp", "Vega"),
    # WISE
    "WISE_W1": ("WISE_W1", "Vega"), "WISE_W2": ("WISE_W2", "Vega"),
}

_registry: dict[str, Band] = {}


def _read_curve(path: Path) -> tuple[np.ndarray, np.ndarray, float]:
    head = open(path).readline()
    zp = float(head.split("ZeroPoint")[1].split()[0]) if "ZeroPoint" in head else np.nan
    w, t = np.loadtxt(path).T
    return w, t, zp


def register(name: str, path: str | Path, system: str = "AB", zp_jy: float | None = None) -> Band:
    """Register a user filter curve (two-column wave_A, T; SVO .dat works as is)."""
    w, t, zp_file = _read_curve(Path(path))
    zp = 3631.0 if system == "AB" else (zp_jy if zp_jy is not None else zp_file)
    if not np.isfinite(zp):
        raise ValueError(f"band {name!r}: Vega system needs zp_jy (no ZeroPoint header found)")
    band = Band(name, w, t, system, zp)
    _registry[name] = band
    return band


def get(name: str) -> Band:
    """Look up a band by name, loading a built-in curve on first use."""
    if name in _registry:
        return _registry[name]
    if name not in BUILTIN:
        raise KeyError(f"unknown band {name!r}; built-ins: {sorted(BUILTIN)}; use filters.register()")
    stem, system = BUILTIN[name]
    path = assets.filter_dir() / f"{stem}.dat"
    w, t, zp_file = _read_curve(path)
    band = Band(name, w, t, system, 3631.0 if system == "AB" else zp_file)
    _registry[name] = band
    return band


def known() -> list[str]:
    return sorted(set(BUILTIN) | set(_registry))

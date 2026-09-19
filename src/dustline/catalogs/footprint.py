"""Survey footprint logic: which photometric catalogs cover a sightline.

The defaults are the simplest all-sky(ish) pair: PS1 DR2 grizy (Dec > -30) and
2MASS JHKs (all sky, joined in the Gaia query). Deeper/redder surveys are used
where their footprints allow:

- DECaPS DR2 (DECam grizY): the southern Galactic plane, |b| < ~10, 6 > l or
  l > 239 (Saydjari+2023 coverage: -124 < l < 6 approximately, |b| <~ 10).
- VVV/VVVX (VISTA JHKs): the bulge -10 < l < 10, -10 < b < 5 and the southern
  disk 295 < l < 350, |b| < 2 (approximate DR footprints).
- DECaLS (DECam grz): Dec < 32, |b| > ~15 (extragalactic program) - a fallback
  for southern high-latitude fields where PS1 is absent (Dec < -30).
"""

from __future__ import annotations

from dataclasses import dataclass

from astropy.coordinates import SkyCoord
import astropy.units as u


@dataclass
class SurveyPlan:
    """Which catalogs to query for a sightline, in slot order."""

    optical: str            # "ps1" | "decaps" | "decals" | "none"
    nir: str                # "2mass" | "vvv"
    in_bulge_window: bool   # red-clump anchor possible (|l| < 10, |b| < 10)
    notes: list[str]


def galactic(ra: float, dec: float) -> tuple[float, float]:
    c = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs").galactic
    return float(c.l.deg), float(c.b.deg)


def plan(ra: float, dec: float, prefer_deep: bool = True) -> SurveyPlan:
    """Choose the survey set for a position."""
    lon, lat = galactic(ra, dec)
    lw = lon if lon <= 180 else lon - 360.0    # wrap to (-180, 180]
    notes: list[str] = []

    # --- optical ---
    optical = "none"
    in_decaps = (abs(lat) < 10.0) and (lw < 6.0) and (lw > -124.0)
    if in_decaps and prefer_deep:
        optical = "decaps"
        notes.append("DECaPS DR2 grizY (southern Galactic plane)")
    elif dec > -30.0:
        optical = "ps1"
        notes.append("PS1 DR2 grizy")
    elif dec <= -30.0 and abs(lat) > 15.0:
        optical = "decals"
        notes.append("DECaLS grz (Dec < -30, high latitude)")
    else:
        notes.append("no optical survey coverage: supply user photometry "
                     "(Dec < -30, low latitude, outside DECaPS)")

    # --- near-infrared ---
    in_vvv = ((-10.0 < lw < 10.5) and (-10.3 < lat < 5.1)) or \
             ((294.7 - 360.0 < lw < 350.0 - 360.0) and (abs(lat) < 2.25))
    if in_vvv and prefer_deep:
        nir = "vvv"
        notes.append("VVV JHKs (2MASS fallback for bright/saturated stars)")
    else:
        nir = "2mass"
        notes.append("2MASS JHKs")

    in_bulge = (abs(lw) < 10.0) and (abs(lat) < 10.0)
    return SurveyPlan(optical=optical, nir=nir, in_bulge_window=in_bulge, notes=notes)

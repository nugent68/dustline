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

Beyond photometry the plan also records
- ``spectro``: whether DESI DR1 Milky Way Survey stellar parameters (T_eff,
  log g, [Fe/H]) are available as per-star template priors (the DR1 footprint:
  high latitude, Dec > -25 approximately);
- ``mode``: "law" (the bulge/plane product: measured R_V + A_V(D) to kpc) or
  "column" (high latitude, |b| > 30: too little dust to measure R_V; the product
  is the total foreground column from the stars behind the dust).
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
    spectro: str = "none"   # "desi" | "none": per-star spectroscopic template priors
    mode: str = "law"       # "law" | "column"
    apogee: bool = True     # APOGEE DR17 ASPCAP parameters (all sky; external T_eff check)
    uv: str = "none"        # "galex" | "none": GALEX GUVcat FUV/NUV (needs the UV-IR model cache)
    mir: str = "wise"       # "wise" | "none": AllWISE W1/W2 (needs the UV-IR model cache)


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

    # red-clump anchor window: the bulge + long bar reach |l| ~ 20 in the plane;
    # clump.find_clump applies its own credibility gates (window population,
    # E(J-Ks) > 0.15, D_RC within 5-12 kpc), so a generous window is safe.
    in_bulge = (abs(lw) < 20.0) and (abs(lat) < 10.0)

    # --- DESI DR1 MWS stellar parameters: the high-latitude DESI footprint ---
    # (Dec > -25, |b| > 15 is a generous cut; an empty query is harmless)
    spectro = "desi" if (dec > -25.0 and abs(lat) > 15.0) else "none"
    if spectro == "desi":
        notes.append("DESI DR1 MWS stellar parameters as template priors")

    # --- product mode: no measurable law at high latitude ---
    mode = "column" if abs(lat) > 30.0 else "law"
    if mode == "column":
        notes.append("high latitude: foreground-column mode (R_V assumed 3.1)")

    # --- UV (GALEX AIS avoids the plane) and mid-IR (AllWISE, confused in the bulge) ---
    uv = "galex" if abs(lat) > 20.0 else "none"
    mir = "wise" if abs(lat) > 10.0 else "none"
    if uv == "galex":
        notes.append("GALEX GUVcat_AIS FUV/NUV (hot stars; UV-IR model cache)")
    if mir == "wise":
        notes.append("AllWISE W1/W2 (UV-IR model cache)")
    return SurveyPlan(optical=optical, nir=nir, in_bulge_window=in_bulge, notes=notes,
                      spectro=spectro, mode=mode, uv=uv, mir=mir)

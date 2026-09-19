"""User-supplied calibrated photometry.

The CSV needs ``ra, dec`` columns (deg, ICRS, Gaia-compatible frame) and pairs
``mag_<key>, magerr_<key>``. ``bands`` maps each key to a filter registry name
(built-in like "DECam_i", or one you registered with dustline.filters.register).

    phot = UserPhotometry("my_phot.csv", bands={"g": "DECam_g", "Ks": "VISTA_Ks"})
    res = dustline.Sightline(ra, dec, photometry=phot).run()
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import filters
from .catalogs.xmatch import match_to_gaia

MATCH_ARCSEC = 0.7


class UserPhotometry:
    def __init__(self, path: str | Path, bands: dict[str, str],
                 match_arcsec: float = MATCH_ARCSEC):
        self.path = Path(path)
        self.bands = dict(bands)
        self.match_arcsec = float(match_arcsec)
        self.table = pd.read_csv(self.path)
        missing = [c for c in ("ra", "dec") if c not in self.table.columns]
        if missing:
            raise ValueError(f"user photometry {self.path}: missing columns {missing}")
        for key, name in self.bands.items():
            filters.get(name)   # raises early on unknown band names
            for c in (f"mag_{key}", f"magerr_{key}"):
                if c not in self.table.columns:
                    raise ValueError(f"user photometry {self.path}: missing column {c}")

    def matched_to_gaia(self, gaia: pd.DataFrame, ra0: float, dec0: float) -> pd.DataFrame:
        """Per-Gaia-source table with columns renamed to the registry band names."""
        m = match_to_gaia(gaia, self.table.ra.values, self.table.dec.values,
                          ra0, dec0, self.match_arcsec)
        p = self.table.iloc[m.icat.values]
        out = pd.DataFrame(dict(source_id=m.source_id.values))
        for key, name in self.bands.items():
            out[f"mag_{name}"] = p[f"mag_{key}"].values
            out[f"magerr_{name}"] = p[f"magerr_{key}"].values
        print(f"user photometry: {len(out)} Gaia matches, bands {sorted(self.bands.values())}")
        return out

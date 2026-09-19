"""Position cross-match helper: nearest survey object per Gaia source within a
tolerance, in the tangent plane around the field centre."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def match_to_gaia(gaia: pd.DataFrame, cat_ra: np.ndarray, cat_dec: np.ndarray,
                  ra0: float, dec0: float, tol_arcsec: float) -> pd.DataFrame:
    """DataFrame(source_id, icat, sep) of unique nearest matches within tol."""
    cosd = np.cos(np.radians(dec0))
    gx = np.column_stack([(gaia.ra.values - ra0) * cosd, gaia.dec.values - dec0]) * 3600
    px = np.column_stack([(np.asarray(cat_ra) - ra0) * cosd, np.asarray(cat_dec) - dec0]) * 3600
    d, j = cKDTree(px).query(gx, distance_upper_bound=tol_arcsec)
    ok = np.isfinite(d)
    m = pd.DataFrame(dict(source_id=gaia.source_id.values[ok], icat=j[ok], sep=d[ok]))
    return m.sort_values("sep").drop_duplicates("icat")

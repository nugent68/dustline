"""Pan-STARRS1 DR2 mean PSF photometry (MAST catalogs API), matched to the Gaia
table at 0.7 arcsec.  Returns per-Gaia-source magnitudes in the PS1 grizy bands
(columns mag_PS1_g .. magerr_PS1_y).

The API's paging is not stable (rows are not ordered between pages: the same
query returned page overlaps of 0, 180 and 6,936 rows in three trials, i.e.
objects silently dropped or duplicated at every page boundary).  So no paging:
a cone that fills a page is split into seven sub-cones, recursively, and the
tiles are deduplicated on objID.
"""

from __future__ import annotations

import io
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

from ..cache import Workspace
from .xmatch import match_to_gaia

API = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv"
COLS = ["objID", "raMean", "decMean", "nDetections", "qualityFlag"] + \
    [f"{b}{s}" for b in "grizy" for s in ("MeanPSFMag", "MeanPSFMagErr", "MeanPSFMagNpt", "QfPerfect")]
MATCH_ARCSEC = 0.7
SYS_FLOOR = 0.02


PAGESIZE = 50000


def _fetch_cone(ra: float, dec: float, radius_deg: float) -> pd.DataFrame:
    q = urllib.parse.urlencode({"ra": ra, "dec": dec, "radius": radius_deg,
                                "nDetections.gte": 3, "pagesize": PAGESIZE, "page": 1,
                                "columns": ",".join(COLS)})
    txt = urllib.request.urlopen(API + "?" + q, timeout=600).read().decode()
    return pd.read_csv(io.StringIO(txt))


def _fetch_all(ra: float, dec: float, radius_deg: float, depth: int = 0) -> pd.DataFrame:
    """All PS1 objects in the cone without paging: a full page means the cone is
    split into a centre tile and six around it (radius 0.65 r at distance 0.5 r
    covers the disk), fetched recursively and deduplicated on objID."""
    df = _fetch_cone(ra, dec, radius_deg)
    if len(df) < PAGESIZE or depth > 6:
        return df.replace(-999.0, np.nan) if len(df) else pd.DataFrame(columns=COLS)
    cosd = max(np.cos(np.radians(dec)), 1e-3)
    tiles = [(ra, dec)] + [(ra + 0.5 * radius_deg * np.cos(t) / cosd,
                            dec + 0.5 * radius_deg * np.sin(t))
                           for t in np.radians(np.arange(0, 360, 60))]
    parts = [_fetch_all(r, d, 0.65 * radius_deg, depth + 1) for r, d in tiles]
    out = pd.concat([p for p in parts if len(p)], ignore_index=True)
    out = out.drop_duplicates("objID").reset_index(drop=True)
    # cut the union of tiles back to the requested cone
    sep = np.hypot((out.raMean - ra) * cosd, out.decMean - dec)
    return out[sep <= radius_deg].reset_index(drop=True)


def fetch(ws: Workspace, gaia: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    """PS1 grizy per Gaia source -> ps1_gaia.csv (cached)."""
    out = ws.path("ps1_gaia.csv")
    if out.exists() and not force:
        return pd.read_csv(out)
    ps = _fetch_all(ws.ra, ws.dec, ws.radius_arcmin / 60.0)
    if not len(ps):
        print("PS1: no coverage / no objects")
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    m = match_to_gaia(gaia, ps.raMean.values, ps.decMean.values,
                      ws.ra, ws.dec, MATCH_ARCSEC)
    p = ps.iloc[m.icat.values]
    res = pd.DataFrame(dict(source_id=m.source_id.values, ps1_sep=m.sep.values,
                            ps1_nDet=p.nDetections.values))
    for b in "grizy":
        good = (p[f"{b}QfPerfect"].values > 0.85) & (p[f"{b}MeanPSFMagNpt"].values >= 2)
        res[f"mag_PS1_{b}"] = np.where(good, p[f"{b}MeanPSFMag"].values, np.nan)
        res[f"magerr_PS1_{b}"] = np.where(
            good, np.sqrt(p[f"{b}MeanPSFMagErr"].values ** 2 + SYS_FLOOR ** 2), np.nan)
    res.round(4).to_csv(out, index=False)
    print(f"PS1: {len(res)} Gaia matches")
    return res

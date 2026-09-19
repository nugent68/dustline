"""Pan-STARRS1 DR2 mean PSF photometry (MAST catalogs API, paged), matched to
the Gaia table at 0.7 arcsec.  Returns per-Gaia-source magnitudes in the PS1
grizy bands (columns mag_PS1_g .. magerr_PS1_y).
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


def _fetch_all(ra: float, dec: float, radius_deg: float) -> pd.DataFrame:
    frames, page = [], 1
    while True:
        q = urllib.parse.urlencode({"ra": ra, "dec": dec, "radius": radius_deg,
                                    "nDetections.gte": 3, "pagesize": 50000, "page": page,
                                    "columns": ",".join(COLS)})
        txt = urllib.request.urlopen(API + "?" + q, timeout=600).read().decode()
        df = pd.read_csv(io.StringIO(txt))
        if not len(df):
            break
        frames.append(df)
        if len(df) < 50000:
            break
        page += 1
    if not frames:
        return pd.DataFrame(columns=COLS)
    return pd.concat(frames, ignore_index=True).replace(-999.0, np.nan)


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

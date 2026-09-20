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
# mean PSF magnitudes brighter than this are saturated (biased faint by 0.1-0.2 mag)
SATURATION = {"g": 14.5, "r": 15.0, "i": 15.0, "z": 14.0, "y": 13.0}


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


def _unsaturate(res: pd.DataFrame) -> pd.DataFrame:
    """NaN the saturated mean-PSF magnitudes (also for tables cached before the cut)."""
    for b, lim in SATURATION.items():
        if f"mag_PS1_{b}" in res:
            sat = res[f"mag_PS1_{b}"] <= lim
            res.loc[sat, [f"mag_PS1_{b}", f"magerr_PS1_{b}"]] = np.nan
    return res


def fetch(ws: Workspace, gaia: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    """PS1 grizy per Gaia source -> ps1_gaia.csv (cached)."""
    out = ws.path("ps1_gaia.csv")
    if out.exists() and not force:
        return _unsaturate(pd.read_csv(out))
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
        good = ((p[f"{b}QfPerfect"].values > 0.85) & (p[f"{b}MeanPSFMagNpt"].values >= 2)
                & (p[f"{b}MeanPSFMag"].values > SATURATION[b]))
        res[f"mag_PS1_{b}"] = np.where(good, p[f"{b}MeanPSFMag"].values, np.nan)
        res[f"magerr_PS1_{b}"] = np.where(
            good, np.sqrt(p[f"{b}MeanPSFMagErr"].values ** 2 + SYS_FLOOR ** 2), np.nan)
    res.round(4).to_csv(out, index=False)
    print(f"PS1: {len(res)} Gaia matches")
    return res


def by_positions(ra: np.ndarray, dec: np.ndarray, radius_arcsec: float = 2.0,
                 threads: int = 4) -> pd.DataFrame:
    """Nearest PS1 DR2 object (nDetections >= 3) for each position: one small cone
    query per star (the crossmatch-upload endpoint is broken), threaded.  Returns
    one row per input index (NaNs when nothing is found) with the mean-PSF columns."""
    from concurrent.futures import ThreadPoolExecutor

    def one(i):
        try:
            q = urllib.parse.urlencode({"ra": ra[i], "dec": dec[i], "radius": radius_arcsec / 3600.0,
                                        "nDetections.gte": 3, "pagesize": 10,
                                        "columns": ",".join(COLS)})
            txt = urllib.request.urlopen(API + "?" + q, timeout=120).read().decode()
            df = pd.read_csv(io.StringIO(txt))
            if not len(df):
                return i, None
            df = df.replace(-999.0, np.nan)
            cosd = np.cos(np.radians(dec[i]))
            sep = np.hypot((df.raMean - ra[i]) * cosd, df.decMean - dec[i]) * 3600
            return i, df.iloc[int(np.argmin(sep.values))][COLS]
        except Exception:  # noqa: BLE001
            return i, None

    rows = [None] * len(ra)
    with ThreadPoolExecutor(threads) as ex:
        for k, (i, r) in enumerate(ex.map(one, range(len(ra))), 1):
            rows[i] = r
            if k % 200 == 0:
                print(f"  ps1 by_positions: {k}/{len(ra)}", flush=True)
    out = pd.DataFrame([r if r is not None else pd.Series(np.nan, index=COLS) for r in rows])
    return out.reset_index(drop=True)

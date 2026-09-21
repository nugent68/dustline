"""VVV JHKs photometry matched to Gaia.  Delivers mag_VISTA_J/H/Ks.

Data Lab has no VVV source table (checked 2026-09-21: only VHS, which excludes
the VVV footprint), so the VVV DR2 JHKs are taken from the DECaPS DR2 brutus
catalogue ``decaps_dr2.stellar_inference`` (Zucker+2025), whose columns mag_6/7/8
(magerr_6/7/8) are the VVV J/H/Ks of each DECaPS object, with its Gaia DR3 id
(``gaia_id``; position match at 0.7" where it is 0).  The original
``vvv_dr4.vvvsource`` query is kept as a fallback should the table appear.
Deeper than 2MASS in the bulge/southern-disk window; 2MASS remains the fallback
for bright stars saturated in VVV (handled at band assembly).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..cache import Workspace
from .datalab import box_query
from .xmatch import match_to_gaia

MATCH_ARCSEC = 0.7
SYS_FLOOR = 0.03
_BANDS = {"j": "VISTA_J", "h": "VISTA_H", "ks": "VISTA_Ks"}


_BRUTUS = {"j": 6, "h": 7, "ks": 8}


def fetch_brutus(ws: Workspace, gaia: pd.DataFrame) -> pd.DataFrame | None:
    """VVV JHKs via decaps_dr2.stellar_inference (mag_6..8); None if the query fails."""
    cols = ["ra", "dec", "gaia_id"] + [f"mag_{k}" for k in _BRUTUS.values()] + \
        [f"magerr_{k}" for k in _BRUTUS.values()] + ["decaps_fracflux_3"]
    try:
        d = box_query("decaps_dr2.stellar_inference", cols, ws.ra, ws.dec, ws.radius_arcmin / 60.0)
    except Exception as e:  # noqa: BLE001
        print(f"VVV (brutus) query failed ({e})")
        return None
    if not len(d):
        return None
    d = d.replace([np.inf, -np.inf], np.nan)
    # match: the catalogue's own Gaia id where present (best fracflux per id), else position
    d["ok_id"] = (d.gaia_id > 0) & d.gaia_id.isin(gaia.source_id.values)
    by_id = d[d.ok_id].sort_values("decaps_fracflux_3", ascending=False).drop_duplicates("gaia_id")
    rest = d[~d.ok_id]
    parts = [pd.DataFrame(dict(source_id=by_id.gaia_id.astype(np.int64).values, vvv_sep=0.0, _row=by_id.index.values))]
    if len(rest):
        g2 = gaia[~gaia.source_id.isin(by_id.gaia_id.values)]
        m = match_to_gaia(g2, rest.ra.values, rest.dec.values, ws.ra, ws.dec, MATCH_ARCSEC)
        parts.append(pd.DataFrame(dict(source_id=m.source_id.values, vvv_sep=m.sep.values, _row=rest.index.values[m.icat.values])))
    mm = pd.concat(parts, ignore_index=True).drop_duplicates("source_id")
    p = d.loc[mm._row.values]
    res = pd.DataFrame(dict(source_id=mm.source_id.values, vvv_sep=mm.vvv_sep.values))
    for b, name in _BANDS.items():
        k = _BRUTUS[b]
        mag = p[f"mag_{k}"].values.astype(float); err = p[f"magerr_{k}"].values.astype(float)
        good = np.isfinite(mag) & np.isfinite(err) & (err < 0.3)
        res[f"mag_{name}"] = np.where(good, mag, np.nan)
        res[f"magerr_{name}"] = np.where(good, np.sqrt(err ** 2 + SYS_FLOOR ** 2), np.nan)
    return res


def fetch(ws: Workspace, gaia: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    out = ws.path("vvv_gaia.csv")
    if out.exists() and not force:
        return pd.read_csv(out)
    res = fetch_brutus(ws, gaia)
    if res is not None and len(res):
        res.round(4).to_csv(out, index=False)
        print(f"VVV (via decaps_dr2.stellar_inference): {len(res)} Gaia matches, "
              f"{int(np.isfinite(res.mag_VISTA_Ks).sum())} with Ks")
        return res
    cols = ["ra2000 AS ra", "dec2000 AS dec"] + \
        [f"{b}apermag3 AS mag_{b}" for b in _BANDS] + \
        [f"{b}apermag3err AS magerr_{b}" for b in _BANDS] + \
        [f"{b}ppErrBits AS bits_{b}" for b in _BANDS]
    try:
        d = box_query("vvv_dr4.vvvsource", cols, ws.ra, ws.dec, ws.radius_arcmin / 60.0,
                      ra_col="ra2000", dec_col="dec2000")
        d.columns = [c.lower() for c in d.columns]
    except Exception as e:  # noqa: BLE001
        print(f"VVV query failed ({e}); falling back to 2MASS only")
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    if not len(d):
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    m = match_to_gaia(gaia, d.ra.values, d.dec.values, ws.ra, ws.dec, MATCH_ARCSEC)
    p = d.iloc[m.icat.values]
    res = pd.DataFrame(dict(source_id=m.source_id.values, vvv_sep=m.sep.values))
    for b, name in _BANDS.items():
        mag = p[f"mag_{b}"].values
        err = p[f"magerr_{b}"].values
        bits = p[f"bits_{b}"].values if f"bits_{b}" in p else np.zeros(len(p))
        good = np.isfinite(mag) & np.isfinite(err) & (err < 0.3) & (bits < 256)
        res[f"mag_{name}"] = np.where(good, mag, np.nan)
        res[f"magerr_{name}"] = np.where(good, np.sqrt(err ** 2 + SYS_FLOOR ** 2), np.nan)
    res.round(4).to_csv(out, index=False)
    print(f"VVV: {len(res)} Gaia matches")
    return res

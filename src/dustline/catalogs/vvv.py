"""VVV DR4 JHKs photometry (NOIRLab Data Lab vhs-like schema), matched to Gaia.

Deeper than 2MASS in the bulge/southern-disk window; 2MASS remains the
fallback for bright stars saturated in VVV (handled at band assembly).
Delivers mag_VISTA_J/H/Ks.
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


def fetch(ws: Workspace, gaia: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    out = ws.path("vvv_gaia.csv")
    if out.exists() and not force:
        return pd.read_csv(out)
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

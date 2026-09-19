"""AllWISE W1/W2 photometry (NOIRLab Data Lab ``allwise.source``) matched to
Gaia: mag_WISE_W1 / mag_WISE_W2 (Vega).

Not an extinction lever (A_W1/A_V ~ 0.06) but an almost extinction-free anchor
on the Rayleigh-Jeans tail, i.e. on the flux scale (R/D)^2 and, through the
radius prior, on T_eff - which helps the stars without a spectroscopic prior
and checks the frozen 2MASS zero points.  W1/W2 zero points are frozen after
the first pass like the other NIR bands (fit.NIR_PREFIXES).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..cache import Workspace
from .datalab import box_query
from .xmatch import match_to_gaia

MATCH_ARCSEC = 1.5
SYS_FLOOR = 0.03
MAX_ERR = 0.3
SAT_MAG = 8.0             # brighter than this W1/W2 are saturated
_BANDS = {"w1": "WISE_W1", "w2": "WISE_W2"}
_COLS = ["ra", "dec", "w1mpro", "w1sigmpro", "w2mpro", "w2sigmpro", "cc_flags", "ext_flg"]


def shape(d: pd.DataFrame, gaia: pd.DataFrame, ra0: float, dec0: float) -> pd.DataFrame:
    m = match_to_gaia(gaia, d.ra.values, d.dec.values, ra0, dec0, MATCH_ARCSEC)
    p = d.iloc[m.icat.values]
    res = pd.DataFrame(dict(source_id=m.source_id.values, wise_sep=m.sep.values))
    cc = p.cc_flags.astype(str).str.ljust(4, "0").values
    ext = p.ext_flg.fillna(0).astype(int).values
    for k, (b, name) in enumerate(_BANDS.items()):
        mag = p[f"{b}mpro"].values.astype(float)
        err = p[f"{b}sigmpro"].values.astype(float)
        good = (np.isfinite(mag) & np.isfinite(err) & (err < MAX_ERR) & (mag > SAT_MAG)
                & (np.array([c[k] for c in cc]) == "0") & (ext == 0))
        res[f"mag_{name}"] = np.where(good, mag, np.nan)
        res[f"magerr_{name}"] = np.where(good, np.sqrt(err ** 2 + SYS_FLOOR ** 2), np.nan)
    return res


def fetch(ws: Workspace, gaia: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    out = ws.path("wise_gaia.csv")
    if out.exists() and not force:
        return pd.read_csv(out)
    try:
        d = box_query("allwise.source", _COLS, ws.ra, ws.dec, ws.radius_arcmin / 60.0)
    except Exception as e:  # noqa: BLE001
        print(f"AllWISE query failed ({e}); continuing without")
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    if not len(d):
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    res = shape(d, gaia, ws.ra, ws.dec)
    res.round(4).to_csv(out, index=False)
    n = {b: int(np.isfinite(res[f"mag_{b}"]).sum()) for b in _BANDS.values()}
    print(f"AllWISE: {len(res)} Gaia matches, usable {n}")
    return res

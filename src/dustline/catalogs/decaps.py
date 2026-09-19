"""DECaPS DR2 mean photometry (NOIRLab Data Lab), matched to Gaia at 0.7".

Columns delivered: mag_DECam_g .. magerr_DECam_Y (grizY; the DECaPS 'y' filter
column is the DECam Y band).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..cache import Workspace
from .datalab import box_query
from .xmatch import match_to_gaia

MATCH_ARCSEC = 0.7
SYS_FLOOR = 0.02
_BANDS = {"g": "DECam_g", "r": "DECam_r", "i": "DECam_i", "z": "DECam_z", "y": "DECam_Y"}


def fetch(ws: Workspace, gaia: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    out = ws.path("decaps_gaia.csv")
    if out.exists() and not force:
        return pd.read_csv(out)
    cols = ["ra", "dec", "posstdev"] + \
        [f"mean_mag_{b}" for b in _BANDS] + [f"err_mag_{b}" for b in _BANDS] + \
        [f"nmag_{b}" for b in _BANDS] + [f"fracflux_avg_{b}" for b in _BANDS]
    try:
        d = box_query("decaps_dr2.object", cols, ws.ra, ws.dec, ws.radius_arcmin / 60.0)
    except Exception as e:  # noqa: BLE001
        print(f"DECaPS query failed ({e}); continuing without")
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    if not len(d):
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    m = match_to_gaia(gaia, d.ra.values, d.dec.values, ws.ra, ws.dec, MATCH_ARCSEC)
    p = d.iloc[m.icat.values]
    res = pd.DataFrame(dict(source_id=m.source_id.values, decaps_sep=m.sep.values))
    for b, name in _BANDS.items():
        good = (p[f"nmag_{b}"].values >= 2) & (p[f"fracflux_avg_{b}"].values > 0.75)
        res[f"mag_{name}"] = np.where(good, p[f"mean_mag_{b}"].values, np.nan)
        res[f"magerr_{name}"] = np.where(
            good, np.sqrt(p[f"err_mag_{b}"].values ** 2 + SYS_FLOOR ** 2), np.nan)
    res.round(4).to_csv(out, index=False)
    print(f"DECaPS: {len(res)} Gaia matches")
    return res

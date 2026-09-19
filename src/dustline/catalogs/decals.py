"""DECaLS / Legacy Surveys DR10 photometry (NOIRLab Data Lab), matched to Gaia.

Used as the optical fallback for Dec < -30 fields outside DECaPS. Delivers
mag_DECam_g/r/z (i where present in DR10 south).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..cache import Workspace
from .datalab import box_query
from .xmatch import match_to_gaia

MATCH_ARCSEC = 0.7
SYS_FLOOR = 0.02
_BANDS = {"g": "DECam_g", "r": "DECam_r", "i": "DECam_i", "z": "DECam_z"}


def fetch(ws: Workspace, gaia: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    out = ws.path("decals_gaia.csv")
    if out.exists() and not force:
        return pd.read_csv(out)
    cols = ["ra", "dec", "type"] + [f"mag_{b}" for b in _BANDS] + \
        [f"snr_{b}" for b in _BANDS] + [f"fracflux_{b}" for b in _BANDS]
    try:
        d = box_query("ls_dr10.tractor", cols, ws.ra, ws.dec, ws.radius_arcmin / 60.0,
                      extra="AND type = 'PSF'")
    except Exception as e:  # noqa: BLE001
        print(f"DECaLS query failed ({e}); continuing without")
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    if not len(d):
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    m = match_to_gaia(gaia, d.ra.values, d.dec.values, ws.ra, ws.dec, MATCH_ARCSEC)
    p = d.iloc[m.icat.values]
    res = pd.DataFrame(dict(source_id=m.source_id.values, decals_sep=m.sep.values))
    for b, name in _BANDS.items():
        snr = p[f"snr_{b}"].values
        good = np.isfinite(p[f"mag_{b}"].values) & (snr > 5) & (p[f"fracflux_{b}"].values < 0.25)
        err = np.where(snr > 0, 1.0857 / np.maximum(snr, 1e-3), np.nan)
        res[f"mag_{name}"] = np.where(good, p[f"mag_{b}"].values, np.nan)
        res[f"magerr_{name}"] = np.where(good, np.sqrt(err ** 2 + SYS_FLOOR ** 2), np.nan)
    res.round(4).to_csv(out, index=False)
    print(f"DECaLS: {len(res)} Gaia matches")
    return res

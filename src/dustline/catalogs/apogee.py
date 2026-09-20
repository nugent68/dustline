"""APOGEE DR17 ASPCAP parameters (NOIRLab Data Lab ``sdss_dr17.apogee2_allstar``),
joined to Gaia by the catalogue's Gaia EDR3 source_id: ``teff_ap``, ``logg_ap``,
``feh_ap`` (+ errors, S/N).

APOGEE's T_eff is calibrated to the IRFM scale (Gonzalez Hernandez & Bonifacio
2009), unlike DESI RVSpecFit's, so it is the external check of the fitted
T_eff in reddened fields; where a star has no DESI parameters its log g and
[Fe/H] serve as the priors.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..cache import Workspace
from .datalab import box_query

_COLS = ["ra", "dec", "gaiaedr3_source_id", "teff", "logg", "fe_h", "fe_h_err", "m_h",
         "alpha_m", "snr", "starflag", "telescope"]
_EXTRA = "AND teff > 0 AND snr > 30"


def shape(d: pd.DataFrame) -> pd.DataFrame:
    d = d[np.isfinite(d.gaiaedr3_source_id) & (d.gaiaedr3_source_id > 0)].copy()
    d["source_id"] = d.gaiaedr3_source_id.astype(np.int64)
    d = d.sort_values("snr", ascending=False).drop_duplicates("source_id")
    return pd.DataFrame(dict(source_id=d.source_id.values, teff_ap=d.teff.values,
                             logg_ap=d.logg.values, feh_ap=d.fe_h.values,
                             feh_ap_err=d.fe_h_err.values, alphafe_ap=d.alpha_m.values,
                             ap_snr=d.snr.values))


def fetch(ws: Workspace, gaia: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    out = ws.path("apogee_gaia.csv")
    if out.exists() and not force:
        return pd.read_csv(out)
    try:
        d = box_query("sdss_dr17.apogee2_allstar", _COLS, ws.ra, ws.dec, ws.radius_arcmin / 60.0,
                      extra=_EXTRA)
    except Exception as e:  # noqa: BLE001
        print(f"APOGEE query failed ({e}); continuing without")
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    res = shape(d) if len(d) else pd.DataFrame(dict(source_id=[]))
    if len(res):
        res = res[res.source_id.isin(gaia.source_id.values)].reset_index(drop=True)
    res.round(4).to_csv(out, index=False)
    print(f"APOGEE DR17: {len(res)} Gaia stars with ASPCAP parameters")
    return res

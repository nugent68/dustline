"""GALEX GUVcat_AIS (Bianchi, Shiao & Thilker 2017; VizieR II/335) FUV/NUV
photometry matched to Gaia: mag_GALEX_FUV / mag_GALEX_NUV (AB).

UV leverage on the extinction: under G23 at R_V 3.1, A_NUV/A_V = 2.9 and
A_FUV/A_V = 2.8 (against 1.2 in g), so A_V = 0.05 is 0.15 mag in the UV.  The
price is the models: the NewEra UV flux of cool and active stars is uncertain
at the 0.1-0.3 mag level (chromospheres, missing line opacity), so the UV bands
are only kept for the hotter stars (Gaia BP-RP colour cuts below, ~T_eff > 5300 K
for NUV and > 6300 K for FUV) and carry a 0.10 mag systematic in the fit
(fit._phot_sys); their zero points are iterated per T_eff bin
(fit.uv_offsets_by_teff), because the model UV bias is T_eff dependent.

GUVcat also carries the SFD E(B-V) at each source ("E(B-V)" column): the
median over the field is recorded as ``sfd_ebv`` in the cached table's
attributes file (galex_sfd.json) for the reference line of the column figure.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from ..cache import Workspace
from .vizier import box_query
from .xmatch import match_to_gaia

MATCH_ARCSEC = 2.5       # GALEX PSF ~ 5"; AIS astrometry ~ 0.5"
SYS_FLOOR = 0.05         # catalogue-level floor (fit adds the model systematic on top)
MAX_ERR = 0.35
BPRP_MAX_NUV = 1.0       # ~T_eff > 5300 K
BPRP_MAX_FUV = 0.6       # ~T_eff > 6300 K
_COLS = ["RAJ2000", "DEJ2000", "FUVmag", "e_FUVmag", "NUVmag", "e_NUVmag",
         "Fafl", "Nafl", "Fexf", "Nexf", '"E(B-V)"']


def shape(d: pd.DataFrame, gaia: pd.DataFrame, ra0: float, dec0: float) -> pd.DataFrame:
    """Match GUVcat rows to Gaia and apply the quality/colour cuts."""
    m = match_to_gaia(gaia, d.RAJ2000.values, d.DEJ2000.values, ra0, dec0, MATCH_ARCSEC)
    p = d.iloc[m.icat.values]
    res = pd.DataFrame(dict(source_id=m.source_id.values, galex_sep=m.sep.values))
    bprp = gaia.set_index("source_id").bp_rp.reindex(res.source_id).values
    for band, col, afl, exf, bmax in (("NUV", "NUVmag", "Nafl", "Nexf", BPRP_MAX_NUV),
                                      ("FUV", "FUVmag", "Fafl", "Fexf", BPRP_MAX_FUV)):
        mag = p[col].values.astype(float)
        err = p["e_" + col].values.astype(float)
        good = (np.isfinite(mag) & np.isfinite(err) & (err < MAX_ERR)
                & (p[afl].values == 0) & (p[exf].values == 0)
                & np.isfinite(bprp) & (bprp < bmax))
        res[f"mag_GALEX_{band}"] = np.where(good, mag, np.nan)
        res[f"magerr_GALEX_{band}"] = np.where(good, np.sqrt(err ** 2 + SYS_FLOOR ** 2), np.nan)
    return res


def fetch(ws: Workspace, gaia: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    out = ws.path("galex_gaia.csv")
    if out.exists() and not force:
        return pd.read_csv(out)
    try:
        d = box_query("II/335/galex_ais", _COLS, ws.ra, ws.dec, ws.radius_arcmin / 60.0)
    except Exception as e:  # noqa: BLE001
        print(f"GALEX query failed ({e}); continuing without")
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    if not len(d):
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    ebv = d["E(B-V)"].astype(float)
    json.dump(dict(sfd_ebv=float(np.nanmedian(ebv)), n=int(np.isfinite(ebv).sum())),
              open(ws.path("galex_sfd.json"), "w"))
    res = shape(d, gaia, ws.ra, ws.dec)
    res.round(4).to_csv(out, index=False)
    n = {b: int(np.isfinite(res[f"mag_GALEX_{b}"]).sum()) for b in ("NUV", "FUV")}
    print(f"GALEX GUVcat_AIS: {len(res)} Gaia matches, usable {n} (SFD E(B-V) {np.nanmedian(ebv):.3f})")
    return res

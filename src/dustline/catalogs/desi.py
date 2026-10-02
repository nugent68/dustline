"""DESI DR1 Milky Way Survey stellar parameters (NOIRLab Data Lab table
``desi_dr1.mws``), joined to Gaia DR3 by source_id, as per-star template priors
for the XP fit.

The MWS value-added catalogue (Koposov et al. 2025, DESI DR1 stellar catalogue)
carries the RVSpecFit parameters T_eff, log g, [Fe/H], [alpha/Fe] with formal
errors for every DESI stellar spectrum; at high latitude the backup/bright
programs reach essentially every Gaia star to G ~ 19, so most Gaia XP stars in
a DESI field have one. We keep ``rr_spectype = STAR``, ``rvs_warn = 0`` and the
primary coadd, and add systematic floors to the formal errors (DR1 paper:
[Fe/H] good to ~0.1 dex at high S/N; T_eff scale differences between template
libraries are ~50-100 K).

[Fe/H] is used directly as the NewEra [M/H] prior (alpha enhancement ignored).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..cache import Workspace
from .datalab import box_query
from .local import MirrorError

TEFF_SYS = 1.0       # K: T_eff is LOCKED to the DESI label (grid point); see calib.TEFF_PRIOR_SIGMA
LOGG_SYS = 0.10      # dex
FEH_SYS = 0.10       # dex
_COLS = ["source_id", "target_ra", "target_dec", "teff", "teff_err", "logg", "logg_err",
         "feh", "feh_err", "alphafe", "snr_med", "survey", "program"]
_WHERE = [("rr_spectype", "==", "STAR"), ("rvs_warn", "==", 0), ("zcat_primary", "==", True)]
_EXTRA = "AND rr_spectype = 'STAR' AND rvs_warn = 0 AND zcat_primary = 't'"   # = render_adql(_WHERE)


def shape(d: pd.DataFrame) -> pd.DataFrame:
    """Reduce raw MWS rows to one row per Gaia source with the prior columns."""
    d = d[np.isfinite(d.source_id) & (d.source_id > 0)].copy()
    d["source_id"] = d.source_id.astype(np.int64)
    d = d.sort_values("snr_med", ascending=False).drop_duplicates("source_id")
    res = pd.DataFrame(dict(
        source_id=d.source_id.values,
        teff_spec=d.teff.values,
        teff_spec_err=np.full(len(d), TEFF_SYS),
        logg_spec=d.logg.values,
        logg_spec_err=np.sqrt(d.logg_err.values ** 2 + LOGG_SYS ** 2),
        feh_spec=d.feh.values,
        feh_spec_err=np.sqrt(d.feh_err.values ** 2 + FEH_SYS ** 2),
        alphafe_spec=d.alphafe.values,
        spec_snr=d.snr_med.values,
        spec_program=d.program.values,
    ))
    ok = np.isfinite(res.teff_spec) & np.isfinite(res.logg_spec) & np.isfinite(res.feh_spec)
    return res[ok].reset_index(drop=True)


ALPHA_TO_MH = True   # use [Fe/H] + log10(0.638 x 10^[a/Fe] + 0.362) (Salaris+93) as the [M/H] prior


def to_mh(res: pd.DataFrame) -> pd.DataFrame:
    """Turn the DESI [Fe/H] prior into total metallicity [M/H], the scaled-solar quantity of the
    NewEra grid: [M/H] = [Fe/H] + log10(0.638 x 10^[alpha/Fe] + 0.362) (Salaris et al. 1993).
    Against dustline's own NewEra fits with the DESI priors removed (benchmarks/desi_labels;
    20 fields, ~5,800 stars) DESI [Fe/H] reads 0.12-0.16 dex low and the alpha-corrected value
    agrees to +/- 0.02. The raw label is kept as feh_desi; applied when the priors are read, so
    cached desi_gaia.csv files (raw [Fe/H] + alphafe_spec) need no refetch."""
    if not ALPHA_TO_MH or "feh_spec" not in res or "feh_desi" in res or not len(res):
        return res
    res = res.copy()
    res["feh_desi"] = res["feh_spec"]
    a = res["alphafe_spec"].astype(float) if "alphafe_spec" in res else pd.Series(np.nan, index=res.index)
    corr = np.log10(0.638 * 10.0 ** a.clip(-0.2, 0.5) + 0.362)   # DESI [a/Fe] is noisy at low S/N
    res["feh_spec"] = res["feh_spec"] + np.where(np.isfinite(corr), corr, 0.0)
    return res


def fetch(ws: Workspace, gaia: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    """DESI MWS priors for the Gaia sources of the field (cached desi_gaia.csv)."""
    out = ws.path("desi_gaia.csv")
    if out.exists() and not force:
        return to_mh(pd.read_csv(out))
    try:
        d = box_query("desi_dr1.mws", _COLS, ws.ra, ws.dec, ws.radius_arcmin / 60.0,
                      ra_col="target_ra", dec_col="target_dec", where=_WHERE)
    except MirrorError:
        raise                             # a broken mirror must not be cached as "no data"
    except Exception as e:  # noqa: BLE001
        print(f"DESI query failed ({e}); continuing without spectroscopic priors")
        pd.DataFrame(dict(source_id=[])).to_csv(out, index=False)
        return pd.read_csv(out)
    res = shape(d) if len(d) else pd.DataFrame(dict(source_id=[]))
    if len(res):
        res = res[res.source_id.isin(gaia.source_id.values)].reset_index(drop=True)
    res.round(4).to_csv(out, index=False)
    print(f"DESI MWS: {len(res)} Gaia stars with spectroscopic parameters")
    return to_mh(res)

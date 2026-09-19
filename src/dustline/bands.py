"""Band assembly: merge the Gaia table, the selected survey photometry and any
user photometry into one per-star table whose photometric columns are named
``mag_<BandName>`` / ``magerr_<BandName>`` with BandName from the filter
registry (e.g. mag_PS1_g, mag_2MASS_Ks, mag_VISTA_J, mag_DECam_i).

The fit uses whatever bands are present; the *band list* of a run is the union
of the columns delivered here.  Spectroscopic template priors (DESI MWS
T_eff/log g/[Fe/H], columns ``teff_spec`` ...) ride on the same table when the
plan asks for them.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .cache import Workspace
from .catalogs import decals, decaps, desi, ps1, vvv
from .catalogs.footprint import SurveyPlan

TMASS_SYS = 0.03    # systematic added to 2MASS errors when used as a VVV fallback


def assemble(ws: Workspace, gaia: pd.DataFrame, plan: SurveyPlan,
             user_phot: pd.DataFrame | None = None,
             force: bool = False) -> tuple[pd.DataFrame, list[str]]:
    """Return (per-star table indexed like gaia, band names present)."""
    df = gaia.copy()

    # --- 2MASS: already joined in the Gaia query as mag_J/H/Ks -> 2MASS band names
    for b in ("J", "H", "Ks"):
        df = df.rename(columns={f"mag_{b}": f"mag_2MASS_{b}", f"magerr_{b}": f"magerr_2MASS_{b}"})

    # --- optical survey ---
    if plan.optical == "ps1":
        p = ps1.fetch(ws, gaia, force=force)
        df = df.merge(p.drop(columns=[c for c in ("ps1_sep", "ps1_nDet") if c in p]),
                      on="source_id", how="left")
    elif plan.optical == "decaps":
        p = decaps.fetch(ws, gaia, force=force)
        df = df.merge(p.drop(columns=[c for c in ("decaps_sep",) if c in p]),
                      on="source_id", how="left")
    elif plan.optical == "decals":
        p = decals.fetch(ws, gaia, force=force)
        df = df.merge(p.drop(columns=[c for c in ("decals_sep",) if c in p]),
                      on="source_id", how="left")

    # --- near-infrared: VVV preferred in its window, 2MASS fills the bright end ---
    if plan.nir == "vvv":
        v = vvv.fetch(ws, gaia, force=force)
        if len(v):
            df = df.merge(v.drop(columns=[c for c in ("vvv_sep",) if c in v]),
                          on="source_id", how="left")
            for b in ("J", "H", "Ks"):
                use2m = (~np.isfinite(df.get(f"mag_VISTA_{b}", np.nan))
                         & np.isfinite(df[f"mag_2MASS_{b}"]))
                df.loc[use2m, f"mag_VISTA_{b}"] = df.loc[use2m, f"mag_2MASS_{b}"]
                df.loc[use2m, f"magerr_VISTA_{b}"] = np.sqrt(
                    df.loc[use2m, f"magerr_2MASS_{b}"] ** 2 + TMASS_SYS ** 2)

    # --- spectroscopic template priors (DESI DR1 MWS), joined by Gaia source_id ---
    if plan.spectro == "desi":
        sp = desi.fetch(ws, gaia, force=force)
        if len(sp):
            df = df.merge(sp, on="source_id", how="left")

    # --- user photometry overrides/augments everything ---
    if user_phot is not None:
        overlap = [c for c in user_phot.columns
                   if c.startswith("mag") and c in df.columns]
        df = df.merge(user_phot, on="source_id", how="left", suffixes=("_survey", ""))
        for c in overlap:   # user value wins; survey fills gaps
            if c + "_survey" in df:
                df[c] = df[c].fillna(df[c + "_survey"])

    bands = sorted({c[4:] for c in df.columns if c.startswith("mag_")
                    and not c.startswith("mag_survey") and f"magerr_{c[4:]}" in df.columns})
    # drop the plain 2MASS duplicates when VVV took over the NIR slots
    if plan.nir == "vvv" and any(b.startswith("VISTA_") for b in bands):
        bands = [b for b in bands if not b.startswith("2MASS_")]
    counts = {b: int(np.isfinite(df[f"mag_{b}"]).sum()) for b in bands}
    print(f"bands assembled: {counts}")
    return df, bands

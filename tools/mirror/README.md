# dustline survey mirror (NERSC)

Root: `/global/cfs/projectdirs/newera/surveys` (project `newera`). The package reads it when
`DUSTLINE_SURVEYS` points here. Otherwise every catalog comes from the network, exactly as before.

```bash
export DUSTLINE_SURVEYS=/global/cfs/projectdirs/newera/surveys
export DUSTLINE_CACHE_DIR=/global/cfs/projectdirs/newera/dustline_cache   # model assets + workspaces + grids
export DUSTLINE_MIRROR=auto      # auto (default) | off (always network) | only (never network: fail instead)
source /global/cfs/projectdirs/newera/surveys/_env/bin/activate
dustline mirror-info             # what resolves locally
dustline run 267.86642 -33.13517 --radius 5 --plx-inflate 1.7
```

Each survey directory has a `manifest.json`. It is keyed by the network table name the package
queries and records:
- the partitioning;
- the column map to the network column names;
- the null conventions to undo;
- provenance;
- the validation against the network service.

The package uses a table only when its `status` is `complete`. A table stays `partial` until
`tools/mirror/check_table.py` (or, for the Gaia tables, the checks listed below) has compared it
with the service.

## Contents

| network table | mirror | rows | status / validation |
|---|---|---|---|
| `gaiadr3.gaia_source` | cosmo `gaia/dr3/healpix` in place (nside 32, sorted by source_id) | 1.81 G | complete. ESA cone for 0095 identical; nulls are stored as 0 (undone by the manifest), booleans as text |
| `gaiadr3.xp_continuous_mean_spectrum` | `gaia_dr3/xp_continuous/` ESA CDN, one FITS per CDN file (source_id range) | 219,197,643 | complete. MD5 plus read-back checks; DataLink rows bit-identical; 0095 `xp_sampled` identical |
| `gaiadr3.tmass_psc_xsc_best_neighbour` | `gaia_dr3/tmass_best_neighbour/hp32` (EDR3 table on the CDN) | 469,051,627 | complete. 0095 join identical to ESA |
| `gaiadr1.tmass_original_valid`, `twomass.psc` | cosmo `2mass/healpix` in place (bbox index; `ph_qual` rebuilt from code columns) | 470,992,970 | complete. 0095 magnitudes, errors and quality identical |
| `sdss_dr17.apogee2_allstar` | `apogee_dr17/allstar/hp32` (from cosmo allStar DR17) | 733,900 | complete. Data Lab: same rows; values to its rounding; PostgreSQL NaN semantics |
| `desi_dr1.mws` | `desi_dr1/mws/hp32` (from the DR1 MWS VAC) | 6,372,607 | complete. Data Lab: same rows and values; 2/1281 `source_id` swapped by Data Lab |
| `II/335/galex_ais` | `galex/guvcat_ais/hp32` (MAST HLSP GUVcat AIS) | 82,992,062 | **partial**. About 1 % of sources in tile overlaps are a different measurement than VizieR II/335 |
| `allwise.source` | `allwise/source/hp32` (from cosmo allwise-catalog) | 747,634,026 | complete. Data Lab: identical rows and values |
| `decaps_dr2.stellar_inference` | `zucker25/stellar_inference/hp32` (Harvard Dataverse K88GFI) | 709,129,917 | complete. Data Lab: identical rows; values to rounding; 0095 `vvv_gaia.csv` identical |
| `vvv.virac2` | — (ESO tap_cat only) | 545 M | not mirrored: one dense bulge pixel took >1 h through ESO TAP; ESO's phase-3 files hold only the light curves (11.2 TB) |
| `decaps_dr2.object` | — | 3.3 G | **not available**: the DECaPS2 dbfits bucket returns 403 |
| PS1 DR2 | — | — | not mirrored (no bulk dump); MAST |
| `ls_dr10.tractor` | — | — | not done yet (cosmo `legacysurvey/dr10` sweeps) |

## 3D dust map: Edenhofer et al. (2023)

`edenhofer_2023/mean_and_std_healpix.fits`: the main map (mean and std of the differential
extinction density, HEALPix, 69 pc to 1.25 kpc), added 2026-10-01. The directory follows the
`dustmaps` layout (`<data_dir>/edenhofer_2023/<file>`), so any `dustmaps` user can read it in place:

```python
from dustmaps.config import config
config["data_dir"] = "/global/cfs/projectdirs/newera/surveys"
from dustmaps.edenhofer2023 import Edenhofer2023Query
q = Edenhofer2023Query(integrated=True)        # integrated E (ZGR23 units; A_V ~ 2.8 E)
```

In dustline it is used only by `tools/build_template_corrections.py` (the `fit` stage, via
`calib.map_extinction`). That code looks in `$DUSTLINE_CACHE_DIR/dustmaps`, so point it here with
`ln -s /global/cfs/projectdirs/newera/surveys/edenhofer_2023 $DUSTLINE_CACHE_DIR/dustmaps/edenhofer_2023`.
`dustmaps` itself is not installed in `_env` (`pip install dustmaps` if needed). This is a map, not a
catalog, so it has no `manifest.json` and is not read by the local-mirror backend.

## Layout

- `<survey>/manifest.json` plus data: position tables are HEALPix nside-32 NESTED files, with rows sorted by the nside-4096 pixel (`_HPX`, or the Gaia `source_id >> 35`).
- `_tools/mirror/`: ingest scripts (copies of the package's `tools/mirror/`).
- `_tools/pkgsrc`, `_tools/dustline-pkg`: the package.
- `_env/`: venv (`--system-site-packages` on NERSC python).
- `_ledger/`: JSONL records of finished work units, plus claim locks.
- `_logs/`, `_staging/`: logs and transient downloads.

Every ingest script is resumable: rerun it and it skips the units in its ledger.

## End-to-end check (2026-09-30)

OGLE-2017-BLG-0095, 5′, run on a Perlmutter node with the mirror. Gaia, XP, 2MASS, APOGEE and VVV
came from the mirror; DECaPS photometry came from Data Lab.

- The data stage took about 2 minutes. Over the network it takes about 2 hours, dominated by the XP fetch.
- `gaia.csv`, `xp_continuous_raw.csv` and `vvv_gaia.csv` are identical to the network-built files, and `xp_sampled.npz` is bit-identical.
- R_V = 3.17 ± 0.25 from 333 stars, the same as the reference run, and the A_V(D) run agrees to 1e-4.
- Per-star differences are floating-point level (machine BLAS): a median ΔA_V of 5e-5, with 2 of 2,995 stars switching best-fit node.

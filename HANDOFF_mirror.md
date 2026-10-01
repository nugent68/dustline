# HANDOFF: dustline survey mirror at NERSC: what data is where

Status as of 2026-09-30. Branch `local-mirror` (commits 54b0dd1, 886bb53, 3a6646b, plus this file). It is pushed to GitHub and is
not merged into `main`.

## 1. What this is

The dustline package normally fetches its catalogs over the network: ESA (Gaia, XP spectra), NOIRLab Data Lab,
MAST and VizieR. The bulk-downloadable ones are now mirrored, or referenced in place, at NERSC. The package
reads them when `DUSTLINE_SURVEYS` is set. With it unset, nothing changes and every catalog comes from the network.

```bash
source /global/cfs/projectdirs/newera/surveys/_env/bin/activate        # package installed (editable) from _tools/dustline-pkg
export DUSTLINE_SURVEYS=/global/cfs/projectdirs/newera/surveys
export DUSTLINE_CACHE_DIR=/global/cfs/projectdirs/newera/dustline_cache
export DUSTLINE_MIRROR=auto     # auto: mirror where complete, else network | off: always network | only: never network (raise)
dustline mirror-info            # shows what resolves locally
```

Effect on OGLE-2017-BLG-0095 (5′): the data stage takes about 2 min instead of about 2 h. The products are identical to the network
run, and R_V = 3.17 ± 0.25 is reproduced (see §6).

## 2. Mirrored tables: `/global/cfs/projectdirs/newera/surveys` (project `newera`, quota 20.5 TB)

The package looks each table up by the **network table name** in `<survey>/manifest.json`. It uses a table only when
the manifest's `status` is `complete`.

| network table (what the package queries) | location | files | size | rows | status / validation |
|---|---|---|---|---|---|
| `gaiadr3.xp_continuous_mean_spectrum` | `gaia_dr3/xp_continuous/XpContinuousMeanSpectrum_AAAAAA-BBBBBB.fits` (one per ESA CDN file; each is a HEALPix level-8 `source_id` range) | 3,386 | 2.7 TB | 219,197,643 (all DR3) | complete. MD5 checked against the CDN; DataLink rows bit-identical |
| `gaiadr3.tmass_psc_xsc_best_neighbour` | `gaia_dr3/tmass_best_neighbour/hp32/NNNNN.fits` (by `source_id >> 49`) | 12,288 | 13 GB | 469,051,627 | complete. 0095 join identical to ESA |
| `gaiadr3.gaia_source` | **in place:** `/global/cfs/cdirs/cosmo/data/gaia/dr3/healpix/healpix-NNNNN.fits` (manifest in `gaia_dr3/`) | 12,288 | (cosmo) | 1.81 G | complete. 0095 cone identical to ESA |
| `gaiadr1.tmass_original_valid`, `twomass.psc` | **in place:** `/global/cfs/cdirs/cosmo/data/2mass/healpix/2mass_hp000–971.fits` (manifest + bbox index in `twomass/`) | 972 | (cosmo, 102 GB) | 470,992,970 | complete. 0095 2MASS identical |
| `decaps_dr2.stellar_inference` (Zucker+25; our VVV JHKs) | `zucker25/stellar_inference/hp32/NNNNN.fits` | 896 | 156 GB | 709,129,917 | complete. Matches Data Lab on 4 fields; 0095 `vvv_gaia.csv` identical |
| `allwise.source` | `allwise/source/hp32/NNNNN.fits` (columns extracted from cosmo `wise/allwise-catalog`) | 12,288 | 33 GB | 747,634,026 | complete. Identical to Data Lab |
| `sdss_dr17.apogee2_allstar` | `apogee_dr17/allstar/hp32/` (from cosmo `sdss/dr17/.../allStar-dr17-synspec_rev1.fits`) | 5,708 | 163 MB | 733,900 | complete. Matches Data Lab |
| `desi_dr1.mws` | `desi_dr1/mws/hp32/` (from `/global/cfs/cdirs/desi/public/dr1/vac/dr1/mws/iron/v1.0/mwsall-pix-iron.fits`) | 4,928 | 888 MB | 6,372,607 | complete. Matches Data Lab except 2 of 1,281 `source_id`s, which Data Lab mis-joined; the mirror follows the VAC |
| `II/335/galex_ais` (GUVcat) | `galex/guvcat_ais/hp32/` (from the MAST HLSP GUVcat AIS CSVs) | 10,068 | 4.8 GB | 82,992,062 | **partial: not used.** About 1 % of sources in tile overlaps are a different measurement than VizieR II/335. Undecided whether to accept it |

Conventions:
- **Position tables:** HEALPix nside 32, NESTED, one file per pixel, rows sorted by the nside-4096 pixel (`_HPX` column).
- **Gaia tables:** sorted by `source_id`, whose `>> 35` is the level-12 pixel. That matches the star's position pixel for 99.4 % of stars, so cones pad by 2′.
- **Column names:** the files use the network (Data Lab / VizieR / ESA) names, or the manifest maps them.

### Quirks of the NERSC copies (handled in the manifests; know them if you read the files directly)
- **cosmo Gaia DR3** stores **every NULL float as 0.0**: parallax and proper motions for 2-parameter solutions, BP/RP magnitudes, GSP-Phot outputs, and more. Booleans are stored as the text `'True'`/`'False'`. The manifest marks float columns `null_zero` and the booleans `bool_text`. This was verified against ESA for all 27,077 stars of the 0095 cone, with zero mismatches.
- **cosmo 2MASS** files are 972 cells of about 10°×8.5°, not HEALPix. The manifest stores a per-file RA/Dec bounding box. Quality flags are byte codes 0–7, mapping to X, U, F, E, A, B, C, D; `ph_qual` is rebuilt from them. X=0 is inferred; the other codes were observed.
- **Data Lab-derived tables** (APOGEE, DESI, AllWISE, stellar_inference) carry `"nan_semantics": "postgres"`. On Data Lab, `x > 0` *keeps* NaN rows, and the local filters reproduce that.
- APOGEE: a missing `gaiaedr3_source_id` is 0 in allStar and −2⁶³ on Data Lab; the manifest replaces it.

## 3. Not mirrored (still network) and open items
| table | why | what would fix it |
|---|---|---|
| `decaps_dr2.object` (DECaPS2 grizY, used by all bulge/plane fields) | `decaps.rc.fas.harvard.edu/.../DR2_REDUX/DATABASE/dbfits/*.fits` returns **403** (only the datamodel example downloads; `decaps.skymaps.info` https is dead). NERSC has only DECaPS **DR1 coadds** (`cosmo/data/decaps/dr1`) and coadd work dirs (`cosmo/work/decaps2_final`), no DR2 catalog | Ask Schlafly / Saydjari for access, or export the 19 needed columns from Data Lab in sky chunks (about 3.3 G rows, ~260 GB). Data Lab queries take seconds, so this is a dependency, not a speed problem |
| `vvv.virac2` (optional alternative VVV source) | ESO `tap_cat` `VVVX_VIRAC_V2_SOURCES` (545 M rows): one dense bulge nside-32 pixel ran over 1 h and was stopped. ESO phase-3 files for this release are the **light curves** only (22,583 files, 11.2 TB) | A bulk export from ESO or the VIRAC team. `tools/mirror/virac2.py` (async TAP per pixel, `EXECUTIONDURATION=3600`, `MAXREC` 15 M) is ready if throughput improves |
| PS1 DR2 | No bulk dump. cosmo `work/ps1/cats/chunks-qz-star-v3` lacks QfPerfect and per-band errors and covers only −22<Dec<32, \|b\|>15 | Stays on MAST (small high-latitude cones) |
| `ls_dr10.tractor` (Dec ≤ −30, \|b\| > 15 fields only) | Not done | Extract PSF-only rows (mag/snr/fracflux) from cosmo `legacysurvey/dr10/south` sweeps, as `extract.py` does for others |
| GUVcat | See §2 | Accept the HLSP version (set status `complete` in `galex/manifest.json`) or keep VizieR |

## 4. Other NERSC locations
| what | where |
|---|---|
| Package (editable install used by `_env`) | `surveys/_tools/dustline-pkg/`. Source-only copy for the ingest scripts: `surveys/_tools/pkgsrc/` |
| Ingest scripts | `surveys/_tools/mirror/` (copies of the repo's `tools/mirror/`) |
| Python env | `surveys/_env` (venv, `--system-site-packages` on NERSC `module load python`, i.e. python 3.13; has gaiaxpy 2.1.4) |
| Work ledgers / claim locks | `surveys/_ledger/*.jsonl`, `surveys/_ledger/locks/`. **Rerunning any ingest script skips the units recorded here** |
| Logs | `surveys/_logs/` |
| Validation references (small) | `surveys/_staging/`: `cmp_0095/` (the Mac network-built 0095 products), `esa_0095_gaia.csv`, comparison scripts |
| dustline cache: model assets, grid caches, workspaces | `/global/cfs/projectdirs/newera/dustline_cache/` (NewEra uvir/full caches, MIST, template corrections, `xp_grid_g23_d5ab7f40b8.npz`, `sightlines/ra+0267.86642_dec-033.13517_r5_b51f4ceb/`: the NERSC 0095 run) |
| Provenance files of deleted raw downloads | `zucker25/dataverse_files.json`, `gaia_dr3/xp_continuous/_MD5SUM.txt`, `gaia_dr3/tmass_best_neighbour/_MD5SUM.txt`, `galex/guvcat_hlsp_readme.txt` |

The raw downloads were deleted after verification: XP CSVs, Zucker+25 batches, GUVcat CSVs and the 2MASS cross-match CSVs.
All can be re-fetched by the scripts below.

## 5. Re-creating or extending (all in `tools/mirror/`, resumable)
| table | command(s) | source |
|---|---|---|
| Gaia XP | `gaia_xp.py list`; `gaia_xp.py run --task-id i --ntasks n` (or `slurm/gaia_xp_node.sh` on a node: 32 workers; a second node with `XP_TAG=_b XP_ARGS=--reverse`); `gaia_xp.py manifest` | ESA CDN, S3 listing at `https://gaia.eu-1.cdn77-storage.com/?prefix=Gaia/gdr3/Spectroscopy/xp_continuous_mean_spectrum/` (~150 MB/s per stream from NERSC) |
| 2MASS best neighbour | `gaia_tmass_bn.py download` (DTN), `convert`, `manifest` | CDN `Gaia/gedr3/cross_match/tmasspscxsc_best_neighbour/` (DR3 reuses EDR3) |
| gaia_source, 2MASS manifests | `reference_manifests.py gaia_source`, `twomass`, `twomass-qual codes.json` | cosmo copies |
| Zucker+25 | `dataverse_download.py doi:10.7910/DVN/K88GFI <staging>` (DTN; needs a User-Agent, already set), `zucker25.py convert` (two streaming passes; one pass on one node ran out of memory), `manifest` | Harvard Dataverse |
| APOGEE / DESI / AllWISE | `extract.py apogee|desi|allwise` | cosmo / DESI public copies |
| GUVcat | `guvcat.py download|convert|manifest` | MAST HLSP |
| Validation | `check_table.py <table> [--mark]`: same box queried over the network and locally; `--mark` sets `complete` on a pass | |

Practicalities:
- Use the **DTNs** (`dtn01.nersc.gov`) for downloads, running NERSC python by full path: `/global/common/software/nersc/pe/conda-envs/26.8.1/python-3.13/nersc-python/bin/python`. They have no `$PSCRATCH`.
- Use **`salloc -q interactive -C cpu -N 1 -t 04:00:00`** for conversions. The batch queues were days deep.
- **Account m2218** was used, about 5 node-hours. The default account **m5071 has no balance**.

## 6. Validation run (0095)
- NERSC workspace: `dustline_cache/sightlines/ra+0267.86642_dec-033.13517_r5_b51f4ceb/`. Mac network workspace: `~/.cache/dustline/sightlines/ra+0267.86642_dec-033.13517_r5_b51f4ceb/`.
- Identical between the two: `gaia.csv` (27,077 stars, all columns), `xp_continuous_raw.csv` (2,995 ids), `vvv_gaia.csv` (11,609).
- `xp_sampled.npz` is bit-identical.
- R_V = 3.17 ± 0.25 (333 stars); the A_V(D) run agrees to 1e-4.
- Per-star fits differ at floating-point level only: median ΔA_V 5e-5; 2 of 2,995 stars change best-fit grid node.

## 7. Package code (branch `local-mirror`)
- `src/dustline/catalogs/local.py`: manifest-driven reader. It handles memory-mapped FITS with O(log n) bisection, header-only opens, TZERO/TSCAL, null and boolean conventions, array-element columns, PostgreSQL NaN semantics, and the service's box or cone geometry, and ends with a CSV round trip so dtypes match. Also `gaia_cone` and the XP writer, which emits DataLink-format CSV.
- `src/dustline/catalogs/_healpix.py`: pure-numpy nested HEALPix. No new dependency.
- **Hooks:** `datalab`, `vizier` and `gaia` try the mirror first. Fetchers pass structured `where` filters, which are rendered back to the identical ADQL for the network. `MirrorError` is never cached as an empty result.
- `dustline mirror-info`; `tests/test_mirror.py` (16 tests).

## 8. Mac-side locations (for reference)
- Worktree for this branch: `~/claude/dustline-pkg-mirror`. The main checkout `~/claude/dustline-pkg` is on `field-curve` and another session works there.
- Mac workspaces: `~/.cache/dustline/sightlines/`. 0095 is in `…r5_b51f4ceb` (official v0.8.0 run) and `…r5_940b3867` (field-curve experiments, including `xp_stars_cut940_*`). COSMOS is in `ra+0150.12000_dec+002.21000_r30_a7380f18`.

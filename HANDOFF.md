# dustline — HANDOFF (state at v0.7.0, 2026-09-21)

For a fresh agent session in this repo (`/Users/nugent/claude/dustline-pkg`, GitHub
`nugent68/dustline`, public, tag `v0.7.0`, HEAD `0d28677`). Self-contained; the science
history is in `docs/method.md` and `CHANGELOG.md`, the research prototype (per-sightline
scripts, the OB170095 results this package must reproduce) is `~/claude/dustline`
(`HANDOFF_dustline.md` there).

## 1. What the package does

`dustline.Sightline(ra, dec, radius_arcmin=…).run()` measures, from Gaia DR3 XP spectra ×
PHOENIX NewEra templates × survey photometry:
- **law mode** (|b| < 30°): R_V (G23 family) from the well-reddened parallax stars, the
  band ratios A_X/A_V, the 3D run A_X(D) as a running median, and on bulge sightlines the
  red-clump bridge to the bulge column;
- **column mode** (|b| > 30°): R_V = 3.1 assumed, the total foreground A_V from stars
  behind the dust (an SFD/Planck cross-check; COSMOS example).

Per-star fit: brute-force grid (T_eff, log g, [M/H]) × A_V (136 values, 0.01 steps below
0.5) × R_V (2.3–5.55, 0.25) with the flux scale profiled, MIST radius prior from the
parallax, DESI DR1 / APOGEE DR17 log g–[Fe/H] priors where available, empirical template
corrections, and (v0.7.0) a fine local (R_V, A_V) refinement around the coarse optimum.
Photometric zero points are iterated (NIR frozen after pass 1; GALEX per T_eff bin).

## 2. Environment and caches

- venv: `.venv` (python 3.14; `source .venv/bin/activate`; `dustline` CLI on PATH). Tests:
  `python -m pytest -q -m "not network"` (49 pass; the two regression tests use the cached
  bulge and COSMOS workspaces).
- Assets (fetched on first use into `~/.cache/dustline/`, sha-verified, registry in
  `src/dustline/assets.py`): `newera_full_cache.npz` + `mist_v1.2_basic.npz` (v0.1.0),
  `newera_uvir_cache.npz` (v0.4.0, 214 MB, 900 Å–6 µm, the default), `template_corrections.npz`
  (v0.7.0, tag `tc26438c9b`, includes the R_V closure table).
- Grid products: `~/.cache/dustline/xp_grid_g23_<key>.npz`, keyed by bands + law + [M/H]
  axis + model cache + corrections tag (+ `b2` since the band-correction fix). A new key
  rebuilds in ~7 min (4366 models) or ~1 min (solar-only, 730 models).
- Workspaces: `~/.cache/dustline/sightlines/ra+…_dec…_r<radius>_<confighash>/` holding
  `gaia.csv`, `xp_continuous_raw.csv`, `xp_sampled.npz` (the expensive part: 1.8 s/star
  from the Gaia archive, resumable), survey csvs, `xp_stars.csv` (per-star fits),
  `phot_offsets.json`, `result.json`, `dust_run_av.csv`. The config hash includes the
  corrections tag, so a new tag → a new directory that **seeds gaia/XP/survey files from a
  sibling of the same position** automatically.
- **Never delete a workspace** to force a rerun (an 80-min XP refetch was lost that way).
  Use `dustline run … --refit` (redo fits + ensemble from cached data) or delete only
  `result.json` + `dust_run_av.csv` to redo the ensemble stage (seconds).
- Timings: fits ~0.5 s/star (4366-model grid) per pass, 3 passes (+1 with GALEX); a 2,600-star
  field ≈ 1–1.5 h alone, ~2× that with two fields in parallel (10 cores, mostly single-thread).
- Gaia archive has 408 outages; `gaia.by_ids`/`cone` fall back to NOIRLab Data Lab.
  PS1 (MAST) is fetched by recursive cone tiling with a saturation cut. DESI (`desi_dr1.mws`),
  APOGEE (`sdss_dr17.apogee2_allstar`), AllWISE via Data Lab; GUVcat via VizieR.
- Gotchas: a long-running process that imported a module you then edit will crash at the
  stage that uses the edit (rerun the CLI; fits are cached). The Claude sandbox blocks
  `rm -r` on the cache (stray config-only workspaces from pytest at intermediate tags can be
  left; harmless). Background jobs: use `nohup … &` with absolute paths.
- NERSC (for rebuilding the NewEra cache only): `salloc -q interactive --no-shell` +
  `srun --jobid`; `tools/build_newera_uvir_cache.py`.

## 3. Code map (`src/dustline/`)

| module | what |
|---|---|
| `api.py` | `Sightline` (options: radius, photometry, prefer_deep, spectro_priors, freeze_offsets, uv, mir, min_av, desi_teff), `run(force, refit)`, `ExtinctionResult` (`.rv`, `.extinction(band)`, `.column`, `.plots()`, `.save()`) |
| `cli.py` | `dustline run RA DEC [--radius] [--min-av] [--no-spectro] [--freeze-offsets] [--no-uv] [--no-mir] [--desi-teff] [--ref-av] [-o] [--plots DIR] [--force \| --refit]`, `dustline fetch-assets` |
| `fit.py` | `build_grid`, `_fit_batch` (χ² cubes), `spec_prior_chi2`, `_summarize`, `_refine`/`_refine_star` (sub-grid), `fit_stars`, `measure_offsets`, `uv_offsets_by_teff`, `run_fit_with_offsets` (passes; `teff_from_colour` in column mode) |
| `calib.py` | template corrections: calibrator selection (DESI dwarfs D<250 pc; APOGEE giants D<1.2 kpc), `teff_from_bprp` (Mamajek locus + [Fe/H] term), `map_extinction` (Edenhofer+2023 via dustmaps), `deredden`, `ratio_spectrum`, `aggregate` (bins by the **best-fit model**, unsmoothed), `load`/`tag`/`apply`, `rv_closure` |
| `ensemble.py` | `law_sample`, `measure_law` (closure-corrected; `rv_raw`), `default_law`, `dust_run`, `bridge_to_clump`, `foreground_column`, `band_ratio_at_rv` |
| `clump.py` | red-clump anchor in (J−Ks, Ks) |
| `models.py`, `extinction.py`, `filters.py`, `bands.py`, `plotting.py`, `cache.py` | NewEra cache + XP LSF + synthetic photometry; G23/F99 curves; 26 packaged filter curves (`data/filters`); band assembly; figures; `Workspace` |
| `catalogs/` | `footprint.plan` (survey selection by position), `gaia`, `ps1`, `decaps`, `decals`, `vvv`, `desi`, `apogee`, `galex`, `wise`, `datalab`, `vizier`, `xmatch` |

Tools: `tools/build_template_corrections.py` (stages `pull`/`fit`/`build` for `dwarfs` /
`giants`; calibration workspaces `ra+0000.00000_dec+000.00000_r0_8659c9cc` (dwarfs) and
`_83240b75` (giants)), `tools/inject_reddening.py` (closure test; `closure` stage writes
`rv_closure_k` into the corrections file), `tools/run_ob240669.py` (regenerates the bulge
example from the research repo's XP cache), `tools/fetch_svo_filters.py`.

## 4. v0.7.0 in one paragraph (details: docs/method.md "Reddening-injection closure")

The injection test (calibrators reddened with a known law, refitted as a field) measured the
T_eff–A_V–R_V coupling (+0.45 R_V, +0.10 A_V per +100 K, symmetric → the ensemble median is
unbiased, but **never split a field's R_V by fitted T_eff**) and fixed: (1) `build_grid`
cancelled the per-band template corrections at every A_V ≠ 0; (2) the corrections table was
smoothed at σ 20 nm and binned by the star's [Fe/H]/log g instead of the model's; (3) high-S/N
posteriors quantised to the grid. A linear 1/A_V residual of the cool templates remains
(dwarfs < 5100 K, giants < 4500 K read R_V 2.7–2.9 for 3.05) and is applied per star as the
closure table. Budget after it: ±0.05 closure, ±0.05 NIR zero points.

Results on this stack (closure-corrected, raw in brackets):

| field | mode | result | N | example |
|---|---|---|---|---|
| OB240669 (l −5°, b −3°, bulge) | law + clump | R_V 2.80 ± 0.54 [2.73]; clump A_V 3.36 at 6.1 kpc | 883 | `examples/ob240669` |
| ZTF19abqmpti (l 11°, b +24°) | law | R_V 3.48 ± 0.48 [3.37]; column 1.07 | 2,228 | `examples/ztf19abqmpti` |
| ZTF20abgaovd (l 16°, b +27°) | law | R_V 3.50 ± 0.84 [3.32]; column 0.57 | 1,394 | `examples/ztf20abgaovd` |
| COSMOS (b +42°) | column | F/G 0.029 ± 0.004; K/M 0.039; SFD 0.059 | 132 | `examples/cosmos` |

## 5. Open issues (in priority order)

1. **K-dwarf T_eff offset**: 4500–5000 K dwarfs come out ~67 K too cool with A_V −0.10 in
   the injection test (independent of the parallax weight), and ZTF20's APOGEE dwarfs read
   −108 K. Candidates: the Mamajek-locus T_eff of the calibration vs the T_eff the XP shape
   + MIST radius prefer at the node (the corrected node templates are then "labelled" hot),
   the [M/H] bin mixing (0.5 dex model steps vs the −0.4 bin edge), the log g prior (DESI
   dwarfs fit at log g 5.5, the grid edge). Tools: `inject_reddening.py … locked` vs `free`,
   the per-node ΔT_eff of `inject_dwarfs_av1.00_rv3.05_free_plx10.csv` in the dwarf
   calibration workspace, ZTF20's `apogee_check`.
2. Closure calibrators are solar-neighbourhood ([Fe/H] −0.5..+0.3, [M/H]-dependent
   residual); bulge/thick-disk fields extrapolate.
3. Free giant fits sit 50–120 K below ASPCAP in reddened fields (no luminosity lever).
4. Column-mode absolute zero point ±0.03 (XP ~0.03 mag too red in g−z for faint red stars).
5. GALEX adds nothing at A_V < 0.1 and drops out where R_V is measurable; WISE is neutral.

## 6. Task for the next agent: OGLE-2017-BLG-0095 on the v0.7.0 stack

**Why**: 0095 is the sightline the method was built for (research prototype, 2026-09-17:
`~/claude/dustline/results/ob170095/`, `HANDOFF_dustline.md` §0). Its deliverable
(`P_SEDdust_XP`, the source-distance prior sent to Natasha) rests on a 2500–25000 Å cache,
solar-[M/H] templates, no corrections, no closure, the 0.25/0.1 grid — everything v0.7.0
changed. A package run is both a regression of the package on the DECaPS/VVV/clump path
(only exercised by the bulge example so far) and an updated deliverable.

**Target**: RA 267.86642, Dec −33.13517 (l 357.1°, b −4.0°). `footprint.plan` gives
`optical=decaps, nir=vvv, in_bulge_window=True, apogee=True, spectro=none` (no DESI here).
There is **no XP cache for 0095 on this machine** (the research repo's `data/` only holds
OB240669); the Gaia fetch is the cost: the research run used 5′ (2,995 XP stars ≈ 1.5 h);
9′ as in the bulge example ≈ 9,000 stars ≈ 4.5 h (resumable). Recommend **5′ first**, so
the star-by-star comparison with the prototype is direct, then 9′ if the clump window needs it.

```bash
source .venv/bin/activate
nohup dustline run 267.86642 -33.13517 --radius 5 -o examples/ob170095/extinction_I.csv \
      --plots examples/ob170095 > /tmp/ob170095.log 2>&1 &
```
(the CLI reports progress; `--plots` writes `law_rv.png`, `extinction_run_I.png`; copy the
workspace's `result.json` to `examples/ob170095/law.json`. If the Gaia archive 408s, the
Data Lab fallback engages by itself; if the DECaPS/VVV fetch fails, `catalogs/decaps.py` /
`vvv.py` are the places to look — both were last exercised in v0.2 on OB240669.)

**Compare against** (`~/claude/dustline/results/ob170095/`, format notes below):
- law: R_V 3.15 (MAD 0.23, 195 stars with A_V ≥ 2); ratios for a 6000 K source at
  R_V 3.15: A_g/A_i 1.885, A_r 1.335, A_z 0.767, A_Y 0.663, A_J 0.442, A_H 0.282,
  A_Ks 0.179 (`xp_law.json`). Expect the v0.7 R_V to differ by up to ~0.1–0.2 (band
  corrections at A_V > 0, closure +0.05–0.1 for the cool giants).
- run: A_i 0.59 / 1.19 / 1.33 / 1.39 / 1.81 at 0.5–1 / 1–1.5 / 1.5–2 / 2–3 / 3–5 kpc
  (`xp_dust_run.csv`, columns `D_kpc, AV_target(=A_i), …`), bridged to the clump column
  **A_i 1.90 ± 0.06 ± 0.19 at 8.0 kpc** (`clump_anchor.json`: 650 stars within 3′,
  (i−Ks, Ks) = (3.41, 13.24)). The package's `extinction_I.csv` is `D_kpc, A_med, A_16,
  A_84, n_stars, bridged, band, ratio_av` in Cousins I (`--filter DECam_i` or
  `res.extinction("DECam_i")` gives DECam i directly). The package's clump is found in
  (J−Ks, Ks) with VVV; the prototype used (i−Ks) — check `result.json["clump_anchor"]`
  (`E_JK`, `AV_column`, `D_RC`) against 1.90 / 8.0 kpc.
- zero points: DECaPS g −0.066, r −0.092, i −0.036, z −0.034, Y −0.041; VVV J −0.007,
  H +0.002, Ks +0.042 (`phot_offsets.json`, obs − synthetic). The package iterates its
  own (`phot_offsets.json` in the workspace; NIR frozen after pass 1).
- APOGEE: 16 giants within 15′, T_eff(XP) − T_eff(ASPCAP) = +49 K for [M/H] > −0.3
  (`apogee_check.txt`); the package now reports `result.json["apogee_check"]` (and has
  [M/H] −2..+0.5 templates, so the metal-poor giants should no longer be +386 K).
- residuals: the prototype's stacked (data − model)/model wiggles (+7 % at 640 nm, −5 % at
  500–540 and 770 nm) should be much reduced by the corrected templates — worth one
  figure (`xp_residuals.txt` has the old numbers).

**Deliverable**: the source-distance prior is produced by the declens snapshot
`~/claude/dustline/ref/ob170095/scripts/sed_distance_prior.py --profile <run.csv>
--profile-unit ai --laws sf11,f99rv25,measured --tag _xp`, which reads the prototype's
`xp_dust_run.csv` profile format (`D_kpc, AV_target, AV_nb_min, AV_nb_med, AV_nb_max, N`
in A_i) and `xp_law.json` (`ratios_source_g23` for the 6000 K source). Write a small
converter from the package outputs (`extinction_I.csv` in DECam i → the profile columns;
`law.json["rv"]` + `ensemble.band_ratio_at_rv(band, rv, teff=6000, logg=4.0)` → the
source ratios) rather than editing the declens script; then rerun the prior and compare
with **5.42 kpc, 68 % [3.95, 6.66], P(D ≥ 6) 0.35**. Port nothing back to declens without
the user's say-so.

**Also worth doing on the way**: pin an `examples/ob170095/README.md` in the style of the
others, add a `tests/test_regression_ob170095.py` (skip-if-workspace-absent pattern), and
put the DECaPS/VVV zero points and the APOGEE check in the changelog.

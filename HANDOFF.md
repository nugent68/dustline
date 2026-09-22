# dustline — HANDOFF (state at v0.8.0, 2026-09-22)

**2026-09-21 evening (0095 session, this Mac `/Users/nugent/claude/dustline-pkg`; venv python
3.14):** §6 DONE — `examples/ob170095` (README with the prototype comparison), `tests/
test_regression_ob170095.py`, `tools/prior_profile.py` (package → declens prior inputs),
CHANGELOG v0.7.1. Two catalog bugs fixed on the way (DECaPS `err_<b>` columns; VVV from the
brutus table — Data Lab has no VVV table), and the prototype's crowded-field parallax
treatment ported (`--plx-inflate`, posterior distances, in-field calibration; see the
changelog). Result: R_V 3.16 ± 0.23 (326 stars), clump A_V 2.71 at 8.4 kpc, prior
P_SEDdust_XP 4.43 kpc [3.35, 6.31] (prototype 4.41). Note the package clump column in DECam i
is 1.67 vs the prototype's five-estimator 1.86: the (J−Ks) anchor is the least leveraged
estimator (0.03 mag of clump colour = 0.1 in A_i); folding the i−Ks / g−i / LF estimators into
`clump.py` is the next improvement for bulge sightlines. The 0095 XP cache is now on this Mac
(`~/.cache/dustline/sightlines/ra+0267.86642_dec-033.13517_r5_*`).

For a fresh agent session in this repo (`/Users/nugent/claude/dustline-pkg`, GitHub
`nugent68/dustline`, public, tag `v0.8.0`). Self-contained; the science
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
corrections (per-model since v0.8.0), and a fine local (R_V, A_V) refinement around the
coarse optimum.
Photometric zero points are iterated (NIR frozen after pass 1; GALEX per T_eff bin).

## 2. Environment and caches

- venv: `.venv` (python 3.14; `source .venv/bin/activate`; `dustline` CLI on PATH). Tests:
  `python -m pytest -q -m "not network"` (49 pass at v0.8.0; the two regression tests use
  the cached bulge and COSMOS workspaces — a new corrections tag makes them skip until the
  fields are rerun and the pins updated).
- Assets (fetched on first use into `~/.cache/dustline/`, sha-verified, registry in
  `src/dustline/assets.py`): `newera_full_cache.npz` + `mist_v1.2_basic.npz` (v0.1.0),
  `newera_uvir_cache.npz` (v0.4.0, 214 MB, 900 Å–6 µm, the default), `template_corrections.npz`
  (v0.8.0, corrections tag `tca602f749`, 123 per-model entries + coarse bins + the R_V
  closure table).
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
| `calib.py` | template corrections: calibrator selection (DESI dwarfs D<250 pc; APOGEE giants D<1.2 kpc), `teff_from_bprp` (Mamajek locus + [Fe/H] term), `map_extinction` (Edenhofer+2023 via dustmaps), `deredden`, `ratio_spectrum`, `aggregate` (per-model entries where >= 6 calibrators share a model, else coarse bins by the **best-fit model**, unsmoothed), `load`/`tag`/`apply` (`MODEL_FALLBACK` = nearest covered model within 150 K / 0.5 dex / 0.5 dex), `rv_closure`, `LABEL_LOCK_SIGMA` |
| `ensemble.py` | `law_sample`, `measure_law` (closure-corrected; `rv_raw`), `default_law`, `dust_run`, `bridge_to_clump`, `foreground_column`, `band_ratio_at_rv` |
| `clump.py` | red-clump anchor in (J−Ks, Ks) |
| `models.py`, `extinction.py`, `filters.py`, `bands.py`, `plotting.py`, `cache.py` | NewEra cache + XP LSF + synthetic photometry; G23/F99 curves; 26 packaged filter curves (`data/filters`); band assembly; figures; `Workspace` |
| `catalogs/` | `footprint.plan` (survey selection by position), `gaia`, `ps1`, `decaps`, `decals`, `vvv`, `desi`, `apogee`, `galex`, `wise`, `datalab`, `vizier`, `xmatch` |

Tools: `tools/build_template_corrections.py` (stages `pull`/`fit`/`build` for `dwarfs` /
`giants`; calibration workspaces `ra+0000.00000_dec+000.00000_r0_8659c9cc` (dwarfs) and
`_83240b75` (giants)), `tools/inject_reddening.py` (closure test; `closure` stage writes
`rv_closure_k` into the corrections file), `tools/run_ob240669.py` (regenerates the bulge
example from the research repo's XP cache), `tools/fetch_svo_filters.py`.

## 4. v0.7.0–v0.8.0 in two paragraphs (details: docs/method.md)

**v0.7.0 — the injection test.** Calibrators reddened with a known law and refitted as a
field measured the T_eff–A_V–R_V coupling (+0.45 R_V, +0.10 A_V per +100 K, symmetric →
the ensemble median is unbiased, but **never split a field's R_V by fitted T_eff**) and
fixed: (1) `build_grid` cancelled the per-band template corrections at every A_V ≠ 0;
(2) the corrections table was smoothed at σ 20 nm and binned by the star's [Fe/H]/log g
instead of the model's; (3) high-S/N posteriors quantised to the grid (`_refine_star`).
A linear 1/A_V residual of the cool templates remained and is applied per star as the
closure table (`rv_closure_k`, `rv_raw` keeps the uncorrected median).

**v0.8.0 — the K-dwarf T_eff offset, closed as an A_V systematic.** Injection controls at
A_V = 0 and 1, T_eff free and locked, corrections on and off (table in `docs/method.md`):
locking T_eff removes the whole A_V bias but leaves R_V at 2.91 for 3.05 injected, so
**the cool-template R_V residual and the T_eff offset are independent and the coupling is
not a bias term in the R_V budget**. The offset survives at A_V = 0 (−25 K), is halved
rather than caused by the corrections (−67 K with them off), is not node binning or
coverage, and is the degeneracy floor of the 0.5 dex [M/H] grid — removing it needs a
finer [M/H] axis or a (T_eff, [M/H]) sub-grid refinement, not a better table. It is
carried as a systematic on **A_V only** (−0.05 at A_V 1, −0.02 at A_V 0, ≤ 0.02 when a
spectroscopic prior pins the star). Calibration changes: log g and [Fe/H] are locked as
well as T_eff when the calibrators are fitted, and the table carries 123 per-model
corrections preferred to the coarse bins. Injection closure on the new table: dwarfs
R_V 2.93 (ΔA_V −0.026), giants 2.98 (−0.014), against 2.7–2.9 before; the closure k
roughly halved.

Results on the v0.8.0 stack (closure-corrected, raw in brackets):

| field | mode | result | N | example |
|---|---|---|---|---|
| OB240669 (l 13°, b −2°, bulge) | law + clump | R_V 2.82 ± 0.50 [2.79]; clump A_V 3.34 at 6.1 kpc | 905 | `examples/ob240669` |
| ZTF19abqmpti (l 11°, b +24°) | law | R_V 3.50 ± 0.45 [3.44]; column 1.11 | 2,229 | `examples/ztf19abqmpti` |
| ZTF20abgaovd (l 16°, b +27°) | law | R_V 3.46 ± 0.79 [3.35]; column 0.60 | 1,559 | `examples/ztf20abgaovd` |
| COSMOS (b +42°) | column | F/G 0.023 ± 0.006; K/M 0.031; SFD 0.059 | 130 | `examples/cosmos` |

Field R_V is stable to ±0.05 across v0.6 → v0.8; the calibration work moved the error
budget (±0.05 closure, ±0.05 NIR zero points), not the answers.

## 5. Open issues (in priority order)

1. **The 0.5 dex [M/H] grid step** is the remaining floor on A_V for free-parameter cool
   dwarfs (v0.8.0, §4): −0.05 at A_V 1 for 4500–5100 K. The fix is a finer [M/H] axis (the
   UV–IR cache holds −2..+0.5 in 0.5 steps only) or a (T_eff, [M/H]) sub-grid refinement
   like `_refine` does in (R_V, A_V). It does **not** affect R_V — that was the v0.8.0
   result — so it is a column-mode / A_V-budget item, not an R_V one.
2. Closure calibrators are solar-neighbourhood ([Fe/H] −0.5..+0.3, [M/H]-dependent
   residual); bulge/thick-disk fields extrapolate.
3. Free giant fits sit 25–90 K below ASPCAP in reddened fields (no luminosity lever;
   was 50–120 K before the label-locked corrections).
4. Column-mode absolute zero point ±0.03 (XP ~0.03 mag too red in g−z for faint red stars);
   COSMOS now reads 0.036 *below* SFD and the injection control rules out a template
   zero point, so the remaining suspects are the maps and the calibrator population.
5. **PS1 band corrections are essentially uncalibrated**: the nearby bright calibrators hit
   the PS1 saturation cut, so `dm_PS1_i` rests on 73 dwarfs, `dm_PS1_z` on 71 and
   `dm_PS1_g` on 603, against 2,386 for 2MASS — `dm_PS1_*` is 0 at most nodes. Fields
   measure their own PS1 offsets (which is why this is not a blocker), but the injection
   closure inherits no PS1 information. A fainter calibrator sample (DESI reaches G 19)
   would fix it.
6. GALEX adds nothing at A_V < 0.1 and drops out where R_V is measurable; WISE is neutral.

## 6. Task for the next agent: OGLE-2017-BLG-0095 on the v0.8.0 stack

**Why**: 0095 is the sightline the method was built for (research prototype, 2026-09-17:
`~/claude/dustline/results/ob170095/`, `HANDOFF_dustline.md` §0). Its deliverable
(`P_SEDdust_XP`, the source-distance prior sent to Natasha) rests on a 2500–25000 Å cache,
solar-[M/H] templates, no corrections, no closure, the 0.25/0.1 grid — everything v0.7/v0.8
changed. A package run is both a regression of the package on the DECaPS/VVV/clump path
(only exercised by the bulge example so far) and an updated deliverable.

**Target**: RA 267.86642, Dec −33.13517 (l 357.0°, b −3.2°). `footprint.plan` gives
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
  corrections at A_V > 0, closure +0.03–0.1 for the cool giants).
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

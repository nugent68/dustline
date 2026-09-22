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
  `python -m pytest -q -m "not network"` (49 pass, 4 skip at v0.8.0; the three regression
  tests use the cached bulge / COSMOS / 0095 workspaces — a new corrections tag makes them
  skip until the fields are rerun and the pins updated, which is why `test_regression_
  ob170095.py` skips: 0095 is still pinned to `tc26438c9b`, §6).
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
| `ensemble.py` | `law_sample`, `measure_law` (closure-corrected; `rv_raw`), `default_law`, `dust_run`, `bridge_to_clump`, `foreground_column`, `band_ratio_at_rv`, `distance_posterior` + `parallax_inflation` (crowded fields, v0.7.1) |
| `clump.py` | red-clump anchor in (J−Ks, Ks) |
| `models.py`, `extinction.py`, `filters.py`, `bands.py`, `plotting.py`, `cache.py` | NewEra cache + XP LSF + synthetic photometry; G23/F99 curves; 26 packaged filter curves (`data/filters`); band assembly; figures; `Workspace` |
| `catalogs/` | `footprint.plan` (survey selection by position), `gaia`, `ps1`, `decaps`, `decals`, `vvv`, `desi`, `apogee`, `galex`, `wise`, `datalab`, `vizier`, `xmatch` |

Tools: `tools/build_template_corrections.py` (stages `pull`/`fit`/`build` for `dwarfs` /
`giants`; calibration workspaces `ra+0000.00000_dec+000.00000_r0_8659c9cc` (dwarfs) and
`_83240b75` (giants)), `tools/inject_reddening.py` (closure test; `closure` stage writes
`rv_closure_k` into the corrections file), `tools/run_ob240669.py` (regenerates the bulge
example from the research repo's XP cache), `tools/prior_profile.py` (package outputs →
the declens source-distance-prior inputs), `tools/fetch_svo_filters.py`.

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
6. **The clump column rests on one colour**: `clump.py` anchors on E(J−Ks) alone, and
   A_i/E(J−Ks) ≈ 3.8, so 0.03 mag of clump colour is 0.1 in A_i. On 0095 it gives A_i 1.67
   against the prototype's five-estimator (g−i, i−Ks, J−Ks, H−Ks, LF peaks) 1.86 ± 0.19.
   Folding the other estimators in is the next improvement for bulge sightlines (v0.7.1).
7. GALEX adds nothing at A_V < 0.1 and drops out where R_V is measurable; WISE is neutral.

## 6. Task for the next agent: 0095 on the v0.8.0 corrections

`examples/ob170095` (OGLE-2017-BLG-0095, RA 267.86642, Dec −33.13517, l 357.0°, b −3.2°, 5′)
was delivered in v0.7.1 — README, `tests/test_regression_ob170095.py`, `tools/prior_profile.py`,
and the declens source-distance prior `P_SEDdust_XP` 4.43 kpc [3.35, 6.31] — but on the
**v0.7.0 corrections (`tc26438c9b`)**, before the label-locked per-model table, so its
regression test skips on this stack.

**What to do**: rerun it on `tca602f749`. **The user is doing this from another machine
(2026-09-22), which holds the 0095 XP cache** — do not start it here without asking. This
Mac has no 0095 data despite the v0.7.1 note above (a run started on 2026-09-22 was stopped
at 1,800/2,995 spectra and its workspace deleted at the user's request) and nothing to seed
from, so here it would be a full fetch: 2,995 XP stars ≈ 90 min, then ~1 h of fits over 3
passes. On the machine with the cache the new corrections tag seeds its workspace from the
sibling, so it is a refit. **Never delete the workspace** on the machine that does the run.

```bash
source .venv/bin/activate
nohup dustline run 267.86642 -33.13517 --radius 5 --plx-inflate 1.7 --filter I \
      -o examples/ob170095/extinction_I.csv --plots examples/ob170095 \
      > /tmp/ob170095_v08.log 2>&1 &
python tools/prior_profile.py 267.86642 -33.13517 --radius 5 --plx-inflate 1.7 \
      -o examples/ob170095/prior
```

Then copy the workspace `result.json` to `examples/ob170095/law.json`, update the README
numbers (and its "Numbers are from the v0.7.0 stack" note), re-pin
`tests/test_regression_ob170095.py`, and rerun the in-field parallax calibration
(`law.json["plx_inflation_clump"]`) in case ×1.7 is no longer the right inflation.

**Expect**: R_V within ~0.05 of 3.16 (the v0.6→v0.8 field-to-field stability, §4) with a
smaller closure term (+0.10 → ~+0.04 for cool giants); the clump anchor and the zero points
should barely move. The prototype comparison to hold against is in the README: R_V 3.15 ±
0.23, clump A_i 1.86 ± 0.19 at 8.3 kpc, prior 4.41 kpc [3.38, 5.90]. Port nothing back to
the declens repo (`~/claude/dustline`) without the user's say-so.

**Also worth doing**: §5.6, the multi-estimator clump column — 0095 is the sightline where
the single-colour anchor is demonstrably 0.19 mag low, so it is the natural test case.

# Changelog

## v0.8.0 — 2026-09-22

### The K-dwarf T_eff offset: located, and shown not to touch R_V
v0.7.0's open item (4500–5100 K dwarfs fitting ~67 K cool with A_V −0.10) was chased with
injection controls at A_V = 0 as well as 1, with T_eff free and locked, and with the
corrections off (table in docs/method.md). Results: locking T_eff removes the whole A_V
bias but leaves R_V at 2.91 for 3.05 injected, so **the cool-template R_V residual is
independent of the T_eff offset and the ±0.45 R_V/100 K coupling is not a bias term** in
the R_V budget; the offset exists at A_V = 0 (−25 K), is halved rather than caused by the
corrections (−67 K with them off), and is the degeneracy floor of the 0.5 dex [M/H] grid
(stars landing on a more metal-rich model than their true [Fe/H] come back unbiased).
It is now carried as a systematic on **A_V only**: −0.05 at A_V = 1, −0.02 at A_V = 0 for
free-parameter K dwarfs, ≤ 0.02 for warm dwarfs and for any star pinned by a
spectroscopic prior (so column mode is unaffected — its template A_V zero point is
−0.007, and the COSMOS–SFD difference is not a template artefact).

### Calibration
- `build_template_corrections.py fit` locks **log g and [Fe/H]** to the grid values
  nearest the DESI / ASPCAP labels as well as T_eff (`calib.LABEL_LOCK_SIGMA`). With only
  T_eff locked, solar K dwarfs were calibrated on log g 5.5 / [M/H] +0.5 models that the
  science fit's soft priors then abandoned; the K-dwarf nodes now sit on log g 5.0 alone.
- `calib.aggregate` writes **per-model corrections** (exact T_eff, log g, [M/H] of the
  best-fit model) wherever ≥ 6 calibrators share a model — 123 entries over 3,474
  calibrators — and `apply` prefers them to the coarse (T_eff, [M/H] ≶ −0.4,
  dwarf/giant) bins, which mix log g 4.5–5.5 and [M/H] 0/+0.5 at a single node. Models
  without their own entry take the nearest covered model within `MODEL_FALLBACK`
  (150 K, 0.5 dex, 0.5 dex) instead of the coarse bin. `calib.tag` now hashes
  `model_key` and `MODEL_FALLBACK` too.
- New corrections asset, tag `tca602f749`. Injection closure on it: dwarfs R_V 2.93
  (ΔA_V −0.026, ΔT_eff −14 K over 2,439 stars), giants 2.98 (−0.014, −11 K over 1,482) —
  against 2.7–2.9 for the cool nodes on v0.7.0. The rebuilt closure table is
  correspondingly smaller: k = +0.017 (dwarfs) / +0.012 (giants) for the cool nodes
  against +0.029 before.
- Asset `template_corrections.npz` re-released (v0.8.0 tag, 494 kB, sha in
  `assets.py`); the bulge and COSMOS regression tests are re-pinned to the new stack.
  `examples/ob170095` (added in v0.7.1) is still measured on `tc26438c9b`, so
  `tests/test_regression_ob170095.py` skips until that field is rerun (HANDOFF §6).

### Examples (all four rerun with `--refit`; closure-corrected R_V, raw in brackets)
- `examples/ob240669` (bulge, A_V 2–5): R_V **2.82 ± 0.50** from 905 stars [2.79];
  clump column 3.34 ± 0.62 at 6.1 kpc; A_I 0.37 / 0.81 / 0.95 / 1.12 / 1.26 at
  1–5 kpc, bridged to 1.94.
- `examples/ztf19abqmpti` (l 11°, b +24°, A_V ≈ 1): R_V **3.50 ± 0.45** from 2,229
  stars [3.44]; column 1.11 ± 0.26 beyond 1 kpc against the rescaled SFD 1.07.
- `examples/ztf20abgaovd` (l 16°, b +27°, A_V ≈ 0.6): R_V **3.46 ± 0.79** from 1,559
  stars [3.35]; column 0.60 ± 0.19, map 0.55. Its APOGEE giants moved from −119 K
  (v0.6) / −50 K (v0.7) to −25 K against ASPCAP.
- `examples/cosmos` (column mode): F/G **0.023 ± 0.006** (MAD 0.046, N 130), K/M 0.031
  — 0.036 below SFD, now demonstrably not a template A_V zero point (above).

The field R_V is stable to ±0.05 across the whole v0.6 → v0.8 calibration overhaul
(bulge 2.76 → 2.80 → 2.82; ZTF19 3.44 → 3.48 → 3.50; ZTF20 3.55 → 3.50 → 3.46), which
is the point of the exercise: the closure work moved the *error budget*, not the answer.

## v0.7.1 — 2026-09-21

### Crowded-field parallaxes (`--plx-inflate`, posterior distances)
Gaia DR3 parallax errors are underestimated in crowded fields — by up to ×4 towards the bulge
above ~300 sources/arcmin² (Luna+2023, arXiv:2307.13719), by up to 80 % for a star with a
neighbour within 4" (El-Badry+2021) — and 1/ϖ is biased for the distant stars, which are
exactly the ones that set the outer run. Ported from the research prototype (2026-09-18):
- `Sightline(plx_inflate=…)` / `dustline run --plx-inflate X`: multiplies `parallax_error` in
  the radius prior and the D errors (config key `plx_inflate`, so a new workspace, seeded).
- `ensemble.distance_posterior`: per star, the inflated-parallax likelihood × the photometric
  distance the fit already implies (log D = log R_MIST(best model) − ½ log C, σ 0.08–0.15 dex;
  `fit._summarize` now writes `logR_iso`, `sig_logR_iso`) × a uniform prior → `D_med`,
  `D_lo`, `D_hi`, `D_frac`. `dust_run` uses `D_med` with inflated S/N > 3 and posterior
  half-width < 25 % (unchanged 1/ϖ path when a fit table lacks the columns); `law_sample`
  and `plot_run` use the inflated S/N.
- `ensemble.parallax_inflation`: the in-field calibration from the red-clump window stars
  (all at D_RC): robust/plain std of (ϖ − zp − 1/D_RC)/σ_ϖ, reported as
  `result.json["plx_inflation_clump"]` on bulge sightlines with the value to rerun with.

### Catalogs
- DECaPS DR2: the Data Lab error columns are `err_<b>`, not `err_mag_<b>` (the query had
  failed silently and the bulge path fell back to 2MASS-only).
- VVV: Data Lab has no VVV table (only VHS, which excludes the VVV footprint); JHKs now
  come from `decaps_dr2.stellar_inference` (mag_6/7/8 = VVV J/H/Ks of each DECaPS object,
  matched by its Gaia id, else by position), with the `vvv_dr4.vvvsource` query as fallback.

### Example: `examples/ob170095` (OGLE-2017-BLG-0095, l 357°, b −3°, 5′, DECaPS + VVV)
R_V 3.16 ± 0.23 from 326 stars [raw 3.07]; zero points DECaPS g −0.078, r −0.092, i −0.033,
z −0.041, Y −0.051, VVV J −0.013, H +0.006, Ks +0.037; in-field parallax-error underestimate
×1.65 (robust) / ×1.88 (std) from 1,657 clump stars; clump E(J−Ks) 0.44 → A_V 2.71 ± 0.44
at 8.4 kpc; A_I 0.67 / 1.23 / 1.32 / 1.40 / 1.55 at 1–5 kpc. `tools/prior_profile.py`
writes the declens source-distance-prior inputs (DECam i profile + source band ratios):
P_SEDdust_XP 4.43 kpc [3.35, 6.31] (prototype 4.41 [3.38, 5.90]). No APOGEE star with
S/N > 30 inside 5′. Regression: `tests/test_regression_ob170095.py`.

## v0.7.0 — 2026-09-21

### Reddening-injection closure test (`tools/inject_reddening.py`)
The template calibrators (2,439 DESI dwarfs, 1,482 APOGEE giants) are reddened with a
known G23 law (A_V 1, R_V 3.05; also 0.5 / 2 and R_V 3.8) and refitted as a field is
(T_eff free, priors, corrections, parallax S/N degraded to 10). It measured the
T_eff–A_V–R_V coupling (+0.45 in R_V, +0.10 in A_V per +100 K; symmetric, so the
ensemble median is unbiased; the R_V–A_V correlation is the per-star error ellipse) and
found three defects, all fixed (docs/method.md):
- **Grid**: the per-band template corrections were cancelled at every A_V ≠ 0
  (`A_band` was measured against the corrected band fluxes). Grids rebuilt (`b2` key).
- **Corrections table**: σ = 20 nm smoothing erased the TiO-band-scale residuals of the
  cool templates (−24 % at 400 nm left in the 3500 K node → XP-only ΔA_V +0.26), and the
  bins were by the star's [Fe/H] / log g while `apply` uses the model's. Now unsmoothed,
  binned by the best-fit model; per-node closure exact. New asset (tag `tc26438c9b`).
- **Estimator**: high-S/N stars' posteriors sit on one grid node and the 0.1 A_V grid
  samples the A_V–R_V valley as a zig-zag; `fit._refine_star` adds a fine local
  (R_V, A_V) grid (0.025 × 0.01) around the best node of each model carrying the
  posterior. Per-star R_V / A_V (and their errors, A_<band>) are now continuous.

What remains is a linear, 1/A_V template residual of the cool nodes (dwarfs < 5100 K,
giants < 4500 K read R_V 2.7–2.9 for 3.05; warm nodes 3.0–3.05), carried as a closure
table in the corrections asset (`rv_closure_k`) that `measure_law` applies per star
(`rv_raw` keeps the uncorrected median; `n_closure` counts the corrected stars).

### Other
- `dustline run … --refit` / `Sightline.run(refit=True)`: redo the fits and ensemble
  products from the cached catalogs and XP spectra (no refetch).
- `calib.aggregate(stat="mean")` option (clipped mean; not used — it did not change the
  closure).
- Asset `template_corrections.npz` re-released (v0.7.0 tag; the corrections tag
  `tc26438c9b` keys the grids and workspaces, the closure table rides along).

### Examples (all rerun on the v0.7.0 stack; closure-corrected R_V, raw in brackets)
- `examples/ob240669` (bulge, A_V 2–5): R_V 2.80 ± 0.54 from 883 stars [2.73]; clump
  column 3.36 ± 0.63 at 6.1 kpc; A_I run within 0.03 of v0.6.
- `examples/cosmos` (column mode): F/G 0.029 ± 0.004 (MAD 0.052, N 132), K/M 0.039.
- `examples/ztf19abqmpti` (new; l 11°, b +24°, A_V ≈ 1): R_V 3.48 ± 0.48 from 2,228 stars
  with A_V ≥ 0.5 [3.37]; column 1.07 ± 0.25 beyond 1 kpc = the rescaled SFD. The per-star
  R_V–A_V–T_eff trend seen here (3.1 → 4.1 across A_V bins at fixed distance) is what
  the injection test was built to explain.
- `examples/ztf20abgaovd` (l 16°, b +27°, A_V ≈ 0.6): R_V 3.50 ± 0.84 from 1,394 stars
  [3.32] (was 754 stars: the band-correction fix raised the field's A_V by ~0.1); column
  0.57 ± 0.18 beyond 1 kpc, on the 3D map (0.55).

## v0.6.0 — 2026-09-20

### Added
- **Giant template corrections**: the corrections table gained a log g class axis
  (dwarfs ≥ 3.5 / giants < 3.5). 1,482 nearby APOGEE DR17 giants with Gaia XP (D < 1.2 kpc,
  |b| > 30°, > 20° below 4500 K; 120 per 100 K node from 4200 to 5200 K), dereddened with
  the Edenhofer+2023 map (median A_V 0.08) and fitted at A_V = 0 with T_eff and log g locked
  to ASPCAP (IRFM scale). Closure: dereddened giants 0.00 ± 0.01 at 4300–5300 K and their
  map extinction recovered to 0.01–0.02; weaker at the sparse ends (< 4300 K: +0.06,
  > 5300 K: −0.06). `tools/build_template_corrections.py {pull,fit} giants`, joint `build`;
  `calib.select_giant_calibrators`, `calib.LOGG_CLASSES`. Legacy 3-D tables still load.
- `gaia.by_ids_datalab` / `gaia.cone_datalab` (+ `_attach_twomass`): NOIRLab Data Lab
  fallbacks (gaia_dr3.gaia_source + twomass.psc) used automatically when the Gaia archive
  fails (it returned 408s for hours on 2026-09-20); `gaia.by_ids` retries.
- `catalogs/apogee.py`: APOGEE DR17 ASPCAP parameters (Data Lab `sdss_dr17.apogee2_allstar`,
  joined on Gaia EDR3 source_id); log g / [Fe/H] priors where DESI has none and an
  IRFM-scale T_eff check recorded in the law dict (`apogee_check`).
- `Sightline(min_av=...)` / `--min-av`: the R_V law threshold (0.5 is workable with priors);
  `Sightline(desi_teff=True)` / `--desi-teff` re-enables the DESI T_eff label as a prior.
- `examples/ztf20abgaovd`: the first mid-latitude SN Ia sightline (l 16°, b +27°, A_V ≈ 0.5):
  R_V = 3.56 ± 0.89 (MAD) from 775 stars with A_V ≥ 0.5; column 0.44 vs SFD 0.47–0.55.

### Changed
- In law mode the DESI T_eff label is no longer used as a prior (its offset from the colour
  scale swings by ±200 K with S/N); DESI/APOGEE log g and [Fe/H] priors remain, and
  `spec_prior_chi2` accepts any subset of the three terms.
- Bulge example on the giant-corrected templates: R_V 2.76 ± 0.50 (857 stars; was 2.80 ±
  0.46), clump 3.41, A_I(D) 5–8 % lower at 3–5 kpc where the stars are giants.

### Known limitation
- A free giant fit in a *reddened* field still lands ~120 K cooler than ASPCAP
  (ZTF20abgaovd: dwarfs −22 K, giants −119 K): giants lack the parallax–luminosity lever
  against the T_eff–A_V degeneracy, so the far giants of a mid-latitude field need an
  external T_eff. The bulge (A_V 3–5) is much less sensitive to it.

## v0.5.0 — 2026-09-20

### Added
- **Empirical dwarf-template corrections** (`dustline.calib`, asset
  `template_corrections.npz`, `tools/build_template_corrections.py`): per T_eff grid node
  and [Fe/H] bin, the XP obs/model ratio spectrum and per-band offsets of 2,451 nearby
  DESI×XP dwarfs, dereddened star by star with the Edenhofer+2023 3D map (`dustmaps`,
  optional dependency, only for building) and fitted at A_V = 0 with T_eff locked;
  multiplied into the dwarf models (log g ≥ 3.5) of the fit grid.
- **T_eff lock on the empirical colour scale** (column mode): each star's T_eff prior is
  the Mamajek-locus T_eff of its dereddened Gaia BP−RP with a [Fe/H] term
  (`calib.teff_from_bprp`, packaged `data/mamajek_dwarf_locus.dat`), locked to the
  nearest grid node; DESI supplies log g/[Fe/H]. The DESI T_eff label proved unstable at
  the 100 K level (S/N- and brightness-dependent) and is no longer used as the T_eff prior.
  `spec_prior_chi2` now takes T_eff-only priors (log g / [M/H] terms optional).
- A_V grid in 0.01 steps below 0.5 and extended to −0.1 (unbiased near zero).
- `fit.fit_stars(av_fixed=...)`, `gaia.by_ids`, `ps1.by_positions`,
  `datalab/vizier.positions_query`, `calib.map_extinction/deredden`.

### Fixed
- PS1 mean-PSF magnitudes brighter than the saturation limits (g 14.5, r 15, i 15, z 14,
  y 13) are dropped (`ps1.SATURATION`), also when reading cached tables.
- The fit-grid cache key now includes the A_V axis.

### Result
- COSMOS F/G column 0.027 ± 0.005 (MAD 0.050; was 0.072 ± 0.009, MAD 0.065); K/M stars
  0.026 (was 0.185): the K/M-dwarf systematic is resolved — it was the DESI T_eff label
  plus NewEra's too-blue K-dwarf SEDs, not dust or XP. Zero point ±0.03 systematic
  (SFD 0.059). Bulge unchanged (R_V 2.80 ± 0.46, clump 3.27).

## v0.4.0 — 2026-09-19

### Added
- **UV and mid-IR photometry**: `catalogs/galex.py` (GUVcat_AIS FUV/NUV via VizieR TAP,
  `catalogs/vizier.py`; colour-gated to hot stars, artefact flags, 0.05 mag floor) and
  `catalogs/wise.py` (AllWISE W1/W2 via Data Lab; cc/ext flags, saturation cut).
  `SurveyPlan.uv` / `.mir`, `Sightline(uv=, mir=)`, CLI `--no-uv` / `--no-mir`. GALEX
  zero points are iterated per T_eff bin (`fit.uv_offsets_by_teff`; the NewEra NUV flux
  is ~0.35 mag too bright for G stars and the bias is T_eff dependent) and carry a
  0.10 mag model systematic; WISE joins the frozen-after-pass-1 NIR set. `foreground_column` splits the
  clean sample by NUV availability. GALEX FUV/NUV filter curves added.
- **UV-IR model cache** `newera_uvir_cache.npz` (`tools/build_newera_uvir_cache.py`):
  4,366 NewEra models from the HSR files' LSR spectra, 900-25000 Å at 2 Å + 25000-60000 Å
  at 20 Å, [M/H] −2..+0.5 — GALEX/WISE coverage and the metal-poor extension.
  `models.DEFAULT_CACHE` / `$DUSTLINE_MODEL_CACHE` select the cache; `models.cache_covers`
  guards every band (UV/MIR are skipped, with a message, when the cache does not reach
  them); the grid cache key and the workspace config carry the cache identity. With
  spectroscopic priors the [M/H] axis now spans whatever the cache holds; stars *without*
  a prior are confined to [M/H] −0.5..+0.5 (`fit.MH_FREE_RANGE`) — free metallicity
  absorbs residuals. **The UV-IR cache is the default model set** (`models.DEFAULT_CACHE`);
  on the OB240669 bulge field it reproduces the old cache to ΔR_V = 0.04 and ΔA_V < 0.05
  at every distance (per-star ΔA_V median −0.009, ΔT_eff −3 K).
- GALEX zero points are iterated per T_eff bin from the previous pass's T_eff
  (`fit.uv_offsets_by_teff`, recorded in `phot_offsets.json`); `UV_SYS = 0.15` mag from the
  measured NUV residual scatter. On COSMOS the UV adds nothing at A_V ≈ 0.06 (per-star
  ΔA_V 0.000 ± 0.006) — see the README for where it should help.
- `Workspace.seed_from_sibling`: a new option set at the same position reuses the Gaia
  cone, XP spectra and per-survey catalog matches of an existing workspace instead of
  refetching.
- `examples/cosmos` regenerated on the v0.4.0 stack (clean F/G column 0.072 ± 0.009,
  MAD 0.065; SFD 0.059); `examples/ob240669` regenerated on the new cache and the fixed
  PS1 fetch.

### Fixed
- **PS1 fetch**: the MAST catalogs API does not order rows between pages, so the paged
  fetch dropped or duplicated objects at every page boundary (three identical queries:
  page overlaps of 0, 180 and 6,936 rows; one refetch of the OB240669 field lost 16 % of
  the matches). `catalogs/ps1.py` now tiles a full cone into seven sub-cones recursively
  and deduplicates on `objID` — no paging.
- The packaged filter curves (`src/dustline/data/filters/*.dat`) were excluded from the
  repository (and from wheels built from it) by the `data/` ignore pattern since v0.1.0;
  the pattern is now root-only and the 26 curves are tracked.

## v0.3.0 — 2026-09-19

### Added
- **Column mode for high-latitude fields** (`plan.mode == "column"`, |b| > 30°): no
  measurable law, so G23 at R_V = 3.1 is assumed (`ensemble.default_law`) and the product
  is the total foreground A_V from the parallax stars behind the dust
  (`ensemble.foreground_column`: median, bootstrap error, MAD, splits by spectroscopic
  prior and T_eff, and a cell map of the field). `ExtinctionResult.mode` / `.column`,
  `column_av.png` (`plotting.plot_column`), CLI `--ref-av` for a reference line.
  Also fixes the latent crash when `measure_law` found < 10 reddened stars: the run now
  falls back to the assumed law instead of raising. In column mode the per-star fits hold
  R_V at the assumed value (`fit_stars(rv_fixed=...)`) — unconstrained at A_V ~ 0.1.
  The headline column uses the T_eff ≥ 5500 K stars: on COSMOS the K/M dwarfs return a
  spurious A_V ≈ 0.2 (cool-template / faint-XP systematic).
- `examples/cosmos/`: the COSMOS column-mode run (0.5°, 592 XP stars, 422 DESI priors):
  F/G column A_V = 0.086 ± 0.006 vs SFD 0.05; the `--no-spectro` A/B gives 0.170 ± 0.011
  with twice the scatter. `tests/test_regression_cosmos.py` pins it.
- **DESI DR1 MWS spectroscopic template priors** (`catalogs/desi.py`, Data Lab
  `desi_dr1.mws` joined on Gaia `source_id`; `plan.spectro == "desi"` for Dec > −25,
  |b| > 15): per-star Gaussian priors on (T_eff, log g, [Fe/H]) added to the fit
  (`fit.spec_prior_chi2`), and the grid's [M/H] axis opened to −0.5/0/+0.5 with the
  MIST radius prior evaluated per metallicity. New per-star columns `mh`, `mh_err`,
  `mh_best`, `spec_prior`, `mh_clamped`, `teff_spec` … `spec_snr`. `Sightline(...,
  spectro_priors=False)` / `--no-spectro` disables them (separate cache key).
- `Sightline(..., freeze_offsets=True)` / `--freeze-offsets`: single fit pass with the
  photometric zero points held at 0 (offsets still measured and recorded); the law dict
  now carries `phot_offsets`, `offsets_frozen`, `n_spec_prior` and the plan's
  `spectro`/`mode`.

### Changed
- `SurveyPlan` gained `spectro` and `mode` fields (defaults "none"/"law"; the workspace
  cache key only changes when they are non-default, so existing law-mode caches are kept).
- `fit.fit_stars` sizes its batch from the number of grid models (memory-bounded).

### Examples and docs
- `examples/ob240669/`: the packaged OGLE-2024-BLG-0669 run (A_I(D) table, law JSON,
  the two figures, README with the numbers); `tools/run_ob240669.py` now regenerates it.
  The README quickstart and API docstring use this sightline.
- `tests/test_regression_ob240669.py` pins the law, clump anchor and bridged run to the
  example (runs when the sightline workspace is in the local cache); it replaces
  `tests/test_regression_ob170095.py`, and the OGLE-2017-BLG-0095 references are gone
  from the README and method notes.
- `docs/method.md`: the "no extrapolation beyond the last bin" statement now describes
  the red-clump bridge on bulge sightlines.

## v0.2.0 — 2026-09-19

### Added
- **Red-clump bulge bridge** (`dustline.clump`, `ensemble.bridge_to_clump`): on bulge
  sightlines the red clump is located in the (J−Ks, Ks) CMD of the full field (2MASS or
  VVV), converted to an A_V column with the measured law ratios and to a distance D_RC
  from the dereddened Ks (NewEra clump model: 4700 K, log g 2.5, [M/H] 0, log L 1.72,
  M_Ks ≈ −1.61). The A_V(D) run is extended past the last parallax bin by a linear
  bridge to the clump column at D_RC and held flat to 12 kpc; bridged rows are flagged
  (`bridged = True`, `n_stars = 0`). When D_RC already lies inside the parallax range
  the run is still extended (0.5 kpc ramp, then flat at the column). Credibility gates:
  ≥ 500 stars with J and Ks, ≥ 200 in the search box, ≥ 60 in the clump window,
  E(J−Ks) > 0.15, 5 < D_RC < 12 kpc — otherwise no bridge is applied.
- `ExtinctionResult.clump_anchor` (the anchor dict, or None) and
  `ExtinctionResult.stars` (the per-star fit table).
- **Diagnostic figures** (`dustline.plotting`, `ExtinctionResult.plots()`,
  `dustline run --plots [dir]`): `law_rv.png` (per-star R_V vs A_V coloured by distance,
  plus the R_V histogram) and `extinction_run_<band>.png` (per-star A_X vs D coloured by
  fitted T_eff, running median with 16–84 % band, dashed clump bridge and anchor).
  matplotlib is an optional dependency (`pip install "dustline[plot]"`).
- `tools/run_ob240669.py`: seeds a workspace for OGLE-2024-BLG-0669 from the research
  repo's cached Gaia/XP data and runs the full v0.2 pipeline.

### Changed
- **Widened bulge window** for the clump anchor from |l| < 10° to |l| < 20° (|b| < 10°
  unchanged) so the long bar is covered; the clump credibility gates make the wider
  window safe.
- `Sightline.run()` skips the Gaia XP fetch when calibrated spectra already exist in the
  workspace (seeded or completed workspaces), unless `force=True`.

## v0.1.0 — 2026-09-19

- First release: pip-installable extinction-law (R_V) + 3D dust-run A_X(D) package.
  Model assets (`newera_full_cache.npz`, `mist_v1.2_basic.npz`) attached to the GitHub
  release and fetched on first use with sha256 verification.

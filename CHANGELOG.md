# Changelog

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

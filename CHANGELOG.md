# Changelog

## Unreleased

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

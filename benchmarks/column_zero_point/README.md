# Column-mode zero point: COSMOS vs the DESI stellar-reddening map (2026-10-01)

## COSMOS (30', column mode, R_V 3.1 assumed)

| estimate | A_V |
|---|---|
| dustline column, F/G beyond 0.5 kpc (N 131) | 0.024 +/- 0.006 (all stars 0.027) |
| dustline, photometric zero points frozen | 0.019 +/- 0.005 (all stars 0.024) |
| DESI E(g-r) (Zhou+24, arXiv:2409.05140), G23 3.1, DECam | 0.046-0.048 +/- 0.003 |
| DESI E(r-z), same | 0.055-0.059 +/- 0.004 |
| SFD x 0.86 (Schlafly & Finkbeiner 2011) | 0.050-0.051 |
| Edenhofer+23 to 1.25 kpc | 0.081 |

dustline reads 0.02-0.03 below DESI/SFD (the package handoff already carried COSMOS 0.036 below SFD
and a +/- 0.03 column zero point). Freezing the zero points lowers it further, so the zero-point /
screen degeneracy is not the cause.

## The calibrators' dust (`tools/calibrator_dust.py`, `calibrator_dust.{txt,csv}`, run at NERSC)

The dwarf template corrections are measured on 2,901 DESI dwarfs (177 pc [111, 229], |b| 56 deg)
dereddened with Edenhofer+23 (A_V in front, median 0.026) and fitted at A_V = 0. At their
positions Edenhofer's total column agrees with DESI (0.077 vs 0.087 g-r / 0.073 r-z; SFD x 0.86
0.061). With Edenhofer's own distance profile (37 % of the column in front) no dust is missing
in front of them (+0.000 vs DESI g-r, -0.005 r-z, -0.007 SFD); only a smooth 125 pc exponential
layer (69 % in front; a poor model inside the Local Bubble) would leave 0.007-0.019 uncounted.
The template corrections' zero point is therefore good to ~0.01 and does not explain most of the
COSMOS deficit.

## Metallicity

COSMOS F/G stars with DESI [Fe/H] (N 90, median -0.55): no monotonic A_V trend with [Fe/H]
(0.026-0.057 across bins, +/- 0.01-0.02 each). Not the cause at this precision.

## Open

The ~0.025 deficit is unexplained. Next: column-mode runs at ~20 high-latitude fields spanning
DESI E(B-V) 0.01-0.1 (cheap at NERSC: few XP stars each) - slope vs intercept of dustline vs DESI
separates a zero-point offset from a scale error.

## 20 high-latitude fields (2026-10-01; `tools/column_fields.py`, `tools/field_batch.py`, NERSC)

19 fields (|b| > 30, Dec > -25, DESI DR1 MWS coverage, one per log bin of DESI E(B-V) 0.011-0.089)
plus COSMOS, 30' each, run at NERSC in column mode. Straight-line fits dustline = a + b x reference:

| T_eff lock | vs DESI E(g-r) | vs DESI E(r-z) | vs SFD x 0.86 |
|---|---|---|---|
| BP-RP colour (v0.5+ default), `column_fields/` | a +0.003, b 0.21 +/- 0.06 | a -0.001, b 0.25 +/- 0.07 | a -0.002, b 0.31 +/- 0.09 |
| DESI label (`--desi-teff`), `column_fields_desi/` | a +0.045 +/- 0.009, b 0.87 +/- 0.09 | a +0.028 +/- 0.007, b 1.01 +/- 0.07 | a +0.027 +/- 0.007, b 1.15 +/- 0.07 |

- The colour lock recovers only 20-30 % of the column: a reddened star is read as a cooler one, and
  the iteration (T_eff from BP-RP dereddened by the previous pass's A_V, snapped to the 100 K grid
  node by the 1 K prior) stalls - on F24038 it sits at A_V 0.063-0.071 over 11 passes while the
  DESI-label lock gives 0.236 (DESI 0.275, SFD x 0.86 0.171). The COSMOS "deficit" was this.
- The DESI-label lock has the right scale (slope ~1) and a +0.03 offset, consistent with the DESI
  T_eff scale running ~30-50 K hot for these F/G stars (~0.09 A_V per 100 K); the package's
  calibration already found DESI 50-85 K hot for K dwarfs.
- Fix in the package (e3aa6e8): `desi_teff` now replaces the colour lock in column mode (it was
  overwritten by it). Not yet the default; the +0.03 offset needs a DESI T_eff-scale correction.

## DESI T_eff on the colour scale (2026-10-01; `tools/desi_teff_scale.py`, `calib.desi_teff_offset`, cd899f0)

The template corrections place every calibrator at the colour T_eff of its dereddened BP-RP, so a
fit locked to the raw DESI label borrows the template of a hotter star. dT = T_DESI - T_colour,
measured on ~70k DESI DR1 dwarfs/subgiants at |b| > 40 within 1 kpc (Edenhofer-dereddened):
+65..+150 K for 4250-5750 K and +5..+30 K for 5750-6250 K at S/N > 60, larger at S/N 7-35,
erratic below S/N 7 (`desi_teff_scale/`). With `desi_teff` the label is now shifted by it (below
S/N 7 only a 250 K prior). Same 20 fields (`column_fields_desi_scaled/`):

| | vs DESI E(g-r) | vs DESI E(r-z) | vs SFD x 0.86 |
|---|---|---|---|
| fit a, b | -0.001 +/- 0.010, 0.95 +/- 0.09 | -0.019 +/- 0.009, 1.10 +/- 0.07 | -0.020 +/- 0.011, 1.28 +/- 0.11 |
| median difference (ratio) | +0.004 (1.03) | +0.000 (1.00) | +0.003 (1.03) |

Field-to-field scatter about DESI 0.015-0.022 in A_V, comparable to DESI's own g-r vs r-z
disagreement. COSMOS: 0.039 +/- 0.008 (DESI 0.045-0.055, SFD x 0.86 0.051). The scale-corrected
DESI lock is the column mode to use where DESI labels exist; the colour lock (still the default)
recovers only 0.2-0.3 of the column.

## Alpha-corrected DESI [M/H] prior (2026-10-02; 220d0fe, `column_fields_alpha/`)

`catalogs.desi.to_mh` now feeds [M/H] = [Fe/H] + log10(0.638 x 10^[a/Fe] + 0.362) as the prior
(benchmarks/desi_labels: raw DESI [Fe/H] reads 0.12-0.16 dex low against the NewEra fits). Same 20
fields, scale-corrected DESI T_eff lock: per-field column change median -0.0007 (range -0.016 ..
+0.002); still unbiased - median difference +0.0016 (DESI g-r), -0.0026 (r-z), -0.0013 (SFD x 0.86),
slope 0.92 +/- 0.09 vs DESI g-r. With T_eff locked the [M/H] prior has little leverage on A_V.

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

# Example: ZTF19abqmpti (mid-latitude SN Ia sightline, A_V ≈ 1, measured R_V)

Sightline of the ZTF DR2 SN Ia ZTF19abqmpti (z = 0.076): RA 251.391, Dec −6.325
(l 11.4°, b +24.3°), **30′ radius**. Law mode with R_V free on the v0.8.0 stack (UV–IR
cache, label-locked per-model template corrections `tca602f749`, sub-grid refinement,
closure-corrected law): 18,974 Gaia sources, 2,629 XP stars fitted;
DESI DR1 log g/[Fe/H] priors for 404 stars (no APOGEE, 15 usable GALEX NUV — too
reddened for the UV). Reference foreground: SFD E(B−V) = 0.403 → A_V 1.25 (1.07 with the
S&F11 rescaling); Edenhofer+2023 3D map 0.88 to 1.25 kpc.

```bash
dustline run 251.391 -6.325 --radius 30 --min-av 0.5 -o extinction_I.csv --plots .
```

## Result

- **R_V = 3.50 ± 0.45** (median ± star-to-star MAD, **2,229 stars** with A_V ≥ 0.5; per-star
  errors 0.16–0.51, ~±0.02 statistical on the median). Uncorrected median 3.44 (dwarfs
  3.41, N 2,019; giants 3.59, N 210); the closure correction adds +0.06 to the median — a
  quarter of the sample is on cool templates (< 5000 K, raw 3.28). No pile-up at the R_V
  grid edges.
- **Column**: A_V = 1.11 ± 0.26 (MAD) beyond 1 kpc, against the rescaled SFD 1.07;
  A_I(D) = 0.74 at 0.5 kpc, 0.68–0.71 at 1–2.5 kpc, rising to 0.85 by 6.5 kpc.
- v0.7.0 (coarse-binned corrections, T_eff-only label lock, larger closure): 3.48 ± 0.48
  from 2,228 stars [raw 3.37], column 1.07. v0.6.0 (no band corrections at A_V > 0,
  unrefined estimator, no closure): 3.44 ± 0.49 from 2,186 stars, column 1.03. The field
  R_V is stable to ±0.03 across all three calibrations.

## The R_V–A_V trend (a fit systematic, not dust)

*(v0.6.0 numbers; the trend is unchanged on the v0.8.0 stack — raw R_V 2.97 / 3.32 /
3.52 / 4.02 for A_V 0.5–0.8 / 0.8–1.1 / 1.1–1.5 / > 1.5. The reddening-injection test
reproduces it on the calibrators and identifies it as the fit's T_eff–A_V–R_V error
ellipse, symmetric about the truth (the injected A_V = 1 calibrators show the same
2.82 → 3.13 rise across their own fitted-A_V quartiles); the cool-star deficit is the
template residual the closure table corrects. See `docs/method.md`.)*

| by A_V bin (all D) | R_V | at fixed D 1–2 kpc | by T_eff (A_V 0.8–1.3) | R_V |
|---|---|---|---|---|
| 0.5–0.8 (N 348) | 3.11 | 3.08 | 3500–4500 K (N 124) | 2.79 |
| 0.8–1.1 (N 862) | 3.39 | 3.40 | 4500–5000 K (N 298) | 3.40 |
| 1.1–1.5 (N 838) | 3.55 | 3.62 | 5000–5500 K (N 454) | 3.46 |
| > 1.5 (N 138) | 3.84 | 4.10 | 5500–6500 K (N 517) | 3.54 |

At fixed A_V there is no trend with distance (3.38 / 3.48 / 3.52 / 3.47 at 0.5–1 / 1–2 /
2–4 / 4–10 kpc), while at fixed distance the trend with A_V is intact and the fitted
T_eff rises along it — so it is the T_eff–A_V–R_V coupling of the per-star fit, strongest
for the cool stars whose template corrections rest on the sparsest calibrator nodes. The
5000–6500 K stars, where the calibration is densest, give 3.46–3.54. With the coupling
measured on the calibrators the quote is **R_V = 3.5 ± 0.02 (stat) ± 0.1 (sys)**
(closure ±0.05, NIR zero points ±0.05; v0.8.0 showed the T_eff offset is not a bias term
in R_V). Compare the bulge (2.80 ± 0.54 at A_V 2–5) and
ZTF20abgaovd (see its example) 5° away on the same side of the inner Galaxy.

## Files

`extinction_I.csv`, `law.json`, `law_rv.png`, `extinction_run_I.png`.

![law](law_rv.png)

![run](extinction_run_I.png)

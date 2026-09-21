# Example: ZTF19abqmpti (mid-latitude SN Ia sightline, A_V ≈ 1, measured R_V)

Sightline of the ZTF DR2 SN Ia ZTF19abqmpti (z = 0.076): RA 251.391, Dec −6.325
(l 11.4°, b +24.3°), **30′ radius**. Law mode with R_V free on the v0.6.0 stack (UV–IR
cache, dwarf + giant template corrections): 18,974 Gaia sources, 2,629 XP stars fitted;
DESI DR1 log g/[Fe/H] priors for 404 stars (no APOGEE, 15 usable GALEX NUV — too
reddened for the UV). Reference foreground: SFD E(B−V) = 0.403 → A_V 1.25 (1.07 with the
S&F11 rescaling); Edenhofer+2023 3D map 0.88 to 1.25 kpc.

```bash
dustline run 251.391 -6.325 --radius 30 --min-av 0.5 -o extinction_I.csv --plots .
```

## Result

- **R_V = 3.44 ± 0.49** (median ± star-to-star MAD, **2,186 stars** with A_V ≥ 0.5; per-star
  errors 0.32–0.50, ~±0.02 statistical on the median). Dwarfs 3.44 (N 1,991), giants 3.48
  (N 195). No pile-up at the R_V grid edges below A_V 1.1 (< 1 %).
- **Column**: A_V = 1.03 ± 0.27 (MAD) beyond 1 kpc, matching the rescaled SFD (1.07);
  A_I(D) = 0.70 at 0.5 kpc, 0.64 at 1.5–2.5 kpc, rising to 0.82 by 6.5 kpc.

## The R_V–A_V trend (a fit systematic, not dust)

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
5000–6500 K stars, where the calibration is densest, give 3.46–3.54. Quote
**R_V = 3.4 ± 0.02 (stat) ± 0.3 (sys)** until the coupling is calibrated out (an external
T_eff for the reddened stars, or the R_V–T_eff correlation measured on calibrators).
Compare the bulge (2.76 ± 0.50 at A_V 2–5) and ZTF20abgaovd (3.55 ± 0.85 at A_V ≈ 0.6, 27°
away on the same side of the inner Galaxy).

## Files

`extinction_I.csv`, `law.json`, `law_rv.png`, `extinction_run_I.png`.

![law](law_rv.png)

![run](extinction_run_I.png)

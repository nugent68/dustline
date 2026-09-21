# Example: ZTF20abgaovd (mid-latitude SN Ia sightline, measured R_V)

Sightline of the ZTF DR2 SN Ia ZTF20abgaovd (z = 0.075): RA 251.192, Dec −1.324
(l 15.9°, b +27.2°), **30′ radius**. Law mode with R_V free: 16,558 Gaia sources, 2,629
XP stars fitted on the v0.7.0 stack (UV–IR cache, rebuilt dwarf + giant template
corrections, sub-grid refinement, closure-corrected law); DESI DR1
log g/[Fe/H] priors for 242 stars, APOGEE DR17 for 49 (external T_eff check); PS1 +
2MASS + AllWISE + GALEX NUV (158). Reference foreground: SFD E(B−V) = 0.178 → A_V 0.55
(0.47 with the S&F11 rescaling); Edenhofer+2023 3D map 0.55 to 1.25 kpc.

```bash
dustline run 251.192 -1.324 --radius 30 --min-av 0.5 -o extinction_I.csv --plots .
```

## Result (v0.7.0)

- **R_V = 3.50 ± 0.84** (median ± star-to-star MAD, **1,394 stars** with A_V ≥ 0.5; per-star
  errors 0.23–0.82, ~±0.03 statistical on the median). Uncorrected median 3.32; raw by
  A_V bin: 3.14 (0.5–0.7, N 833), 3.55 (0.7–1.0, N 471), 4.00 (> 1, N 90) — the
  T_eff–A_V–R_V error ellipse of the per-star fit (see `docs/method.md`), not dust.
  Edge pile-up at R_V 2.3 / 5.55 is down to 1.4 % (sub-grid refinement).
- **Column**: A_V = 0.57 ± 0.18 (MAD) beyond 1 kpc, now on the Edenhofer+2023 map (0.55)
  and the unrescaled SFD (0.55); A_I(D) = 0.34 at 0.5 kpc, 0.32–0.36 from 1 to 2 kpc,
  0.42 at 4 kpc, 0.48–0.53 at 6–7.5 kpc.
- **APOGEE check** (46 stars): T_eff(fit) − T_eff(ASPCAP) = −82 K overall (MAD 121);
  −108 K for the 26 dwarfs, −50 K for the 20 giants.
- v0.5/v0.6 (grid without the band corrections at A_V > 0): R_V 3.55 ± 0.85 from 754
  stars, column 0.44 — the fix raised every A_V by ~0.1 here, which is why twice as many
  stars now pass the A_V ≥ 0.5 cut.

## Systematics (this is a first look)

- The DESI T_eff label is not used as a prior (it varies by ±200 K against the colour
  scale with S/N); T_eff comes from XP + photometry + the parallax radius prior, with
  DESI/APOGEE log g and [Fe/H]. A ~65–110 K scale tension with APOGEE remains, worth
  ~0.05–0.1 in A_V and ~0.1–0.2 in R_V.
- Since v0.6.0 the template corrections cover giants too (calibrated on 1,482 nearby
  APOGEE giants with T_eff locked to ASPCAP; closure 0.00 ± 0.01 in the red-clump
  range). But a *free* giant fit in a reddened field still lands ~120 K cooler than
  ASPCAP: giants have no parallax–luminosity lever against the T_eff–A_V degeneracy, so
  the far giants (D > 4 kpc, where A(D) still rises) need an external T_eff (APOGEE here
  covers 20). The dwarfs, which carry the law and the column, are on the APOGEE scale.
- The mean R_V of the G23 family over 0.5–1 mag of extinction is what is measured. With
  the v0.7.0 closure calibration the budget is ±0.05 (closure) ± 0.05 (NIR zero points),
  but at A_V ≈ 0.6 the per-star scatter (MAD 0.84) and the 1/A_V template residual are
  both at their largest; quote R_V = 3.5 ± 0.03 (stat) ± 0.15 (sys).

## Files

`extinction_I.csv`, `law.json` (law, ratios, `apogee_check`, offsets), `law_rv.png`,
`extinction_run_I.png`.

![law](law_rv.png)

![run](extinction_run_I.png)

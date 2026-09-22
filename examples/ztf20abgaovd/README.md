# Example: ZTF20abgaovd (mid-latitude SN Ia sightline, measured R_V)

Sightline of the ZTF DR2 SN Ia ZTF20abgaovd (z = 0.075): RA 251.192, Dec −1.324
(l 15.9°, b +27.2°), **30′ radius**. Law mode with R_V free: 16,558 Gaia sources, 2,629
XP stars fitted on the v0.8.0 stack (UV–IR cache, label-locked per-model template
corrections `tca602f749`, sub-grid refinement, closure-corrected law); DESI DR1
log g/[Fe/H] priors for 242 stars, APOGEE DR17 for 46 (external T_eff check); PS1 +
2MASS + AllWISE + GALEX NUV (158). Reference foreground: SFD E(B−V) = 0.178 → A_V 0.55
(0.47 with the S&F11 rescaling); Edenhofer+2023 3D map 0.55 to 1.25 kpc.

```bash
dustline run 251.192 -1.324 --radius 30 --min-av 0.5 -o extinction_I.csv --plots .
```

## Result (v0.8.0)

- **R_V = 3.46 ± 0.79** (median ± star-to-star MAD, **1,559 stars** with A_V ≥ 0.5; per-star
  errors 0.24–0.81, ~±0.03 statistical on the median). Uncorrected median 3.35 (dwarfs
  3.32, N 1,410; giants 3.49, N 149); the closure correction adds +0.12 — a quarter of the
  sample is on cool templates (< 5000 K, raw 3.15). Raw by A_V bin: 3.18 (0.5–0.8,
  N 1,164), 3.86 (0.8–1.1, N 328), 4.42 (1.1–1.5, N 62) — the T_eff–A_V–R_V error ellipse
  of the per-star fit, not dust (see `docs/method.md`; the injection test reproduces it on
  the calibrators). Edge pile-up at R_V 2.3 / 5.55 is 1.5 % (sub-grid refinement).
- **Column**: A_V = 0.60 ± 0.19 (MAD) beyond 1 kpc, against the Edenhofer+2023 map (0.55)
  and the unrescaled SFD (0.55); A_I(D) = 0.35 at 0.5 kpc, 0.33–0.40 from 1 to 2.5 kpc,
  0.41 at 3–4 kpc, 0.50–0.54 at 6–7.5 kpc.
- **APOGEE check** (46 stars): T_eff(fit) − T_eff(ASPCAP) = −50 K overall (MAD 124); on the
  quality-cut subset −88 K for the 21 dwarfs and −25 K for the 15 giants — the giants were
  −119 K on v0.6.0 and −50 K on v0.7.0, so the label-locked per-model corrections took most
  of the giant T_eff offset out.
- v0.7.0 (coarse-binned corrections, larger closure): R_V 3.50 ± 0.84 from 1,394 stars
  [raw 3.32], column 0.57. v0.5/v0.6 (grid without the band corrections at A_V > 0): R_V
  3.55 ± 0.85 from 754 stars, column 0.44 — the fix raised every A_V by ~0.1 here, which is
  why twice as many stars now pass the A_V ≥ 0.5 cut. The field R_V is stable to ±0.05
  across the three calibrations; the sample has doubled.

## Systematics (this is a first look)

- The DESI T_eff label is not used as a prior (it varies by ±200 K against the colour
  scale with S/N); T_eff comes from XP + photometry + the parallax radius prior, with
  DESI/APOGEE log g and [Fe/H]. The residual ~50–90 K scale tension with APOGEE is worth
  ~0.05 in A_V; v0.8.0's injection controls showed it is *not* a bias term in R_V (locking
  T_eff on the calibrators removes the A_V bias and leaves R_V unchanged — `docs/method.md`).
- Since v0.6.0 the template corrections cover giants too (calibrated on 1,482 nearby
  APOGEE giants with T_eff locked to ASPCAP; injection closure now R_V 2.98 for 3.05,
  ΔA_V −0.014). Giants have no parallax–luminosity lever against the T_eff–A_V degeneracy,
  so the far giants (D > 4 kpc, where A(D) still rises) are the weakest part of the run;
  the dwarfs, which carry the law and the column, are within ~90 K of the APOGEE scale.
- The mean R_V of the G23 family over 0.5–1 mag of extinction is what is measured. With
  the v0.8.0 closure calibration the budget is ±0.05 (closure) ± 0.05 (NIR zero points),
  but at A_V ≈ 0.6 the per-star scatter (MAD 0.79) and the 1/A_V template residual are
  both at their largest; quote R_V = 3.5 ± 0.03 (stat) ± 0.15 (sys).

## Files

`extinction_I.csv`, `law.json` (law, ratios, `apogee_check`, offsets), `law_rv.png`,
`extinction_run_I.png`.

![law](law_rv.png)

![run](extinction_run_I.png)

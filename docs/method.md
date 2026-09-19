# dustline: method notes

## What it measures

For a sightline (RA, Dec), dustline fits every Gaia DR3 star that has an XP
spectrum (G ≲ 17.65) with PHOENIX NewEra model atmospheres reddened by a
Gordon et al. (2023, G23) extinction curve:

    F_obs(λ) = C · F_NewEra(T_eff, log g, [M/H]=0; λ) · 10^(−0.4 · A_V · ext(λ; R_V))

- The XP spectrum (336–1020 nm, R ≈ 50) constrains the SED shape.
- Broadband photometry (PS1 grizy + 2MASS JHKs by default; DECaPS grizY, DECaLS,
  VVV JHKs where available; or your own calibrated photometry) extends the
  wavelength lever to ~2.2 µm — **without it R_V is biased high and T_eff cool**;
  the fit refuses to run XP-only.
- The flux scale C = (R/D)² is profiled analytically; a MIST-isochrone radius
  prior with D = 1/ϖ (Gaia parallax, Lindegren+2021 zero point −0.017 mas)
  breaks the T_eff–A_V degeneracy of hot stars.
- Grid: NewEra [M/H] = 0, T_eff ≥ 3200 K, all log g × A_V 0–8 (0.1) ×
  R_V 2.3–5.55 (0.25).

Per star: (T_eff, log g, A_V, R_V, band extinctions A_b, D). Ensemble:

- **The law**: median R_V and star-to-star MAD over the well-reddened stars
  (A_V ≥ 2 by default), plus per-star band ratios A_b/A_V.
- **The run**: running median (±25 % window in D, ≥ 12 stars) and 16/84 %
  envelope of A_V over the parallax stars (plx S/N > 5), 0.1 kpc steps.
- **Any output filter**: A_X(D) = ratio × A_V(D), where the ratio is either the
  measured per-star ratio (bands in the fit) or the band-integrated G23 ratio
  at the measured R_V through the requested curve on a reference red-giant SED.

## Photometric zero points

Survey zero points are tied to the Gaia XP flux scale by iterating: fit → per-band
median (observed − synthetic) offset → refit. Optical offsets iterate to
convergence (< 0.02 mag); **NIR offsets are frozen after the first pass** — a
common J/H/Ks shift is degenerate with the fit and iterating it drifts R_V by
+0.04/pass (a ±0.05 systematic recorded in the result).

## Known systematics and caveats

- XP blue end (< 500 nm) of faint (G > 16), heavily reddened stars is
  calibration-limited; residual stacks in the research work showed −8/−11/+10 %
  features at 500/550/640 nm for A_V ≥ 3.5 stars. Masking < 500 nm moved R_V by
  +0.05 on one sightline (robust); masking < 600 nm removes the R_V lever.
- Solar metallicity is assumed. Freeing [M/H] lets the extra models absorb the
  XP/NewEra residuals (Teff +225 K, R_V +0.35 — a failure mode, not a
  systematic; validated against APOGEE giants: T_eff bias +50 K at [M/H] ≳ −0.3).
- The G23 curve family is adequate to < 3 % over 450–950 nm on the validated
  sightlines; the R_V you get is the G23 R_V.
- The run covers the parallax-supported range (typically ≲ 5–7 kpc at plx
  S/N > 5 in the plane). On bulge sightlines (|l| < 20°, |b| < 10°) it is
  extended by a linear bridge to the red-clump column at D_RC and held flat
  behind it (rows flagged `bridged`); elsewhere no extrapolation is applied
  beyond the last bin.
- Quality cuts throughout: ruwe < 1.4, ipd_frac_multi_peak ≤ 10, χ²/n < 2.5,
  ≥ 4 photometric bands, A_V < 7.9 (grid edge).

## Validated against

- OGLE-2024-BLG-0669 (l 13.2, b −1.6; PS1 + 2MASS, 9′): R_V = 2.88, MAD 0.53
  (986 stars); red-clump column A_V = 3.28 ± 0.61 at D_RC = 6.1 kpc;
  A_I = 0.42/0.89/1.09/1.21/1.32 at 1/2/3/4/5 kpc, bridged to 1.92 beyond
  the clump. The packaged run is in examples/ob240669 and pinned by
  tests/test_regression_ob240669.py (runs when the sightline workspace is in
  the local cache).

## References

- NewEra models: Hauschildt et al. 2025 (PHOENIX/NewEra)
- Extinction: Gordon et al. 2023 (G23), via `dust_extinction`
- Gaia XP: Montegriffo et al. 2023; De Angeli et al. 2023 (gaiaxpy)
- MIST isochrones: Choi et al. 2016, v1.2
- Parallax zero point: Lindegren et al. 2021

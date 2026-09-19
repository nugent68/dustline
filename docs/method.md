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
  The exception is a field with DESI MWS parameters (high latitude): there each
  star with a DESI spectrum carries a Gaussian prior on (T_eff, log g, [Fe/H])
  and the [M/H] axis (−0.5/0/+0.5) is opened — the prior, not the XP residuals,
  then sets the metallicity.
- The G23 curve family is adequate to < 3 % over 450–950 nm on the validated
  sightlines; the R_V you get is the G23 R_V.
- The run covers the parallax-supported range (typically ≲ 5–7 kpc at plx
  S/N > 5 in the plane). On bulge sightlines (|l| < 20°, |b| < 10°) it is
  extended by a linear bridge to the red-clump column at D_RC and held flat
  behind it (rows flagged `bridged`); elsewhere no extrapolation is applied
  beyond the last bin.
- Quality cuts throughout: ruwe < 1.4, ipd_frac_multi_peak ≤ 10, χ²/n < 2.5,
  ≥ 4 photometric bands, A_V < 7.9 (grid edge).

## Column mode (|b| > 30°)

No A_V ≥ 2 stars, so R_V is not measurable: G23 at R_V = 3.1 is assumed (and
held fixed in the per-star fits) and the product is the total foreground
column, the median A_V of the plx S/N > 5 F/G stars (T_eff ≥ 5500 K) beyond
0.5 kpc (behind all the local dust), with a bootstrap error, the star-to-star
MAD, the cool-star / spectroscopic-prior / G-magnitude splits, and a cell map
of the field. What the COSMOS test (examples/cosmos, SFD A_V ≈ 0.05) showed:
- K/M dwarfs (mostly G > 15.5) return A_V ≈ 0.2–0.3: a cool-template / faint-XP
  systematic, not dust — hence the T_eff cut on the headline number.
- Without DESI priors the F/G fits run ~155 K hot and absorb it as +0.08 mag of
  A_V (the T_eff–A_V degeneracy at R ≈ 50); with the priors the F/G column is
  0.086 ± 0.006 (MAD 0.048) instead of 0.170 ± 0.011 (MAD 0.105).
- The residual ~0.03–0.04 mag above SFD for the DESI-anchored F/G stars (the
  brightest, G < 14, give 0.04) is the method's floor at low A_V, set by the
  XP/NewEra residuals once T_eff is pinned. The zero-point offsets converged at
  |Δ| ≤ 0.016 mag, so the screen / zero-point degeneracy (`freeze_offsets`
  control) is not what sets it.
- Since v0.4.0 the model cache (built from the NewEra HSR files' low-sampling
  spectra) spans 900 Å–6 µm and [M/H] −2..+0.5: the metal-poor halo stars no
  longer sit at a grid edge, GALEX FUV/NUV and WISE W1/W2 can enter the fit.
  Stars without a spectroscopic prior are confined to [M/H] ±0.5 (free
  metallicity absorbs residuals). Net effect on the COSMOS F/G column:
  0.086 → 0.072 (SFD 0.059); the DESI-anchored F/G stars stay at 0.093.
- GALEX: the NewEra NUV flux is too bright by +0.15/+0.18/+0.04 mag at
  5000–5500/5500–6000/6000–6500 K (a per-T_eff zero point is iterated), and the
  residual star-to-star scatter is 0.15–0.23 mag, i.e. σ(A_V) ≈ 0.065 per
  star through A_NUV/A_V = 2.9 — no gain at A_V ≈ 0.06 (per-star ΔA_V
  0.000 ± 0.006 with vs without NUV). WISE W1 agrees with 2MASS and the models
  to a millimag and moves the column by 0.002.

## Validated against

- OGLE-2024-BLG-0669 (l 13.2, b −1.6; PS1 + 2MASS, 9′): R_V = 2.84, MAD 0.48
  (946 stars); red-clump column A_V = 3.33 ± 0.62 at D_RC = 6.1 kpc;
  A_I = 0.41/0.84/1.00/1.17/1.31 at 1/2/3/4/5 kpc, bridged to 1.94 beyond
  the clump (v0.4.0 cache; the 2500–25000 Å cache gave 2.88 ± 0.53). The packaged run is in examples/ob240669 and pinned by
  tests/test_regression_ob240669.py (runs when the sightline workspace is in
  the local cache).
- COSMOS (l 237, b +42; column mode, 30′): F/G foreground column A_V = 0.072 ±
  0.009 vs SFD 0.059 (examples/cosmos, tests/test_regression_cosmos.py).

## References

- NewEra models: Hauschildt et al. 2025 (PHOENIX/NewEra)
- Extinction: Gordon et al. 2023 (G23), via `dust_extinction`
- Gaia XP: Montegriffo et al. 2023; De Angeli et al. 2023 (gaiaxpy)
- MIST isochrones: Choi et al. 2016, v1.2
- Parallax zero point: Lindegren et al. 2021

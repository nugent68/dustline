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

## Empirical template corrections and the T_eff lock (v0.5.0)

The K-dwarf systematic of the COSMOS test (A_V ≈ 0.2 for T_eff < 5000 K against a
0.06 foreground) decomposed into three pieces, none of them dust or XP
calibration (the XP-independent colour excesses E(g−Ks), E(g−W1), E(g−z) show
the same 0.16–0.28 mag):
1. The DESI RVSpecFit T_eff label is not a stable coordinate at the 100 K level:
   it is 50–120 K too hot against the empirical (IRFM) Mamajek BP−RP locus, the
   offset depends on S/N and brightness, and calibrators of one DESI label
   differ by 0.1 mag in BP−RP between G 12 and G 14.
2. NewEra's dwarf SEDs at a given T_eff are too blue: 18 % too bright at 400 nm
   at 4000–4500 K, 7–11 % at 4500–5250 K, 2–6 % at 5250–6000 K.
3. A soft T_eff prior (σ 50–100 K) cannot hold either: with A_V free the fit
   keeps the label and pays with A_V (~0.05 per 100 K); with A_V held at 0 it
   drifts 150–200 K cooler instead — so corrections measured one way did not
   transfer to the other.

The design that closes (`dustline.calib`, `tools/build_template_corrections.py`):
- **T_eff lock**: every dwarf's T_eff is the empirical Mamajek locus value of its
  dereddened Gaia BP−RP (`calib.teff_from_bprp`, packaged locus, a NewEra-based
  [Fe/H] term of 0.03–0.06 mag/dex — metal-poor stars read 150–200 K too hot
  otherwise), locked (σ = 1 K) to the nearest 100 K grid node. In column mode
  the dereddening uses the previous pass's A_V (0 in pass 1; converges in two
  passes at A_V ~ 0.05); DESI supplies log g and [Fe/H] (field median when a
  star has none). Stars without a BP−RP T_eff keep the old behaviour.
- **Calibrators**: 2,451 DESI DR1 dwarfs (log g > 4, S/N > 20) within 250 pc at
  |b| > 40° with XP, PS1, 2MASS, AllWISE, GALEX. They are *not* dust-free (a
  100–250 pc star at |b| > 40° sits at z = 80–190 pc, behind most of a
  ~100 pc dust layer): each is dereddened with the Edenhofer+2023 3D map
  (median A_V 0.026, A_V = 2.8 E) before an A_V = 0 fit at its locked node.
- **Corrections**: per grid node, [M/H] bin and log g class *of the best-fit
  model*, the median XP obs/model ratio spectrum (unsmoothed, not renormalised)
  and per-band offsets, multiplied into the models of the fit grid
  (`template_corrections.npz`); the band extinctions stay those of the
  uncorrected SED, so the band offsets hold at every A_V (v0.7.0 fixes).
- **A_V grid**: 0.01 steps below 0.5 and down to −0.1, because with the lock a
  0.1 grid quantises the posterior and a grid starting at 0 biases every
  near-zero star positive.

Giants (v0.6.0): the same machinery with 1,482 nearby APOGEE giants (D < 1.2 kpc,
|b| > 30°, map-dereddened, median A_V 0.08), T_eff and log g locked to ASPCAP
(IRFM scale for giants); the table carries a log g class axis. Closure: 0.00 ±
0.01 at 4300–5300 K, map extinction recovered to 0.01–0.02; the sparse ends
(< 4300 K, > 5300 K) are ±0.06. Caveat: a *free* giant fit in a reddened field
still sits ~120 K below ASPCAP (no parallax–luminosity lever for giants), so
far giants of mid-latitude fields need an external T_eff.

Closure on the dwarf calibrators (corrected templates, A_V free): dereddened G/F
stars 0.00, K stars 0.02–0.04; their own map extinction is recovered to 0.01.
COSMOS: F/G column 0.027 ± 0.005 (MAD 0.050, from 0.072/0.065), K/M stars
0.026 (from 0.185), G-magnitude trend 0.01/0.03/0.03/0.05 (from
0.02/0.06/0.10/0.20); the lock alone without corrections gives 0.039 (F/G) and
0.069 (K/M). SFD 0.059: the absolute zero point carries a ±0.03 systematic
(calibrators are bright, solar-metallicity, 100–250 pc stars dereddened with
the 3D map at A_V = 2.8 E; XP is ~0.03 mag too red in g−z against PS1 for faint
red stars). PS1 mean-PSF magnitudes brighter than g 14.5 / r 15 / i 15 / z 14 /
y 13 are saturated and now dropped everywhere. The bulge field is unaffected
(R_V 2.80 ± 0.46 vs 2.84 ± 0.48). `DUSTLINE_TEMPLATE_CORR=none` disables.

## Reddening-injection closure and the R_V systematics (v0.7.0)

The two SN Ia sightlines showed the per-star R_V rising with the fitted A_V and
T_eff (3.1 → 4.1 across A_V quartiles at fixed distance; cool stars 2.8, hot
3.5). `tools/inject_reddening.py` measures this directly on the template
calibrators: the map-dereddened DESI dwarfs and APOGEE giants get a known G23
reddening (A_V = 1, R_V = 3.05 on the grid) in their XP spectra and photometry
and are refitted exactly as a field is (T_eff free, log g / [Fe/H] priors,
template corrections), with the parallax errors inflated to S/N 10 to mimic
1–2 kpc field stars. What it found, in order:

1. **A grid bug** (v0.5–v0.6): `build_grid` measured the band extinctions
   against the *corrected* band fluxes, so the per-band template corrections
   cancelled at every A_V ≠ 0 — a 0.03–0.05 mag NIR correction is a third of
   A_Ks at A_V = 1, and the coolest dwarfs came back at the R_V grid floor
   (2.31). Fixed; grids are rebuilt (cache key `b2`).
2. **The correction table**: the σ = 20 nm smoothing of the median ratio
   spectra erased the TiO-band-scale structure of the cool residuals (−24 % at
   400 nm left in the 3500 K node, which a free A_V read as extra reddening:
   XP-only ΔA_V +0.26), and the table was binned by the calibrator's own
   [Fe/H] / log g while applied by the model's — at 3500 K the −0.5 and 0.0
   templates differ by 30 % in the blue, so the mixed medians fitted neither.
   Now unsmoothed and binned by the best-fit model; per-node closure is exact.
3. **Grid quantisation**: for a high-S/N star the posterior is confined to one
   node (0.25 in R_V, 0.1 in A_V above 0.5) and the 0.1 A_V grid samples the
   diagonal A_V–R_V valley at discrete points (a zig-zag χ² profile), so the
   posterior mean quantised to the nodes and every cool node's median sat at
   exactly 2.80. `fit._refine_star` now evaluates a fine local (R_V, A_V) grid
   (0.025 × 0.01) around the best node of every model carrying the posterior,
   with the profiled scale and the radius prior; the R_V histogram is smooth.
4. **The coupling** (giants, T_eff free): the fit trades +0.45 in R_V and +0.10
   in A_V per +100 K of T_eff error; with T_eff locked the R_V–A_V correlation
   remains — it is the tilt of the per-star error ellipse. Both are symmetric:
   the ensemble median is unbiased by them, and splitting a field by *fitted*
   T_eff selects on the error (the bulge's 3.8 vs 2.55 warm/cool split), so
   no product does that.
5. **What is left**: a template residual, linear (mirroring each star's XP
   about its template mirrors the offset: 2.92 ↔ 3.17), scaling as 1/A_V
   (Δ(1/R_V) = 0.060 / 0.029 / 0.018 at A_V 0.5 / 1 / 2), [Fe/H]-dependent
   within the 0.5 dex model bins, and confined to the cool templates: dwarfs
   < 5100 K read R_V 2.7–2.9 for 3.05 injected, giants < 4500 K 2.8, warmer
   nodes 3.0–3.05 (ΔA_V −0.02 to −0.09, the K dwarfs also 67 K cool). A star
   whose XP is replaced by its own template recovers 3.05 exactly.

The residual is carried as a **closure table** in `template_corrections.npz`
(`rv_closure_k[node, class]` = A_V·median(1/R_V_fit − 1/R_V_inj), `stage
closure`): `ensemble.measure_law` subtracts k/A_V from each star's 1/R_V
(`calib.rv_closure`), moves its band ratios along the law, and reports the
uncorrected median as `rv_raw`. The correction is +0.03 in 1/R_V for the cool
nodes (R_V 2.8 → 3.05 at A_V 1; half that at A_V 2) and < 0.01 for the warm
ones; with the cool stars a quarter of a mid-latitude law sample it moves the
field R_V by +0.05–0.1. Honest budget after it: ±0.05 from the closure
(calibrators are solar-neighbourhood, [Fe/H] −0.5..+0.3), plus the
NIR zero-point freeze ±0.05.

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

- OGLE-2024-BLG-0669 (l 13.2, b −1.6; PS1 + 2MASS, 9′): R_V = 2.76, MAD 0.50
  (857 stars); red-clump column A_V = 3.41 ± 0.64 at D_RC = 6.1 kpc;
  A_I = 0.34/0.81/0.95/1.07/1.15 at 1/2/3/4/5 kpc, bridged to 1.96 beyond
  the clump (v0.6.0; v0.5.0 gave 2.80 ± 0.46, v0.4.0 2.84 ± 0.48, the
  2500–25000 Å cache 2.88 ± 0.53).
- ZTF20abgaovd (l 15.9, b +27.2; mid-latitude SN Ia sightline, 30′): R_V = 3.55,
  MAD 0.85 (754 stars with A_V ≥ 0.5); column A_V = 0.44 ± 0.11 beyond 1 kpc vs
  SFD 0.47–0.55 (examples/ztf20abgaovd). The packaged run is in examples/ob240669 and pinned by
  tests/test_regression_ob240669.py (runs when the sightline workspace is in
  the local cache).
- COSMOS (l 237, b +42; column mode, 30′): F/G foreground column A_V = 0.027 ±
  0.005 (K/M 0.026) vs SFD 0.059, ±0.03 systematic (examples/cosmos,
  tests/test_regression_cosmos.py).

## References

- NewEra models: Hauschildt et al. 2025 (PHOENIX/NewEra)
- Extinction: Gordon et al. 2023 (G23), via `dust_extinction`
- Gaia XP: Montegriffo et al. 2023; De Angeli et al. 2023 (gaiaxpy)
- MIST isochrones: Choi et al. 2016, v1.2
- Parallax zero point: Lindegren et al. 2021

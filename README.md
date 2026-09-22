# dustline

**Measure the extinction law and the 3D dust run along any sightline.**

Give `dustline` an RA and Dec; it returns the extinction A_X(D) as a function of
distance (kpc) in any filter you ask for (Cousins **I** by default), together with the
**measured** total-to-selective extinction ratio R_V for that line of sight — no assumed
law.

It works by fitting Gaia DR3 XP spectra (330–1050 nm, R ≈ 50) of every field star with
PHOENIX NewEra model atmospheres × a G23(R_V) extinction curve, anchored by broadband
photometry (PS1 grizy + 2MASS JHKs by default; DECaPS/DECaLS/VVV where available, or
your own calibrated photometry) and Gaia parallaxes. Per-star (T_eff, A_V, R_V, D)
estimates are combined into the sightline law and a running extinction–distance profile.

## Install

```bash
pip install git+https://github.com/nugent68/dustline
pip install "dustline[plot] @ git+https://github.com/nugent68/dustline"   # + matplotlib for the figures
```

Model assets (the PHOENIX NewEra spectral cache, 214 MB: 4,366 models, 900 Å–6 µm,
[M/H] −2..+0.5; MIST v1.2 isochrones) are downloaded automatically from the GitHub
release on first use and cached in `~/.cache/dustline/`.

## Quickstart

```python
import dustline

sl = dustline.Sightline(ra=275.089125, dec=-18.269389)  # degrees, ICRS (OGLE-2024-BLG-0669)
res = sl.run()          # first run per sightline: ~1-3 h (Gaia XP fetch + fits); cached after

print(res.rv)           # measured R_V: median, scatter (MAD), N stars
tab = res.extinction()  # A_I(D): DataFrame with D_kpc, A_med, A_16, A_84, n_stars
tab = res.extinction("DECam_g")   # any filter with a curve
print(res.law)          # band ratios A_X/A_I under the measured law
res.clump_anchor        # red-clump bulge column {AV_column, D_RC, E_JK, ...} or None
res.plots("I")          # law_rv.png + extinction_run_I.png (needs matplotlib)
```

Command line:

```bash
dustline run 275.089125 -18.269389 --filter I -o extinction.csv
dustline run 275.089125 -18.269389 --phot my_calibrated_phot.csv
dustline run 275.089125 -18.269389 --plots figs/     # also write the two diagnostic figures
```

### Your own photometry

Point the pipeline at a calibrated CSV with columns `ra, dec, mag_<band>, magerr_<band>`
and tell it which filter curve each band is (built-in names or your own curve files):

```python
phot = dustline.UserPhotometry("my_phot.csv", bands={"g": "DECam_g", "Ks": "VISTA_Ks"})
res = dustline.Sightline(ra, dec, photometry=phot).run()
```

## What you get

- **R_V** — median and star-to-star scatter of the well-reddened (A_V ≥ 2, or `--min-av`)
  parallax stars, closure-corrected per star (the reddening-injection test on the
  template calibrators: cool templates read R_V ~0.1–0.25 low at A_V 1; `rv_raw`
  keeps the uncorrected value), with systematic notes (closure ±0.05, NIR zero-point
  sensitivity ±0.05, blue-end XP caveats).
- **A_X(D)** — running median and 16/84 % envelope of the per-star extinctions of the
  Gaia-parallax stars in 0.1 kpc steps, for any filter.
- **Red-clump bridge** — on bulge/long-bar sightlines (|l| < 20°, |b| < 10°) the red
  clump is located in the field's (J−Ks, Ks) CMD (2MASS or VVV), its colour excess is
  turned into an A_V column with the *measured* law ratios, and its distance D_RC comes
  from the dereddened Ks against the NewEra clump model. The run is then extended past
  the last parallax bin by a linear bridge to the clump column at D_RC and held flat
  behind it (rows flagged `bridged`, `n_stars = 0`). The anchor is only accepted when the
  clump window is well populated, E(J−Ks) > 0.15 and 5 < D_RC < 12 kpc; otherwise the
  run simply stops at the last parallax bin.
- **The law** — band-integrated A_X/A_I ratios under the measured R_V for any curve.
- **Figures** — `res.plots()` / `dustline run ... --plots [dir]` write `law_rv.png`
  (per-star R_V vs A_V and the R_V histogram) and `extinction_run_<band>.png` (per-star
  A_X vs D coloured by T_eff, the running median with its 16–84 % band, and the dashed
  clump bridge when present).

## High-latitude fields: column mode

Above |b| ≈ 30° there is no A_V ≥ 2 star to measure R_V with, and all the dust sits within
a few hundred pc. The pipeline then switches to **column mode** (`plan.mode == "column"`):
G23 at R_V = 3.1 is assumed (and held in the per-star fits), the A_V(D) run is still
produced, and the product is the **total foreground column** — the median A_V of the
F/G parallax stars (T_eff ≥ 5500 K) behind the dust (D > 0.5 kpc), its bootstrap error
and star-to-star MAD, the same for the cool stars, the stars with / without a
spectroscopic prior and per G-magnitude bin, and a coarse map of the column across the
field (`res.column`, `column_av.png`). Meant for extragalactic fields, as an independent
check of SFD/Planck-type foregrounds. On the COSMOS test field
([examples/cosmos](examples/cosmos)) the F/G column is 0.023 ± 0.006 (MAD 0.046) and the
K/M dwarfs agree to 0.008 (0.031); SFD gives 0.05–0.06. The absolute zero point carries a
±0.03 systematic from the template calibration — the v0.8.0 injection controls show it is
not a template A_V zero point (−0.007 for pinned F/G stars).

**T_eff lock and template corrections (v0.5.0).** In column mode every dwarf's T_eff is
locked to the empirical (Mamajek) main-sequence locus of its dereddened Gaia BP−RP, with a
[Fe/H] term, and the NewEra templates carry empirical corrections per T_eff node built
from 2,451 nearby DESI×XP dwarfs and 1,482 nearby APOGEE×XP giants (`dustline.calib`;
asset `template_corrections.npz`). Where the DESI DR1 Milky Way Survey covers the field
(Dec > −25, |b| > 15) it supplies log g and [Fe/H] (Data Lab `desi_dr1.mws`, joined on
Gaia `source_id`); its T_eff is *not* used as the prior — it proved unstable at the 100 K
level, and it was the main cause of a 0.2 mag K-dwarf excess on COSMOS. Stars without a
DESI [Fe/H] stay on the [M/H] ±0.5 range. `--no-spectro` drops the DESI priors;
`DUSTLINE_TEMPLATE_CORR=none` drops the corrections; `--freeze-offsets` holds the
photometric zero points at 0 (the control for the screen / zero-point degeneracy).

**GALEX and WISE.** Off the plane (|b| > 20° for GALEX AIS, > 10° for AllWISE) the fit
also takes GALEX FUV/NUV (GUVcat_AIS) and AllWISE W1/W2. The UV is the extinction lever
at low A_V (A_NUV/A_V = 2.9 vs 1.2 in g), but the models' UV flux of cool/active stars is
uncertain at 0.1–0.3 mag, so NUV is kept only for BP−RP < 1.0 (~T_eff > 5300 K), FUV for
BP−RP < 0.6, both with a 0.10 mag systematic. Their zero points are iterated **per T_eff
bin** (the NewEra UV flux is ~0.35 mag too bright for G stars and the bias depends on
T_eff; a global or frozen offset would leak into A_V), so what the UV adds is per-star
precision and field structure, not an independent mean column. W1/W2 carry no dust signal
(A_W1/A_V ≈ 0.06) — they anchor the Rayleigh–Jeans tail, i.e. the flux scale and, through
the radius prior, T_eff, which is what the stars without a DESI prior lack. Both need the
UV–IR model cache (`newera_uvir_cache.npz`, 900 Å–6 µm, [M/H] −2..+0.5, the default
since v0.4.0); with the older 2500–25000 Å cache (`DUSTLINE_MODEL_CACHE=newera_full_cache.npz`)
they are skipped with a message. `--no-uv` / `--no-mir` turn them off.

What the COSMOS test found: WISE agrees with 2MASS and the models to a millimag (W1
offset +0.001) and shifts the column by 0.002; GALEX NUV, after its per-T_eff zero point
(+0.15/+0.18/+0.04 mag at 5000–5500/5500–6000/6000–6500 K: the models are too bright),
leaves a star-to-star residual of 0.15–0.23 mag — σ(A_V) ≈ 0.065 per star through
A_NUV/A_V = 2.9, worse than the 0.035 the XP + DESI fit already reaches. At A_V ≈ 0.06 the
UV therefore changes per-star A_V by 0.000 ± 0.006 and the column not at all. It should
start to pay at A_V ≳ 0.3, where the UV signal (≈ 0.9 mag) dwarfs the model noise.

## Example

[examples/ob240669](examples/ob240669) is a complete run on the OGLE-2024-BLG-0669
sightline (l 13.2°, b −1.6°; PS1 + 2MASS): R_V = 2.82 ± 0.50 from 905 stars, the A_I(D)
table bridged to the red-clump column (A_V = 3.34 at 6.1 kpc), the law JSON and the two
figures. [examples/cosmos](examples/cosmos) is the high-latitude column-mode run, and
[examples/ztf19abqmpti](examples/ztf19abqmpti) (R_V = 3.50 ± 0.45 from 2,229 stars at
A_V ≈ 1) and [examples/ztf20abgaovd](examples/ztf20abgaovd) (3.46 ± 0.79 at A_V ≈ 0.6) are
mid-latitude SN Ia sightlines with a measured R_V.

![extinction run](examples/ob240669/extinction_run_I.png)

## Validation: reddening injection

`tools/inject_reddening.py` reddens the template calibrators (nearby DESI dwarfs and
APOGEE giants, map-dereddened) with a known G23 law and refits them as a field is fitted.
It measured the fit's T_eff–A_V–R_V coupling (+0.45 in R_V and +0.10 in A_V per +100 K of
T_eff error — symmetric, so the ensemble median is unbiased, but never split a field by
*fitted* T_eff), caught a grid bug and two flaws in the correction table (v0.7.0), and
left a per-node closure table that `measure_law` applies. Warm templates (dwarfs ≥ 5100 K,
giants ≥ 4500 K) close to 0.05 in R_V; cool ones needed the correction. The v0.8.0 controls
(A_V = 0 as well as 1, T_eff locked as well as free) traced the residual K-dwarf T_eff
offset to the 0.5 dex [M/H] grid step and showed it biases **A_V only** — locking T_eff
removes the A_V bias and leaves R_V unchanged — so the coupling is not a bias term in the
R_V budget. See [docs/method.md](docs/method.md).

## Caveats

- The first run per sightline fetches XP spectra star-by-star from the Gaia archive
  (~1.8 s each, typically 2,000–6,000 stars) — hours. It is resumable and fully cached.
- Everything rests on Gaia: sightlines need enough XP + parallax stars (crowded
  low-latitude fields are fine; very high latitudes have little dust to measure —
  there the product is the foreground column, see above).
- XP spectra of faint (G > 16), heavily reddened stars carry blue-end calibration
  systematics; the per-star law cut and residual stacking guard against them, but see
  the method paper for the shape caveats.

## Method / citation

The method was developed for microlensing source-distance work, in particular on the
OGLE-2024-BLG-0669 sightline (see the example above). Paper reference to follow.
NewEra models: Hauschildt et al. (2025). Extinction curves via `dust_extinction`
(G23: Gordon et al. 2023).

## License

MIT.

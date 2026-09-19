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
```

Model assets (PHOENIX NewEra spectral cache ~91 MB, MIST v1.2 isochrones) are downloaded
automatically from the GitHub release on first use and cached in `~/.cache/dustline/`.

## Quickstart

```python
import dustline

sl = dustline.Sightline(ra=267.86642, dec=-33.13517)   # degrees, ICRS
res = sl.run()          # first run per sightline: ~1-3 h (Gaia XP fetch + fits); cached after

print(res.rv)           # measured R_V: median, scatter (MAD), N stars
tab = res.extinction()  # A_I(D): DataFrame with D_kpc, A_med, A_16, A_84, n_stars
tab = res.extinction("DECam_g")   # any filter with a curve
print(res.law)          # band ratios A_X/A_I under the measured law
```

Command line:

```bash
dustline run 267.86642 -33.13517 --filter I -o extinction.csv
dustline run 275.089125 -18.269389 --phot my_calibrated_phot.csv
```

### Your own photometry

Point the pipeline at a calibrated CSV with columns `ra, dec, mag_<band>, magerr_<band>`
and tell it which filter curve each band is (built-in names or your own curve files):

```python
phot = dustline.UserPhotometry("my_phot.csv", bands={"g": "DECam_g", "Ks": "VISTA_Ks"})
res = dustline.Sightline(ra, dec, photometry=phot).run()
```

## What you get

- **R_V** — median and star-to-star scatter of the well-reddened (A_V ≥ 2.5) parallax
  stars, with systematic notes (NIR zero-point sensitivity ±0.05, blue-end XP caveats).
- **A_X(D)** — running median and 16/84 % envelope of the per-star extinctions of the
  Gaia-parallax stars in 0.1 kpc steps, for any filter; in bulge sightlines
  (VVV window) the run is bridged to the red-clump bulge column.
- **The law** — band-integrated A_X/A_I ratios under the measured R_V for any curve.

## Caveats

- The first run per sightline fetches XP spectra star-by-star from the Gaia archive
  (~1.8 s each, typically 2,000–6,000 stars) — hours. It is resumable and fully cached.
- Everything rests on Gaia: sightlines need enough XP + parallax stars (crowded
  low-latitude fields are fine; very high latitudes have little dust to measure).
- XP spectra of faint (G > 16), heavily reddened stars carry blue-end calibration
  systematics; the per-star law cut and residual stacking guard against them, but see
  the method paper for the shape caveats.

## Method / citation

The method was developed for microlensing source-distance work on the sightlines
OGLE-2017-BLG-0095 (measured R_V = 3.15 ± 0.23) and OGLE-2024-BLG-0669
(R_V = 2.73 ± 0.34). Paper reference to follow. NewEra models: Hauschildt et al. (2025).
Extinction curves via `dust_extinction` (G23: Gordon et al. 2023).

## License

MIT.

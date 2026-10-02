# DESI DR1 log g and [Fe/H] vs dustline's NewEra fits (2026-10-02)

`tools/desi_label_check.py`, run at NERSC on the DESI-labelled stars of the 20 high-latitude
fields of `benchmarks/column_zero_point` (DESI S/N >= 10, parallax S/N > 5, good fits): each star
refitted with the DESI log g / [Fe/H] priors removed (one pass, the field's zero points, R_V 3.1,
[M/H] free over -2..+0.5), T_eff either locked to the scale-corrected DESI label (`tlock`) or
free (`free`). `label_check_summary.{txt,csv}`.

| | free (5,870 stars) | tlock (5,697) |
|---|---|---|
| log g: fit - DESI | -0.005 (sd 0.242) | +0.010 (sd 0.256) |
| [M/H]_fit - [Fe/H]_DESI | +0.164 (sd 0.287) | +0.123 (sd 0.243) |
| [M/H]_fit - [M/H]_DESI, alpha-corrected (Salaris+93) | +0.016 (sd 0.297) | -0.019 (sd 0.260) |
| [M/H]_fit vs [Fe/H]_DESI slope | 0.72 | 0.87 |

- log g agrees to the fit precision (our errors ~0.25 dex); cool dwarfs (4000-4750 K) +0.24.
- The [Fe/H] offset is the alpha enhancement: NewEra [M/H] is scaled-solar total metallicity, so
  DESI should enter as [Fe/H] + log10(0.638 x 10^[a/Fe] + 0.362). Then the median offset is
  +/- 0.02; residual structure: hot stars (> 5750 K) +0.14..+0.22, metal-poor compression.
- Consequence: `catalogs.desi.shape` feeds DESI [Fe/H] (no alpha term) as the [M/H] prior, pulling
  alpha-enhanced high-latitude stars ~0.14 dex metal-poor.

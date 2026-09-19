# Example: COSMOS (column mode, DESI priors, GALEX + WISE)

Sightline RA 150.12, Dec +2.21 (l 237°, b +42°), **30′ radius**. High latitude, so
column mode: no measurable R_V (G23 at R_V = 3.1 assumed and held in the per-star fits),
and the product is the total foreground A_V. Photometry PS1 grizy + 2MASS JHKs + AllWISE
W1/W2 (2,281 / 1,793 usable) + GALEX GUVcat_AIS NUV/FUV (294 / 17 usable after the
hot-star colour cuts); 592 Gaia XP stars fitted on the v0.4.0 model cache (900 Å–6 µm,
[M/H] −2..+0.5), **422 of them with DESI DR1 MWS T_eff / log g / [Fe/H] priors**.
Reference: the SFD E(B−V) carried by GUVcat, 0.0186 over the field → A_V = 0.059.

```bash
dustline run 150.12 2.21 --radius 30 --ref-av 0.059 -o extinction_I.csv --plots .
dustline run 150.12 2.21 --radius 30 --no-spectro --ref-av 0.059 --plots nospec/   # A/B
```

(~18 min for the XP fetch, 5 min for the 4,366-model grid, ~10 min of fits.)

## Result

| sample (plx S/N > 5, D > 0.5 kpc) | N | A_V median | ± (bootstrap) | MAD |
|---|---|---|---|---|
| **F/G stars, T_eff ≥ 5500 K** (the recommended column) | 149 | **0.072** | 0.009 | 0.065 |
| … with a DESI prior | 99 | 0.092 | 0.005 | 0.042 |
| … with a GALEX NUV point | 101 | 0.064 | 0.006 | 0.059 |
| K/M stars, T_eff < 5500 K | 210 | 0.185 | 0.010 | 0.107 |
| G < 14 | 40 | 0.015 | 0.011 | |
| 14 ≤ G < 15.5 | 89 | 0.064 | 0.011 | |
| 15.5 ≤ G < 16.5 | 100 | 0.105 | 0.008 | |
| G ≥ 16.5 | 130 | 0.202 | 0.005 | |
| all stars | 359 | 0.110 | | |

**A/B without DESI priors** (`--no-spectro`): the F/G column is 0.145 ± 0.015 with MAD
0.106 — twice the value and scatter. The free fit runs hot on F/G stars and compensates
with A_V (correlation of ΔA_V with ΔT_eff: 0.98 in the v0.3 test); the DESI prior breaks
that degeneracy. The template scales agree: median T_eff(DESI) − T_eff(fit) = 33 K.

**What each ingredient did** (clean F/G column): old 2500–25000 Å cache 0.086 → new cache
with [M/H] −2..+0.5 0.074 → + WISE 0.072 → + GALEX 0.072. The DESI-anchored F/G stars
are 0.092–0.095 in every variant; the shift came from the stars without a prior (now
confined to [M/H] ±0.5). WISE W1 agrees with 2MASS and the models to +0.001 mag. GALEX
NUV needs a per-T_eff zero point (+0.12 / +0.18 / +0.06 mag at 5000–5500 / 5500–6000 /
6000–6500 K: the models are too bright) and leaves a 0.15–0.23 mag star-to-star residual,
so it changes per-star A_V by 0.000 ± 0.006 here — no information at A_V ≈ 0.06.

**What limits the column.** The K/M dwarfs (mostly G > 15.5) return A_V ≈ 0.2 against a
0.06 foreground — a cool-template / faint-XP systematic, not dust — so the headline uses
the T_eff ≥ 5500 K stars. Among those, the DESI-anchored value (0.092) sits ~0.03 above
SFD while the brightest stars (G < 14) give 0.015: the floor of the method at low A_V is
**~0.03 mag**, set by the XP/NewEra residuals once T_eff is pinned. The zero-point offsets
converged at |Δ| ≤ 0.01 mag in the optical (2MASS H +0.03, W2 +0.03).

## Files

- `extinction_I.csv` — `res.extinction("I")` (A_I/A_V = 0.604 at R_V 3.1)
- `law.json` — `res.law`: the assumed law, `column` (all the numbers above, plus the 4×4
  cell map), `phot_offsets` (incl. the GALEX per-T_eff table), `sfd_ebv`, `model_cache`
- `column_av.png` — per-star A_V vs D (DESI-prior stars as squares, coloured by T_eff),
  the column and its MAD, the cool-star level, the SFD reference; cell map of the F/G column
- `column_av_nospectro.png` — the same for the `--no-spectro` A/B run
- `extinction_run_I.png` — the A_I(D) running median (flat beyond ~0.5 kpc)

![column](column_av.png)

![column, no DESI priors](column_av_nospectro.png)

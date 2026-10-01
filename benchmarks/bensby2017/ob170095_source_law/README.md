# 0095 source-distance prior vs the source's optical extinction law

Stage A showed the bulge optical law is steeper between V and I than G23 at the R_V that fits I
and the NIR (A_I/E(V-I) 1.23 vs G23 3.17's 1.63), and that the J-Ks anchor's A_I is ~4 % high.
The declens prior fits the source's g r i z Y Kp photometry with A_X/A_i ratios from
`xp_law.json` (G23 at the field R_V 3.17, 6000 K source). Variants of those ratios
(`xp_law_*.json`), with the official dust profile and with its bridged far end (N = 0 rows,
D >= 5.6 kpc) x 0.96 (`xp_dust_run_anchor096.csv`); declens `sed_distance_prior.py
--profile-unit ai --laws measured`, IMF weights:

| source law (g / r / z / Y per A_i) | official profile | far end x 0.96 |
|---|---|---|
| A G23 3.17 (1.877 / 1.333 / 0.769 / 0.665), official | 4.39 [3.41, 6.33] | 4.51 [3.44, 6.53] |
| B 0095 clump-calibrated: E(g-i) x1.128, E(r-i) x1.021, E(i-z) x1.129, E(i-Y) x1.035 (1.989 / 1.340 / 0.739 / 0.654) | 4.30 [3.46, 6.29] | 4.43 [3.49, 6.60] |
| C bulge-ensemble bracket: E(g-i), E(r-i) x1.329 (Nataf A_I/E(V-I)) (2.166 / 1.442 / 0.769 / 0.665) | 4.11 [3.41, 5.65] | 4.19 [3.44, 6.33] |

(P_SEDdust_XP median, 68 % interval, kpc.) The best-supported combination (B + anchor x0.96)
gives 4.43 [3.49, 6.60], indistinguishable from the official 4.39 [3.41, 6.33]. The law
systematics (these and the clump-curve dust profile, 4.66) span 4.1-4.7 kpc in the median,
about a tenth of the prior's 68 % width.

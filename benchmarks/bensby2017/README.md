# Microlensed-dwarf benchmark (Bensby et al. 2017)

Independent truth for the extinction toward microlensing SOURCES: the 91 microlensed bulge
dwarfs/subgiants with high-resolution spectra of Bensby et al. (2017, A&A 605, A89; VizieR
J/A+A/605/A89). Built by `tools/bensby_events.py` (see its docstring for every definition).

## What the truth is

Each star has a dust-free intrinsic colour and absolute magnitude from spectroscopy, and the
microlensing-technique values, which deredden the source relative to its field's red clump
and place it at the clump distance. The differences give, per source,

- `dE_VI`  = E(V-I)_source - E(V-I)_clump
- `dM_I`   = (mu_source - mu_clump) + (A_I,source - A_I,clump)

and, with the clump's A_I, E(V-I) and distance from the OGLE-III map (Nataf et al. 2013),
the source's `EVI_src`, `AI_src` (Cousins I), `mu_src` / `D_src_kpc`, with errors
`e_EVI_src`, `e_AI_src`, `e_mu_src`.

85 of 91 sources have the full truth (5 lack microlensing colours, 3 are outside OGLE-III;
overlapping). The table's own l, b has typos for two events (OGLE-2014-BLG-1122 by 0.19 deg,
confirmed against the OGLE EWS position); RA/Dec is used.

## What it can and cannot test (first look, 2026-09-30)

- The sources are in the bulge: D_src median 8.4 kpc [6.6, 11.7], and A_I,src / A_I,clump
  = 0.97 [0.90, 1.05]. The benchmark tests the FAR end of an extinction-vs-distance run: the
  red-clump anchor and the conversion to the I band. It does not test the near (< 5 kpc) run.
- Per-star truth error ~10 % in E(V-I) (spectroscopic colour + 0.05 mag microlensing colour +
  map scatter; possibly conservative) and ~0.4 mag in mu_src.
- Baselines on the 77 sources with truth and DECaPS-map coverage:
  - DECaPS 3D map (Zucker et al. 2025) E(B-V) at D_src x its median ratio: 9 % scatter in the
    source E(V-I) (median E(V-I)_src / E(B-V)_map = 1.56);
  - "the source sees the clump column" (Nataf): 8 % scatter, -2 % bias.
  Both are at the truth's noise floor: per star the benchmark cannot rank methods better than
  ~8 %; with ~80 sources it measures ensemble BIASES to ~1 %.

## Plan

- Stage A (all fields, no XP): the dustline clump anchor (J-Ks; and the prototype's
  five-estimator column) converted to A_I, vs `AI_src` and vs Nataf's `AI_RC`: bias and
  scatter of the anchor that sets the far end of P(A_I | D).
- Stage B (`stage_b` = True: 15 fields, 3 per A_I quintile, ~48k XP spectra, ~24 h of XP fetch):
  full `dustline run`s; A_I(D_src) from the run + bridge vs truth, next to the DECaPS 3D map;
  the measured law's I-band ratio; and, with the source photometry, the distance prior's
  coverage of `D_src_kpc`.

`n_xp_5arcmin` gives each field's XP cost (~1.8 s per spectrum); MOA-2011-BLG-445 and
MACHO-1999-BLG-022 are 4.4' apart and share a workspace.

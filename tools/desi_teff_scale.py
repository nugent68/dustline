"""DESI DR1 MWS T_eff against the colour scale the template corrections are built on.

The dwarf template corrections place every calibrator at the empirical T_eff of its dereddened
Gaia BP-RP (calib.teff_from_bprp, Mamajek locus + [Fe/H] term). A column-mode fit that locks T_eff
to the raw DESI label (``desi_teff``) therefore borrows the template of a star dT hotter, where
dT = T_DESI - T_colour. On the high-S/N calibrators dT is +120..+165 K for 4500-6000 K dwarfs;
the 20-field column test (benchmarks/column_zero_point) implies less for the fainter field stars,
so dT is measured here as a function of T_DESI and the DESI spectrum S/N.

Sample: DESI DR1 MWS stars (RVSpecFit clean, log g > 3.5) at |b| > 40 deg within 1 kpc (parallax
S/N > 5), BP-RP dereddened by Edenhofer+23 at their distance (A_V = 2.8 E; at |b| > 40 this is
nearly the whole column, and Edenhofer's total agrees with the DESI stellar-reddening map there).

  python tools/desi_teff_scale.py --out DIR          -> desi_teff_scale.json (+ .csv sample)
"""
import argparse
import io
import json
import os
import urllib.parse
import urllib.request
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
TEFF_EDGES = np.arange(4250.0, 7251.0, 250.0)
SNR_EDGES = np.array([3.0, 7.0, 12.0, 20.0, 35.0, 60.0, 100.0, 1e5])
MIN_PER_CELL = 25


def query(subsample: int):
    from dustline.catalogs.datalab import TAP, UA
    cols = ("m.source_id, m.teff, m.logg, m.feh, m.snr_med, g.ra, g.dec, g.b, g.parallax, g.parallax_over_error, "
            "g.phot_g_mean_mag, g.bp_rp")
    adql = (f"SELECT {cols} FROM desi_dr1.mws m JOIN gaia_dr3.gaia_source g ON g.source_id = m.source_id "
            f"WHERE m.rr_spectype = 'STAR' AND m.rvs_warn = 0 AND m.zcat_primary = 't' AND m.logg > 3.5 "
            f"AND m.snr_med > 3 AND m.teff BETWEEN 4250 AND 7250 AND g.parallax > 1.0 AND g.parallax_over_error > 5 "
            f"AND ABS(g.b) > 40 AND g.bp_rp IS NOT NULL AND MOD(m.source_id, {subsample}) = 1")
    q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql, MAXREC=500000))
    txt = urllib.request.urlopen(urllib.request.Request(TAP + "?" + q, headers=UA), timeout=1800).read().decode()
    return pd.read_csv(io.StringIO(txt)).drop_duplicates("source_id")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--subsample", type=int, default=3)
    a = ap.parse_args()
    from dustline import calib
    import calibrator_dust as CD          # Edenhofer reader (same directory)
    d = query(a.subsample)
    print(f"DESI dwarfs/subgiants at |b| > 40, < 1 kpc: {len(d)} (1/{a.subsample} of source_ids)", flush=True)
    q = CD.edenhofer()
    d["av_E"] = CD.eden_av(q, d.ra.values, d.dec.values, 1000.0 / d.parallax.values)
    d["teff_col"] = calib.teff_from_bprp(d.bp_rp.values, d.av_E.values, d.feh.values)
    d["dT"] = d.teff - d.teff_col
    d = d[np.isfinite(d.dT)]
    d.round(4).to_csv(os.path.join(a.out, "desi_teff_scale_sample.csv"), index=False)
    nt, ns = len(TEFF_EDGES) - 1, len(SNR_EDGES) - 1
    med = np.full((nt, ns), np.nan); n = np.zeros((nt, ns), int); err = np.full((nt, ns), np.nan)
    it = np.digitize(d.teff, TEFF_EDGES) - 1; isn = np.digitize(d.snr_med, SNR_EDGES) - 1
    for i in range(nt):
        for j in range(ns):
            m = (it == i) & (isn == j)
            n[i, j] = m.sum()
            if m.sum() >= MIN_PER_CELL:
                v = d.dT.values[m]
                med[i, j] = np.median(v)
                err[i, j] = 1.2533 * 1.4826 * np.median(np.abs(v - med[i, j])) / np.sqrt(m.sum())
    out = dict(teff_edges=TEFF_EDGES.tolist(), snr_edges=SNR_EDGES.tolist(), dT=np.round(med, 1).tolist(),
               dT_err=np.round(err, 1).tolist(), n=n.tolist(), n_total=int(len(d)), subsample=a.subsample,
               definition="T_DESI (RVSpecFit) - T_colour(BP-RP dereddened by Edenhofer+23, calib.teff_from_bprp)",
               sample="DESI DR1 MWS, rvs_warn 0, log g > 3.5, |b| > 40, plx > 1 mas, plx S/N > 5")
    json.dump(out, open(os.path.join(a.out, "desi_teff_scale.json"), "w"), indent=1)
    print("dT = T_DESI - T_colour [K], median per (T_DESI, S/N) cell (N):")
    print("  T_DESI \\ S/N " + "".join(f"{f'{SNR_EDGES[j]:.0f}-{SNR_EDGES[j + 1]:.0f}':>13s}" for j in range(ns)))
    for i in range(nt):
        print(f"  {TEFF_EDGES[i]:.0f}-{TEFF_EDGES[i + 1]:.0f}  " + "".join(
            f"{(f'{med[i, j]:+.0f} ({n[i, j]})' if np.isfinite(med[i, j]) else f'- ({n[i, j]})'):>13s}" for j in range(ns)))


if __name__ == "__main__":
    main()

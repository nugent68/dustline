"""Mirror GALEX GUVcat_AIS (Bianchi+2017; VizieR II/335, which has no flat dump) from the MAST HLSP
CSVs (12 Galactic-latitude slices), as the columns the package queries from VizieR:

  surveys/galex/guvcat_ais/hp32/NNNNN.fits   rows sorted by nside-4096 pixel (_HPX)

HLSP -> VizieR II/335 names: ra->RAJ2000, dec->DEJ2000, fuv_mag->FUVmag, fuv_magerr->e_FUVmag,
nuv_mag->NUVmag, nuv_magerr->e_NUVmag, fuv_artifact->Fafl, nuv_artifact->Nafl, fuv_flags->Fexf,
nuv_flags->Nexf, e_bv->E(B-V).  The HLSP writes missing values as -999 (or -99); they become NaN
(VizieR blanks), otherwise `e_FUVmag < 0.35` would pass them.

  python guvcat.py download     # on a DTN
  python guvcat.py convert      # ~40 GB RAM
  python guvcat.py manifest
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import SURVEYS, Ledger, atomic_path, write_manifest  # noqa: E402
import _healpix as H  # noqa: E402

HLSP = "https://archive.stsci.edu/hlsps/guvcat"
SLICES = ["00n-15n", "15n-30n", "15s-00n", "30n-45n", "30s-15s", "45n-60n", "45s-30s", "60n-75n", "60s-45s",
          "75n-90n", "75s-60s", "90s-75s"]
RAW = SURVEYS / "_staging" / "guvcat"
OUT = SURVEYS / "galex" / "guvcat_ais"
LEDGER = SURVEYS / "_ledger" / "guvcat.jsonl"
MAP = {"ra": ("RAJ2000", "D"), "dec": ("DEJ2000", "D"), "fuv_mag": ("FUVmag", "E"), "fuv_magerr": ("e_FUVmag", "E"),
       "nuv_mag": ("NUVmag", "E"), "nuv_magerr": ("e_NUVmag", "E"), "fuv_artifact": ("Fafl", "J"),
       "nuv_artifact": ("Nafl", "J"), "fuv_flags": ("Fexf", "J"), "nuv_flags": ("Nexf", "J"), "e_bv": ("E(B-V)", "E")}
SENTINEL_COLS = ("fuv_mag", "fuv_magerr", "nuv_mag", "nuv_magerr", "e_bv")


def fname(s):
    return f"hlsp_guvcat_galex_imaging_ais-glat-{s}_fuv-nuv_v1_cat.csv"


def cmd_download(a):
    RAW.mkdir(parents=True, exist_ok=True)
    for s in SLICES:
        subprocess.run(["curl", "-sSf", "--retry", "5", "-C", "-", "-o", str(RAW / fname(s)), f"{HLSP}/{fname(s)}"], check=True)
        print("ok", fname(s), flush=True)


def cmd_convert(a):
    import pyarrow.csv as pacsv
    from astropy.io import fits

    parts = {k: [] for k in MAP}
    n_in = 0
    for s in SLICES:
        t = pacsv.read_csv(RAW / fname(s), convert_options=pacsv.ConvertOptions(include_columns=list(MAP)))
        for k in MAP:
            v = t.column(k).to_numpy(zero_copy_only=False)
            parts[k].append(v)
        n_in += t.num_rows
        print(f"  {s}: {t.num_rows:,}", flush=True)
    d = {k: np.concatenate(v) for k, v in parts.items()}
    for k in SENTINEL_COLS:
        v = d[k].astype("f8")
        d[k] = np.where(v <= -98.9, np.nan, v)
    fine = H.ang2pix_nest(4096, d["ra"], d["dec"])
    o = np.argsort(fine, kind="stable")
    fine = fine[o]
    d = {k: v[o] for k, v in d.items()}
    fp = fine >> 14
    edges = np.searchsorted(fp, np.arange(12 * 32 * 32 + 1))
    led = Ledger(LEDGER)
    for p in range(12 * 32 * 32):
        lo, hi = edges[p], edges[p + 1]
        if hi == lo:
            continue
        cols = [fits.Column(name=MAP[k][0], format=MAP[k][1], array=d[k][lo:hi]) for k in MAP]
        cols.append(fits.Column(name="_HPX", format="K", array=fine[lo:hi]))
        dest = OUT / "hp32" / f"{p:05d}.fits"
        tmp = atomic_path(dest)
        fits.BinTableHDU.from_columns(cols, name="GUVCAT_AIS").writeto(tmp, overwrite=True)
        os.replace(tmp, dest)
        led.add(f"hp{p:05d}", nrows=int(hi - lo))
    print(f"{n_in:,} rows written", flush=True)


def cmd_manifest(a):
    done = Ledger(LEDGER).done()
    table = {
        "status": "partial", "format": "fits", "hdu": 1, "files": "guvcat_ais/hp32/{pix:05d}.fits",
        "partition": "position", "file_nside": 32, "nest": True, "pixels": sorted(int(u[2:]) for u in done),
        "sort_key": {"col": "_HPX", "shift": 0, "pad_arcsec": 0},
        "columns": {v[0]: {"src": v[0]} for v in MAP.values()},
        "nrows_total": int(sum(r["nrows"] for r in done.values())),
        "provenance": {"source": HLSP, "tool": "tools/mirror/guvcat.py", "created": time.strftime("%Y-%m-%d"),
                       "note": "status set to complete by check_table.py after comparison with VizieR II/335"},
    }
    p = write_manifest(SURVEYS / "galex", {"II/335/galex_ais": table}, "galex")
    print(f"{p}: {len(table['pixels'])} pixels, {table['nrows_total']:,} rows")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for c in ("download", "convert", "manifest"):
        sub.add_parser(c)
    a = ap.parse_args()
    {"download": cmd_download, "convert": cmd_convert, "manifest": cmd_manifest}[a.cmd](a)


if __name__ == "__main__":
    main()

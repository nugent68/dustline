"""Mirror the Zucker et al. (2025) DECaPS2 stellar-inference catalog (709 M stars; Data Lab table
decaps_dr2.stellar_inference; Harvard Dataverse doi:10.7910/DVN/K88GFI, 100 FITS batches).

The package reads VVV JHKs from it (mag_6..8 / magerr_6..8, with gaia_id and decaps_fracflux_3).
The batches store the photometry as arrays (mag[13], magerr[13], decaps_fracflux[5]); Data Lab
flattens them to mag_1..mag_13 etc., which the manifest reproduces with per-column array indices.
Kept: ids, astrometry, all photometry, chi2, and the 2.5/16/50/84/97.5 percentiles of distance,
A_V and R_V (float32) - not the posterior samples or the model-grid parameters.

  surveys/zucker25/stellar_inference/hp32/NNNNN.fits   rows sorted by nside-4096 pixel (_HPX)

  python dataverse_download.py doi:10.7910/DVN/K88GFI $SURVEYS/_staging/zucker25   # on a DTN
  python zucker25.py convert        # whole-node memory (~150 GB): gather, sort, write per pixel
  python zucker25.py manifest
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import SURVEYS, Ledger, atomic_path, write_manifest  # noqa: E402
import _healpix as H  # noqa: E402

RAW = SURVEYS / "_staging" / "zucker25"
OUT = SURVEYS / "zucker25" / "stellar_inference"
LEDGER = SURVEYS / "_ledger" / "zucker25.jsonl"
KEEP = {  # batch column -> (FITS format, output dtype)
    "decaps_id": ("K", "i8"), "gaia_id": ("K", "i8"), "parallax": ("E", "f4"), "parallax_error": ("E", "f4"),
    "ra": ("D", "f8"), "dec": ("D", "f8"), "decaps_fracflux": ("5E", "f4"), "mag": ("13E", "f4"),
    "magerr": ("13E", "f4"), "chi2": ("E", "f4"), "dist": ("5E", "f4"), "extinction": ("5E", "f4"), "rv": ("5E", "f4"),
}


FRAG = Path(os.environ.get("PSCRATCH", "/tmp")) / "dustline_mirror_staging" / "zucker25_frag"


def _pass1(name):
    """One batch -> per-pixel fragments (needed columns + fine pixel), sorted by fine pixel."""
    import fitsio
    led = Ledger(LEDGER)
    if f"p1/{name}" in led.done():
        return name, 0
    d = fitsio.read(RAW / name, columns=list(KEEP))
    rec = np.empty(len(d), dtype=[(k.upper(), dt, (int(f[:-1]),)) if f[:-1].isdigit() else (k.upper(), dt)
                                  for k, (f, dt) in KEEP.items()] + [("_HPX", "i8")])
    for k in KEEP:
        rec[k.upper()] = d[k]
    rec["_HPX"] = H.ang2pix_nest(4096, d["ra"], d["dec"])
    del d
    rec = rec[np.argsort(rec["_HPX"], kind="stable")]
    fp = rec["_HPX"] >> 14
    stem = name.replace(".fits.gz", "")
    for p in np.unique(fp):
        lo, hi = np.searchsorted(fp, [p, p + 1])
        dd = FRAG / f"{p:05d}"
        dd.mkdir(parents=True, exist_ok=True)
        tmp = atomic_path(dd / f"{stem}.npy")
        with open(tmp, "wb") as fh:
            np.save(fh, rec[lo:hi])
        os.replace(tmp, dd / f"{stem}.npy")
    led.add(f"p1/{name}", nrows=int(len(rec)))
    return name, len(rec)


def _pass2(pixdir):
    from astropy.io import fits
    p = int(pixdir.name)
    led = Ledger(LEDGER)
    if f"hp{p:05d}" in led.done():
        return p, 0
    rec = np.concatenate([np.load(f) for f in sorted(pixdir.glob("*.npy"))])
    rec = rec[np.argsort(rec["_HPX"], kind="stable")]
    cols = [fits.Column(name=k.upper(), format=KEEP[k][0], array=rec[k.upper()]) for k in KEEP]
    cols.append(fits.Column(name="_HPX", format="K", array=rec["_HPX"]))
    dest = OUT / "hp32" / f"{p:05d}.fits"
    tmp = atomic_path(dest)
    fits.BinTableHDU.from_columns(cols, name="STELLAR_INFERENCE").writeto(tmp, overwrite=True)
    os.replace(tmp, dest)
    led.add(f"hp{p:05d}", nrows=int(len(rec)))
    return p, len(rec)


def cmd_convert(a):
    """Two streaming passes (the 709 M rows do not fit one node's memory with workers):
    1) per batch, in parallel: needed columns -> per-pixel fragments on $PSCRATCH;
    2) per pixel, in parallel: merge its fragments, sort by fine pixel, write the FITS file."""
    from multiprocessing import Pool
    files = json.load(open(RAW / "_files.json"))
    names = sorted(f["filename"] for f in files)
    n = 0
    with Pool(a.nproc) as pool:
        for name, k in pool.imap_unordered(_pass1, names):
            n += k
            print(f"  pass 1 {name}: {k:,}", flush=True)
    done1 = [u for u in Ledger(LEDGER).done() if u.startswith("p1/")]
    if len(done1) != len(names):
        raise SystemExit(f"pass 1 incomplete: {len(done1)}/{len(names)} batches")
    pix = sorted(FRAG.iterdir())
    with Pool(a.nproc2) as pool:
        tot = sum(k for _, k in pool.imap_unordered(_pass2, pix, chunksize=4))
    print(f"pass 2: {len(pix)} pixels, {tot:,} rows written this run", flush=True)


def cmd_manifest(a):
    done = {u: r for u, r in Ledger(LEDGER).done().items() if u.startswith("hp")}
    cols = {"ra": {"src": "RA"}, "dec": {"src": "DEC"}, "gaia_id": {"src": "GAIA_ID"}, "decaps_id": {"src": "DECAPS_ID"},
            "parallax": {"src": "PARALLAX"}, "parallax_error": {"src": "PARALLAX_ERROR"}, "chi2": {"src": "CHI2"}}
    for k in range(1, 14):
        cols[f"mag_{k}"] = {"src": "MAG", "index": k - 1}
        cols[f"magerr_{k}"] = {"src": "MAGERR", "index": k - 1}
    for k in range(1, 6):
        cols[f"decaps_fracflux_{k}"] = {"src": "DECAPS_FRACFLUX", "index": k - 1}
    for q, pct in enumerate(("2p5", "16", "50", "84", "97p5")):
        for c in ("dist", "extinction", "rv"):
            cols[f"{c}_p{pct}"] = {"src": c.upper(), "index": q}
    table = {
        "status": "partial", "format": "fits", "hdu": 1, "files": "stellar_inference/hp32/{pix:05d}.fits",
        "partition": "position", "file_nside": 32, "nest": True, "pixels": sorted(int(u[2:]) for u in done),
        "sort_key": {"col": "_HPX", "shift": 0, "pad_arcsec": 0}, "columns": cols,
        "nrows_total": int(sum(r["nrows"] for r in done.values())), "nrows_expected": 709_129_917,
        "provenance": {"source": "Harvard Dataverse doi:10.7910/DVN/K88GFI (Zucker et al. 2025)",
                       "tool": "tools/mirror/zucker25.py", "created": time.strftime("%Y-%m-%d"),
                       "note": "status set to complete by check_table.py after comparison with Data Lab"},
    }
    p = write_manifest(SURVEYS / "zucker25", {"decaps_dr2.stellar_inference": table}, "zucker25")
    print(f"{p}: {len(table['pixels'])} pixels, {table['nrows_total']:,} rows")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("convert"); c.add_argument("--nproc", type=int, default=12)
    c.add_argument("--nproc2", type=int, default=32)
    sub.add_parser("manifest")
    a = ap.parse_args()
    {"convert": cmd_convert, "manifest": cmd_manifest}[a.cmd](a)


if __name__ == "__main__":
    main()

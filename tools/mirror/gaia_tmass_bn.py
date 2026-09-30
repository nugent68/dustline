"""Mirror the Gaia–2MASS best-neighbour cross-match (gaiadr3.tmass_psc_xsc_best_neighbour; DR3
reuses the EDR3 table, published on the ESA CDN under gedr3/cross_match/tmasspscxsc_best_neighbour).

The package's Gaia cone joins gaia_source -> best neighbour (source_id) -> 2MASS PSC
(designation = original_ext_source_id).  Here the best-neighbour rows are regrouped by Gaia
HEALPix nside-32 pixel (source_id >> 49), the same partition as the cosmo gaia_source copy, so a
cone reads one small file per Gaia pixel:

  surveys/gaia_dr3/tmass_best_neighbour/hp32/NNNNN.fits   (source_id, original_ext_source_id, angular_distance)

  python gaia_tmass_bn.py download      # on a DTN: md5-checked, into STAGING
  python gaia_tmass_bn.py convert       # needs ~40 GB RAM: read all, sort by source_id, write per pixel
  python gaia_tmass_bn.py manifest
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import STAGING, SURVEYS, Ledger, atomic_path, md5sum, write_manifest  # noqa: E402

CDN = "https://cdn.gea.esac.esa.int/Gaia/gedr3/cross_match/tmasspscxsc_best_neighbour"
RAW = STAGING / "gaia_tmass_bn"
OUT = SURVEYS / "gaia_dr3" / "tmass_best_neighbour"
LEDGER = SURVEYS / "_ledger" / "gaia_tmass_bn.jsonl"
SHIFT = 49                      # source_id >> 49 = HEALPix nside-32 nested index
NPIX = 12 * 32 * 32


def md5_list() -> list[tuple[str, str]]:
    p = RAW / "_MD5SUM.txt"
    if not p.exists():
        RAW.mkdir(parents=True, exist_ok=True)
        subprocess.run(["curl", "-sSf", "-o", str(p), f"{CDN}/_MD5SUM.txt"], check=True)
    return [tuple(line.split()[::-1]) for line in p.read_text().splitlines() if line.strip()]


def cmd_download(a):
    for name, md5 in md5_list():
        dest = RAW / name
        for k in range(3):
            if dest.exists() and md5sum(dest) == md5:
                break
            subprocess.run(["curl", "-sSf", "--retry", "5", "-o", str(dest), f"{CDN}/{name}"], check=False)
        else:
            raise SystemExit(f"{name}: md5 mismatch after 3 attempts")
        print(f"  {name} ok", flush=True)


def cmd_convert(a):
    import pyarrow as pa
    import pyarrow.csv as pacsv
    from astropy.io import fits

    files = md5_list()
    sid, des, dist = [], [], []
    n_in = 0
    for name, md5 in files:
        if md5sum(RAW / name) != md5:
            raise SystemExit(f"{name}: md5 mismatch; rerun download")
        with __import__("gzip").open(RAW / name, "rt") as fh:
            nhead = 0
            for line in fh:
                if not line.startswith("#"):
                    break
                nhead += 1
        t = pacsv.read_csv(RAW / name, read_options=pacsv.ReadOptions(skip_rows=nhead),
                           convert_options=pacsv.ConvertOptions(
                               include_columns=["source_id", "original_ext_source_id", "angular_distance"],
                               column_types={"original_ext_source_id": pa.string()}))
        sid.append(t.column("source_id").to_numpy())
        des.append(np.asarray(t.column("original_ext_source_id").to_pylist(), dtype="S17"))
        dist.append(t.column("angular_distance").to_numpy().astype("f4"))
        n_in += t.num_rows
        print(f"  read {name}: {t.num_rows:,}", flush=True)
    sid = np.concatenate(sid); des = np.concatenate(des); dist = np.concatenate(dist)
    o = np.argsort(sid, kind="stable")
    sid, des, dist = sid[o], des[o], dist[o]
    if np.any(np.diff(sid) == 0):
        raise SystemExit("duplicate source_id in the best-neighbour table")
    pix = sid >> SHIFT
    edges = np.searchsorted(pix, np.arange(NPIX + 1))
    led = Ledger(LEDGER)
    done = led.done()
    written = 0
    for p in range(NPIX):
        lo, hi = edges[p], edges[p + 1]
        if hi == lo or f"hp{p:05d}" in done:
            continue
        dest = OUT / "hp32" / f"{p:05d}.fits"
        cols = [fits.Column(name="source_id", format="K", array=sid[lo:hi]),
                fits.Column(name="original_ext_source_id", format="17A", array=des[lo:hi]),
                fits.Column(name="angular_distance", format="E", array=dist[lo:hi])]
        tmp = atomic_path(dest)
        fits.BinTableHDU.from_columns(cols, name="TMASS_BN").writeto(tmp, overwrite=True)
        os.replace(tmp, dest)
        led.add(f"hp{p:05d}", nrows=int(hi - lo))
        written += 1
    print(f"{n_in:,} rows -> {written} pixel files", flush=True)


def cmd_manifest(a):
    done = Ledger(LEDGER).done()
    pixels = sorted(int(u[2:]) for u in done)
    total = sum(r["nrows"] for r in done.values())
    table = {
        "status": "complete", "format": "fits", "hdu": 1, "partition": "source_id",
        "files": "tmass_best_neighbour/hp32/{pix:05d}.fits", "file_nside": 32, "nest": True, "id_shift": SHIFT,
        "pixels": pixels, "id_col": "source_id", "sorted_by": "source_id",
        "columns": {c: {"src": c} for c in ("source_id", "original_ext_source_id", "angular_distance")},
        "nrows_total": total,
        "provenance": {"source": CDN, "tool": "tools/mirror/gaia_tmass_bn.py", "created": time.strftime("%Y-%m-%d"),
                       "note": "Gaia DR3 archive table gaiadr3.tmass_psc_xsc_best_neighbour is the EDR3 cross-match"},
    }
    p = write_manifest(SURVEYS / "gaia_dr3", {"gaiadr3.tmass_psc_xsc_best_neighbour": table}, "gaia_dr3")
    print(f"{p}: {len(pixels)} pixels, {total:,} rows")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for c in ("download", "convert", "manifest"):
        sub.add_parser(c)
    a = ap.parse_args()
    {"download": cmd_download, "convert": cmd_convert, "manifest": cmd_manifest}[a.cmd](a)


if __name__ == "__main__":
    main()

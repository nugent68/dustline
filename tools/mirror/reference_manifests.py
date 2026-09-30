"""Manifests for catalogs NERSC already holds (read in place, not copied).

  python reference_manifests.py gaia_source   # cosmo/data/gaia/dr3/healpix (nside 32 nested by position,
                                              #   rows sorted by source_id; nulls stored as 0, booleans as text)
  python reference_manifests.py twomass [--nproc 32]   # cosmo/data/2mass/healpix: per-file RA/Dec bounding boxes
  python reference_manifests.py twomass-qual CODES.json   # set the 2MASS quality-code -> letter map, mark complete

Each manifest records the column map from the NETWORK table's column names (what the fetchers
ask for) to the file columns, plus the null conventions the reader must undo.
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
from common import SURVEYS, write_manifest  # noqa: E402

GAIA_DIR = "/global/cfs/cdirs/cosmo/data/gaia/dr3/healpix"
TMASS_DIR = "/global/cfs/cdirs/cosmo/data/2mass/healpix"

# checked against the ESA archive for the 27,077 stars of the OGLE-2017-BLG-0095 5' cone (2026-09-30):
# non-null values identical; every ESA NULL of a float column is exactly 0.0 in these files and no
# non-null value is 0.0; integer columns are never NULL there; booleans are 'True'/'False' text
GAIA_VALIDATION = {"against": "ESA gaiadr3.gaia_source", "field": "267.86642,-33.13517 r=5'", "n": 27077,
                   "result": "pass", "date": "2026-09-30"}


def cmd_gaia_source(a):
    import fitsio
    f = fitsio.FITS(f"{GAIA_DIR}/healpix-07182.fits")[1]
    dt = f.read(rows=[0]).dtype
    coords = {"RA", "DEC", "L", "B", "ECL_LON", "ECL_LAT"}
    cols = {}
    for n in dt.names:
        k = dt[n].kind
        spec = {"src": n}
        if k == "f" and n not in coords:
            spec["null_zero"] = True
        if k in "US" and n.startswith(("HAS_", "IN_")):
            spec["bool_text"] = True
        cols[n.lower()] = spec
    pixels = sorted(int(p[8:13]) for p in os.listdir(GAIA_DIR) if p.startswith("healpix-") and p.endswith(".fits"))
    table = {
        "status": "complete", "format": "fits", "hdu": 1, "root": GAIA_DIR, "files": "healpix-{pix:05d}.fits",
        "partition": "position", "file_nside": 32, "nest": True, "pixels": pixels,
        "sort_key": {"col": "SOURCE_ID", "shift": 35, "pad_arcsec": 120,
                     "note": "source_id>>35 is the level-12 pixel of the source's early position (99.4% equal to "
                             "the DR3 position's); pad covers the rest"},
        "ra_col": "ra", "dec_col": "dec", "columns": cols, "validation": GAIA_VALIDATION,
        "provenance": {"source": GAIA_DIR, "tool": "tools/mirror/reference_manifests.py", "created": time.strftime("%Y-%m-%d")},
    }
    p = write_manifest(SURVEYS / "gaia_dr3", {"gaiadr3.gaia_source": table}, "gaia_dr3")
    print(f"{p}: gaia_source, {len(pixels)} files, {len(cols)} columns")


def _bbox(fn):
    import fitsio
    d = fitsio.FITS(fn)[1].read(columns=["RA", "DEC"])
    ra, dec = d["RA"].astype(float), d["DEC"].astype(float)
    if not len(ra):
        return None
    lo, hi = float(ra.min()), float(ra.max())
    if hi - lo > 180:                                    # straddles RA 0: store as a wrapped interval
        r2 = np.where(ra > 180, ra - 360, ra)
        lo, hi = float(r2.min()), float(r2.max())
    return [os.path.basename(fn), lo, hi, float(dec.min()), float(dec.max()), int(len(ra))]


def cmd_twomass(a):
    from multiprocessing import Pool
    files = sorted(str(p) for p in Path(TMASS_DIR).glob("2mass_hp*.fits"))
    with Pool(a.nproc) as pool:
        boxes = [b for b in pool.map(_bbox, files, chunksize=4) if b]
    cols = {"ra": {"src": "RA"}, "dec": {"src": "DEC"}, "designation": {"src": "DESIGNATION"},
            "j_m": {"src": "J_MAG"}, "j_msigcom": {"src": "J_MSIGCOM"}, "h_m": {"src": "H_MAG"},
            "h_msigcom": {"src": "H_MSIGCOM"}, "ks_m": {"src": "K_MAG"}, "ks_msigcom": {"src": "K_MSIGCOM"},
            "k_m": {"src": "K_MAG"}, "k_msigcom": {"src": "K_MSIGCOM"},
            "ph_qual": {"chars": ["J_QUALITY", "H_QUALITY", "K_QUALITY"], "map": {}}}
    base = {"status": "partial", "format": "fits", "hdu": 1, "root": TMASS_DIR, "partition": "bbox",
            "bboxes": [b[:5] for b in boxes], "ra_col": "ra", "dec_col": "dec", "columns": cols,
            "nrows_total": int(sum(b[5] for b in boxes)),
            "provenance": {"source": TMASS_DIR, "tool": "tools/mirror/reference_manifests.py",
                           "created": time.strftime("%Y-%m-%d"),
                           "note": "status stays partial until the quality-code map is set (twomass-qual)"}}
    p = write_manifest(SURVEYS / "twomass", {"gaiadr1.tmass_original_valid": base, "twomass.psc": base}, "twomass")
    print(f"{p}: {len(boxes)} files, {base['nrows_total']:,} rows")


def cmd_twomass_qual(a):
    codes = json.load(open(a.codes))
    p = SURVEYS / "twomass" / "manifest.json"
    d = json.loads(p.read_text())
    for t in d["tables"].values():
        t["columns"]["ph_qual"]["map"] = codes
        t["status"] = "complete"
    write_manifest(SURVEYS / "twomass", d["tables"], "twomass")
    print(f"{p}: quality map {codes}; complete")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("gaia_source")
    t = sub.add_parser("twomass"); t.add_argument("--nproc", type=int, default=16)
    q = sub.add_parser("twomass-qual"); q.add_argument("codes")
    a = ap.parse_args()
    {"gaia_source": cmd_gaia_source, "twomass": cmd_twomass, "twomass-qual": cmd_twomass_qual}[a.cmd](a)


if __name__ == "__main__":
    main()

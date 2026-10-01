"""Extract the columns the package queries from catalogs NERSC already holds into HEALPix nside-32
FITS files named with the NETWORK (Data Lab) column names, so the manifest maps them 1:1:

  apogee   sdss_dr17.apogee2_allstar   <- cosmo sdss/dr17 allStar-dr17-synspec_rev1.fits (733 k rows)
  desi     desi_dr1.mws                <- desi public dr1 VAC mwsall-pix-iron.fits (RVTAB+SPTAB+GAIA, 6.4 M)
  allwise  allwise.source              <- cosmo wise/allwise-catalog part01..50 (747 M rows, 1.5 TB read)

The source rows are wide (allStar ~5 kB, AllWISE ~2 kB), so reading a few columns in place would read
the whole file per query; the extracts are a few GB.  Each writes surveys/<survey>/<table>/hp32/NNNNN.fits
(rows sorted by nside-4096 pixel, _HPX) and a manifest with status "partial" until check_table passes.

  python extract.py apogee|desi|allwise [--nproc 16]
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import SURVEYS, Ledger, atomic_path, write_manifest  # noqa: E402
import _healpix as H  # noqa: E402

ALLSTAR = "/global/cfs/cdirs/cosmo/data/sdss/dr17/apogee/spectro/aspcap/dr17/synspec_rev1/allStar-dr17-synspec_rev1.fits"
MWS = "/global/cfs/cdirs/desi/public/dr1/vac/dr1/mws/iron/v1.0/mwsall-pix-iron.fits"
ALLWISE = "/global/cfs/cdirs/cosmo/data/wise/allwise-catalog"


def _fmt(v: np.ndarray) -> str:
    k = v.dtype.kind
    if k == "b":
        return "L"
    if k in "iu":
        return {1: "I", 2: "I", 4: "J", 8: "K"}[v.dtype.itemsize]
    if k == "f":
        return "D" if v.dtype.itemsize == 8 else "E"
    if k in "SU":
        return f"{max(1, int(np.char.str_len(v.astype(str)).max()) if len(v) else 1)}A"
    raise ValueError(f"unsupported dtype {v.dtype}")


def write_hp32(outdir: Path, cols: dict[str, np.ndarray], ra: np.ndarray, dec: np.ndarray, ledger: Ledger) -> int:
    from astropy.io import fits
    fine = H.ang2pix_nest(4096, ra, dec)
    o = np.argsort(fine, kind="stable")
    fine = fine[o]
    fp = fine >> 14
    edges = np.searchsorted(fp, np.arange(12 * 32 * 32 + 1))
    fmts = {k: _fmt(v) for k, v in cols.items()}
    n = 0
    for p in range(12 * 32 * 32):
        lo, hi = edges[p], edges[p + 1]
        if hi == lo:
            continue
        idx = o[lo:hi]
        fc = [fits.Column(name=k, format=fmts[k], array=(v[idx].astype(str).astype("S") if v.dtype.kind == "U" else v[idx]))
              for k, v in cols.items()]
        fc.append(fits.Column(name="_HPX", format="K", array=fine[lo:hi]))
        dest = outdir / "hp32" / f"{p:05d}.fits"
        tmp = atomic_path(dest)
        fits.BinTableHDU.from_columns(fc).writeto(tmp, overwrite=True)
        os.replace(tmp, dest)
        ledger.add(f"{outdir.name}/hp{p:05d}", nrows=int(hi - lo))
        n += 1
    return n


def manifest(survey: str, table: str, sub: str, cols: dict, ledger: Ledger, ra_col="ra", dec_col="dec", source=""):
    done = {u: r for u, r in ledger.done().items() if u.startswith(f"{sub}/")}
    spec = {"status": "partial", "format": "fits", "hdu": 1, "files": f"{sub}/hp32/{{pix:05d}}.fits",
            "partition": "position", "file_nside": 32, "nest": True,
            "pixels": sorted(int(u.split("/hp")[1]) for u in done), "sort_key": {"col": "_HPX", "shift": 0, "pad_arcsec": 0},
            "ra_col": ra_col, "dec_col": dec_col, "columns": {k: {"src": k} for k in cols},
            "nrows_total": int(sum(r["nrows"] for r in done.values())),
            "provenance": {"source": source, "tool": "tools/mirror/extract.py", "created": time.strftime("%Y-%m-%d"),
                           "note": "status set to complete by check_table.py after comparison with Data Lab"}}
    p = write_manifest(SURVEYS / survey, {table: spec}, survey)
    print(f"{p}: {table} {len(spec['pixels'])} pixels, {spec['nrows_total']:,} rows")


def do_apogee(a):
    import fitsio
    names = ["ra", "dec", "gaiaedr3_source_id", "teff", "logg", "fe_h", "fe_h_err", "m_h", "alpha_m", "snr",
             "starflag", "telescope", "apogee_id"]
    d = fitsio.read(ALLSTAR, ext=1, columns=[n.upper() for n in names])
    cols = {n: np.asarray(d[n.upper()]).astype(d[n.upper()].dtype.newbyteorder("=")) if d[n.upper()].dtype.kind in "iuf"
            else np.char.strip(np.asarray(d[n.upper()]).astype(str)) for n in names}
    ok = np.isfinite(cols["ra"]) & np.isfinite(cols["dec"])
    cols = {k: v[ok] for k, v in cols.items()}
    led = Ledger(SURVEYS / "_ledger" / "apogee.jsonl")
    out = SURVEYS / "apogee_dr17" / "allstar"
    write_hp32(out, cols, cols["ra"], cols["dec"], led)
    manifest("apogee_dr17", "sdss_dr17.apogee2_allstar", "allstar", cols, led, source=ALLSTAR)


def do_desi(a):
    import fitsio
    f = fitsio.FITS(MWS)
    rv = f["RVTAB"].read(columns=["TARGET_RA", "TARGET_DEC", "TEFF", "TEFF_ERR", "LOGG", "LOGG_ERR", "FEH", "FEH_ERR",
                                  "ALPHAFE", "RVS_WARN", "RR_SPECTYPE", "SURVEY", "PROGRAM", "PRIMARY", "TARGETID"])
    sp = f["SPTAB"].read(columns=["SNR_MED", "TARGETID"])
    ga = f["GAIA"].read(columns=["SOURCE_ID"])
    if not np.array_equal(rv["TARGETID"], sp["TARGETID"]):
        raise SystemExit("RVTAB / SPTAB rows not aligned")
    nat = lambda x: np.asarray(x).astype(x.dtype.newbyteorder("=")) if x.dtype.kind in "iuf" else np.char.strip(np.asarray(x).astype(str))  # noqa: E731
    cols = {"source_id": nat(ga["SOURCE_ID"]), "target_ra": nat(rv["TARGET_RA"]), "target_dec": nat(rv["TARGET_DEC"]),
            "teff": nat(rv["TEFF"]), "teff_err": nat(rv["TEFF_ERR"]), "logg": nat(rv["LOGG"]), "logg_err": nat(rv["LOGG_ERR"]),
            "feh": nat(rv["FEH"]), "feh_err": nat(rv["FEH_ERR"]), "alphafe": nat(rv["ALPHAFE"]), "snr_med": nat(sp["SNR_MED"]),
            "survey": nat(rv["SURVEY"]), "program": nat(rv["PROGRAM"]), "rr_spectype": nat(rv["RR_SPECTYPE"]),
            "rvs_warn": nat(rv["RVS_WARN"]), "zcat_primary": np.asarray(rv["PRIMARY"]).astype(bool),
            "targetid": nat(rv["TARGETID"])}
    ok = np.isfinite(cols["target_ra"]) & np.isfinite(cols["target_dec"])
    cols = {k: v[ok] for k, v in cols.items()}
    led = Ledger(SURVEYS / "_ledger" / "desi_mws.jsonl")
    out = SURVEYS / "desi_dr1" / "mws"
    write_hp32(out, cols, cols["target_ra"], cols["target_dec"], led)
    manifest("desi_dr1", "desi_dr1.mws", "mws", cols, led, ra_col="target_ra", dec_col="target_dec", source=MWS)


def _allwise_part(fn):
    import fitsio
    d = fitsio.read(fn, ext=1, columns=["RA", "DEC", "W1MPRO", "W1SIGMPRO", "W2MPRO", "W2SIGMPRO", "CC_FLAGS", "EXT_FLG"])
    return {"ra": d["RA"].astype("f8"), "dec": d["DEC"].astype("f8"), "w1mpro": d["W1MPRO"].astype("f4"),
            "w1sigmpro": d["W1SIGMPRO"].astype("f4"), "w2mpro": d["W2MPRO"].astype("f4"), "w2sigmpro": d["W2SIGMPRO"].astype("f4"),
            "cc_flags": np.char.strip(np.asarray(d["CC_FLAGS"]).astype(str)), "ext_flg": d["EXT_FLG"].astype("i2")}


def do_allwise(a):
    from multiprocessing import Pool
    parts = sorted(str(p) for p in Path(ALLWISE).glob("wise-allwise-cat-part[0-9][0-9].fits"))
    with Pool(a.nproc) as pool:
        chunks = pool.map(_allwise_part, parts, chunksize=1)
    cols = {k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]}
    del chunks
    led = Ledger(SURVEYS / "_ledger" / "allwise.jsonl")
    out = SURVEYS / "allwise" / "source"
    write_hp32(out, cols, cols["ra"], cols["dec"], led)
    manifest("allwise", "allwise.source", "source", cols, led, source=ALLWISE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=["apogee", "desi", "allwise"])
    ap.add_argument("--nproc", type=int, default=16)
    a = ap.parse_args()
    {"apogee": do_apogee, "desi": do_desi, "allwise": do_allwise}[a.which](a)


if __name__ == "__main__":
    main()

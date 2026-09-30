"""Mirror Gaia DR3 XP continuous spectra (xp_continuous_mean_spectrum) from the ESA CDN.

The CDN holds 3,386 gzipped ECSV files (3.61 TB), each one contiguous HEALPix level-8 range of
source_id (the file name XpContinuousMeanSpectrum_AAAAAA-BBBBBB.csv.gz gives the range;
level-8 index = source_id >> 43).  The package only ever looks XP spectra up by source_id, so
each input file is converted on its own into one FITS binary table covering the same range,
sorted by source_id:

  surveys/gaia_dr3/xp_continuous/XpContinuousMeanSpectrum_AAAAAA-BBBBBB.fits

Columns: every column of the DataLink XP_CONTINUOUS RAW table.  Coefficients are kept as float64
(the ESA text carries doubles), coefficient errors and correlations as float32 (their declared
subtype), so the rows can be written back exactly as the DataLink CSV the package caches.
has_bp / has_rp flag a missing band (its scalar fields are then written back empty).

  python gaia_xp.py list                       # fetch _MD5SUM.txt (the unit list)
  python gaia_xp.py run [--task-id i --ntasks n] [--limit k]   # download, verify md5, convert, verify, delete raw
  python gaia_xp.py manifest                   # ledger -> surveys/gaia_dr3/manifest.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import STAGING, SURVEYS, Ledger, atomic_path, claim, md5sum, my_units, release, write_manifest  # noqa: E402

CDN = "https://cdn.gea.esac.esa.int/Gaia/gdr3/Spectroscopy/xp_continuous_mean_spectrum"
OUT = SURVEYS / "gaia_dr3" / "xp_continuous"
LEDGER = SURVEYS / "_ledger" / "gaia_xp.jsonl"
NROWS_DR3 = 219_197_643          # Gaia DR3 sources with XP continuous spectra
NCOEF, NCORR = 55, 55 * 54 // 2

SCALARS = {  # name: (numpy dtype, FITS format)
    "source_id": ("i8", "K"), "solution_id": ("i8", "K"),
}
for _b in ("bp", "rp"):
    SCALARS.update({
        f"{_b}_basis_function_id": ("i2", "I"), f"{_b}_degrees_of_freedom": ("i2", "I"),
        f"{_b}_n_parameters": ("i2", "I"), f"{_b}_n_measurements": ("i2", "I"),
        f"{_b}_n_rejected_measurements": ("i2", "I"), f"{_b}_standard_deviation": ("f4", "E"),
        f"{_b}_chi_squared": ("f4", "E"), f"{_b}_n_relevant_bases": ("i2", "I"),
        f"{_b}_relative_shrinking": ("f4", "E"),
    })
ARRAYS = {}
for _b in ("bp", "rp"):
    ARRAYS.update({f"{_b}_coefficients": ("f8", NCOEF, "D"), f"{_b}_coefficient_errors": ("f4", NCOEF, "E"),
                   f"{_b}_coefficient_correlations": ("f4", NCORR, "E")})
# the DataLink/ECSV column order
ORDER = ["source_id", "solution_id"] + [
    f"{b}_{c}" for b in ("bp", "rp") for c in (
        "basis_function_id", "degrees_of_freedom", "n_parameters", "n_measurements", "n_rejected_measurements",
        "standard_deviation", "chi_squared", "coefficients", "coefficient_errors", "coefficient_correlations",
        "n_relevant_bases", "relative_shrinking")]


def units() -> list[tuple[str, str]]:
    p = OUT / "_MD5SUM.txt"
    if not p.exists():
        raise SystemExit(f"{p} missing: run `gaia_xp.py list` first")
    out = []
    for line in p.read_text().split("\n"):
        if line.strip():
            md5, name = line.split()
            out.append((name, md5))
    return out


def hpx8_range(name: str) -> tuple[int, int]:
    a, b = name.split("_")[1].split(".")[0].split("-")
    return int(a), int(b)


def download(name: str, md5: str, dest: Path, tries: int = 3) -> None:
    for k in range(tries):
        if not (dest.exists() and md5sum(dest) == md5):
            subprocess.run(["curl", "-sSf", "--retry", "5", "--retry-delay", "20", "-C", "-", "-o", str(dest),
                            f"{CDN}/{name}"], check=False)
        if dest.exists() and md5sum(dest) == md5:
            return
        print(f"  {name}: md5 mismatch after attempt {k + 1}; refetching", flush=True)
        dest.unlink(missing_ok=True)
        time.sleep(30 * (k + 1))
    raise RuntimeError(f"{name}: could not download a file matching md5 {md5}")


def _parse_arrays(strs: list, length: int, dtype: str) -> tuple[np.ndarray, np.ndarray]:
    """'[a,b,...]' strings (None/'' for a null array; 'null' elements allowed) -> (n, length) array
    and a has-value mask."""
    n = len(strs)
    out = np.full((n, length), np.nan, dtype=dtype)
    have = np.array([bool(s) and s != "null" for s in strs])
    idx = np.flatnonzero(have)
    for k0 in range(0, len(idx), 5000):
        sel = idx[k0:k0 + 5000]
        joined = ",".join(strs[i][1:-1] for i in sel).replace("null", "nan")
        vals = np.fromstring(joined, sep=",", dtype=np.float64)
        if vals.size != len(sel) * length:
            # rare ragged row: parse one by one to find it
            for i in sel:
                v = np.fromstring(strs[i][1:-1].replace("null", "nan"), sep=",", dtype=np.float64)
                if v.size != length:
                    raise ValueError(f"array of length {v.size} (expected {length}) in row {i}")
                out[i] = v
            continue
        out[sel] = vals.reshape(len(sel), length)
    return out, have


def convert(src: Path, dest: Path, name: str, md5: str) -> dict:
    import pyarrow.csv as pacsv
    from astropy.io import fits

    with gzip.open(src, "rt") as fh:            # ECSV: the YAML header lines start with '#'
        nhead = 0
        for line in fh:
            if not line.startswith("#"):
                break
            nhead += 1
    ropt = pacsv.ReadOptions(skip_rows=nhead, block_size=1 << 26)
    copt = pacsv.ConvertOptions(column_types={c: __import__("pyarrow").string() for c in ARRAYS},
                                strings_can_be_null=True)
    tab = pacsv.read_csv(src, read_options=ropt, convert_options=copt)
    if tab.column_names != ORDER:
        raise ValueError(f"{name}: unexpected columns {tab.column_names}")
    n = tab.num_rows
    sid = tab.column("source_id").to_numpy()
    order = np.argsort(sid, kind="stable")
    sid = sid[order]
    if np.any(np.diff(sid) == 0):
        raise ValueError(f"{name}: duplicate source_id")
    lo, hi = hpx8_range(name)
    h8 = sid >> 43
    if h8.min() < lo or h8.max() > hi:
        raise ValueError(f"{name}: source_id outside its level-8 range {lo}-{hi}: {h8.min()}-{h8.max()}")
    cols, have = [], {}
    for c, (dt, fmt) in SCALARS.items():
        a = tab.column(c)
        null = a.is_null().to_numpy(zero_copy_only=False)
        v = a.to_numpy(zero_copy_only=False)
        fill = -1 if dt.startswith("i") else np.nan
        v = np.where(null, fill, v).astype(dt)[order]
        cols.append(fits.Column(name=c, format=fmt, array=v))
    for c, (dt, length, fmt) in ARRAYS.items():
        strs = tab.column(c).to_pylist()
        arr, h = _parse_arrays(strs, length, dt)
        cols.append(fits.Column(name=c, format=f"{length}{fmt}", array=arr[order]))
        have[c] = h[order]
    for b in ("bp", "rp"):
        h = have[f"{b}_coefficients"]
        cols.append(fits.Column(name=f"has_{b}", format="L", array=h))
    hdu = fits.BinTableHDU.from_columns(cols, name="XP_CONTINUOUS")
    hdr = hdu.header
    hdr["SRCFILE"] = name; hdr["SRCMD5"] = md5; hdr["HPX8LO"] = lo; hdr["HPX8HI"] = hi
    hdr["SIDMIN"] = int(sid[0]); hdr["SIDMAX"] = int(sid[-1]); hdr["ORIGIN"] = "dustline tools/mirror/gaia_xp.py"
    tmp = atomic_path(dest)
    hdu.writeto(tmp, overwrite=True)
    # verify what landed on disk before publishing it
    with fits.open(tmp, memmap=True) as f:
        d = f[1].data
        if len(d) != n or int(d["source_id"][0]) != int(sid[0]) or int(d["source_id"][-1]) != int(sid[-1]):
            raise RuntimeError(f"{name}: written table does not read back")
        k = np.linspace(0, n - 1, min(n, 7)).astype(int)
        for c in ("bp_coefficients", "rp_coefficient_correlations"):
            ref, _ = _parse_arrays([tab.column(c)[int(order[i])].as_py() for i in k], ARRAYS[c][1], ARRAYS[c][0])
            if not np.array_equal(np.asarray(d[c][k]), ref, equal_nan=True):
                raise RuntimeError(f"{name}: {c} read-back mismatch")
    os.replace(tmp, dest)
    return dict(nrows=int(n), sid_min=int(sid[0]), sid_max=int(sid[-1]), hpx8=[lo, hi],
                n_no_bp=int((~have["bp_coefficients"]).sum()), n_no_rp=int((~have["rp_coefficients"]).sum()),
                bytes=dest.stat().st_size)


def cmd_list(a):
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "_MD5SUM.txt"
    subprocess.run(["curl", "-sSf", "-o", str(dest), f"{CDN}/_MD5SUM.txt"], check=True)
    print(f"{dest}: {len(units())} files")


def cmd_run(a):
    done = Ledger(LEDGER).done()
    todo = [u for u in units() if u[0] not in done]
    # split the FIXED unit list (tasks start at different times, so splitting `todo` would overlap)
    mine = [u for u in my_units(units(), a.task_id, a.ntasks) if u[0] not in done]
    if a.reverse:
        mine = mine[::-1]
    if a.limit:
        mine = mine[:a.limit]
    stage = STAGING / "gaia_xp"
    stage.mkdir(parents=True, exist_ok=True)
    print(f"gaia_xp: {len(done)} done, {len(todo)} to do; this task {len(mine)}", flush=True)
    led = Ledger(LEDGER)
    locks = LEDGER.parent / "locks" / "gaia_xp"
    for name, md5 in mine:
        if name in led.done() or not claim(locks, name):
            continue                      # finished or being processed by another runner
        t0 = time.time()
        raw = stage / name
        dest = OUT / name.replace(".csv.gz", ".fits")
        try:
            download(name, md5, raw)
            t1 = time.time()
            rec = convert(raw, dest, name, md5)
        except Exception as e:  # noqa: BLE001 - log and move on; the unit stays undone for the next run
            print(f"  FAILED {name}: {e!r}", flush=True)
            release(locks, name)
            continue
        led.add(name, md5=md5, out=str(dest), dl_s=round(t1 - t0, 1), conv_s=round(time.time() - t1, 1), **rec)
        raw.unlink(missing_ok=True)
        release(locks, name)
        print(f"  {name}: {rec['nrows']} rows, download {t1 - t0:.0f} s, convert {time.time() - t1:.0f} s", flush=True)


def cmd_manifest(a):
    us = units()
    done = Ledger(LEDGER).done()
    recs = [done[n] for n, _ in us if n in done]
    total = sum(r["nrows"] for r in recs)
    complete = len(recs) == len(us) and total == NROWS_DR3
    table = {
        "status": "complete" if complete else "partial",
        "format": "fits", "hdu": 1, "partition": "source_id_range",
        "files": "xp_continuous/{file}", "range_shift": 43,
        "ranges": [[r["hpx8"][0], r["hpx8"][1], Path(r["out"]).name] for r in recs],
        "id_col": "source_id", "sorted_by": "source_id",
        "columns": {c: {"src": c} for c in ORDER},
        "flags": {"has_bp": "bp_* fields null when false", "has_rp": "rp_* fields null when false"},
        "nrows_total": total, "nrows_expected": NROWS_DR3, "n_files": len(recs), "n_files_expected": len(us),
        "provenance": {"source": CDN, "md5_list": "xp_continuous/_MD5SUM.txt", "tool": "tools/mirror/gaia_xp.py",
                       "created": time.strftime("%Y-%m-%d")},
    }
    p = write_manifest(SURVEYS / "gaia_dr3", {"gaiadr3.xp_continuous_mean_spectrum": table}, "gaia_dr3")
    print(f"{p}: {len(recs)}/{len(us)} files, {total:,} rows ({table['status']})")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    r = sub.add_parser("run")
    r.add_argument("--task-id", type=int, default=None); r.add_argument("--ntasks", type=int, default=None)
    r.add_argument("--limit", type=int, default=0)
    r.add_argument("--reverse", action="store_true", help="process this share from the end (a second runner)")
    sub.add_parser("manifest")
    a = ap.parse_args()
    {"list": cmd_list, "run": cmd_run, "manifest": cmd_manifest}[a.cmd](a)


if __name__ == "__main__":
    main()

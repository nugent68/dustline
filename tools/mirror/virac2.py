"""Mirror VVV/VVVX VIRAC2 sources (ESO catalogue VVVX_VIRAC_V2_SOURCES, 545 M rows; Smith et al.)
through the ESO catalogue TAP service (tap_cat; async jobs, FITS output, <= 15 M rows per job).

One job per HEALPix nside-32 pixel over the VVV/VVVX area: an RA/Dec box around the pixel, cut
to the pixel locally, sorted by nside-4096 pixel and written as

  surveys/vvv/virac2/hp32/NNNNN.fits

Columns: identifiers, astrometry, per-band mean magnitudes / scatter / epoch counts and the Ks
variability summary needed to judge a mean magnitude.  Not a default input of the package: the
VVV photometry the fits use comes from decaps_dr2.stellar_inference; VIRAC2 is the switchable
alternative (DUSTLINE_VVV=virac2).

  python virac2.py plan                 # the pixel list (Galactic footprint, generous)
  python virac2.py run [--task-id i --ntasks n] [--pixels 7182 ...]
  python virac2.py manifest
"""

from __future__ import annotations

import argparse
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import SURVEYS, Ledger, atomic_path, claim, my_units, release, write_manifest  # noqa: E402
import _healpix as H  # noqa: E402

TAP = "https://archive.eso.org/tap_cat"
TABLE = "VVVX_VIRAC_V2_SOURCES"
OUT = SURVEYS / "vvv" / "virac2"
LEDGER = SURVEYS / "_ledger" / "virac2.jsonl"
PLAN = SURVEYS / "_ledger" / "virac2_pixels.json"
UA = {"User-Agent": "dustline-mirror"}
COLS = ["sourceid", "ra", "de", "ra_error", "de_error", "ref_epoch", "parallax", "parallax_error", "pmra", "pmra_error",
        "pmde", "pmde_error", "chisq", "uwe", "duplicate", "astfit_params"] + \
       [f"phot_{b}_{q}" for b in ("z", "y", "j", "h", "ks") for q in ("mean_mag", "std_mag", "n_epochs")] + \
       [f"{b}_n_det" for b in ("z", "y", "j", "h", "ks")] + ["ks_mad", "ks_med_err"]


def gal2eq(l, b):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(l=np.atleast_1d(l) * u.deg, b=np.atleast_1d(b) * u.deg, frame="galactic").icrs
    return c.ra.deg, c.dec.deg


def cmd_plan(a):
    """nside-32 pixels whose centre lies in a generous VVV+VVVX envelope (-130 < l < 20, |b| < 16)."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    ra, dec = H.pix2ang_nest(32, np.arange(12 * 32 * 32))
    g = SkyCoord(ra=ra * u.deg, dec=dec * u.deg).galactic
    lw = (g.l.deg + 180) % 360 - 180
    sel = (lw > -130) & (lw < 20) & (np.abs(g.b.deg) < 16)
    pix = np.flatnonzero(sel).tolist()
    PLAN.parent.mkdir(parents=True, exist_ok=True)
    PLAN.write_text(json.dumps(pix))
    print(f"{len(pix)} pixels -> {PLAN}")


def _box(pix: int):
    """RA/Dec box enclosing the pixel (corners sampled from its nside-4096 children)."""
    kids = (np.int64(pix) << 14) + np.arange(0, 1 << 14, 37, dtype=np.int64)
    ra, dec = H.pix2ang_nest(4096, kids)
    pad = 0.02
    if ra.max() - ra.min() > 180:
        ra = np.where(ra > 180, ra - 360, ra)
    return ra.min() - pad, ra.max() + pad, dec.min() - pad, dec.max() + pad


def _async(adql: str, dest, poll=10, timeout=4 * 3600):
    data = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="fits", MAXREC=15_000_000,
                                       PHASE="RUN", EXECUTIONDURATION=3600, QUERY=adql)).encode()
    req = urllib.request.Request(f"{TAP}/async", data=data, headers=UA)
    resp = urllib.request.urlopen(req, timeout=300)
    job = resp.geturl()
    t0 = time.time()
    while True:
        phase = urllib.request.urlopen(urllib.request.Request(f"{job}/phase", headers=UA), timeout=120).read().decode().strip()
        if phase == "COMPLETED":
            break
        if phase in ("ERROR", "ABORTED"):
            err = urllib.request.urlopen(urllib.request.Request(f"{job}/error", headers=UA), timeout=120).read().decode()[:500]
            raise RuntimeError(f"ESO job {phase}: {err}")
        if time.time() - t0 > timeout:
            raise RuntimeError("ESO job timed out")
        time.sleep(poll)
    with urllib.request.urlopen(urllib.request.Request(f"{job}/results/result", headers=UA), timeout=3600) as r, open(dest, "wb") as fh:
        while True:
            b = r.read(1 << 24)
            if not b:
                break
            fh.write(b)
    return job


def do_pixel(pix: int) -> dict:
    from astropy.io import fits
    r0, r1, d0, d1 = _box(pix)
    cond = f"ra BETWEEN {r0 % 360} AND {r1 % 360}" if r0 >= 0 and r1 < 360 else \
        f"(ra >= {r0 % 360} OR ra <= {r1 % 360})"
    adql = f"SELECT {', '.join(COLS)} FROM {TABLE} WHERE {cond} AND de BETWEEN {d0} AND {d1}"
    stage = SURVEYS / "_staging" / "virac2"
    stage.mkdir(parents=True, exist_ok=True)
    raw = stage / f"{pix:05d}.fits"
    t0 = time.time()
    _async(adql, raw)
    t_q = time.time() - t0
    with fits.open(raw) as h:
        d = h[1].data
        n_box = len(d)
        if n_box >= 15_000_000:
            raise RuntimeError(f"pixel {pix}: result hit the 15 M row cap")
        fine = H.ang2pix_nest(4096, np.asarray(d["ra"], float), np.asarray(d["de"], float))
        keep = (fine >> 14) == pix
        o = np.argsort(fine[keep], kind="stable")
        sub = d[keep][o]
        cols = [fits.Column(name=c.name, format=c.format, array=sub[c.name]) for c in h[1].columns]
        cols.append(fits.Column(name="_HPX", format="K", array=fine[keep][o]))
    n = int(keep.sum())
    if n:
        dest = OUT / "hp32" / f"{pix:05d}.fits"
        tmp = atomic_path(dest)
        fits.BinTableHDU.from_columns(cols, name="VIRAC2").writeto(tmp, overwrite=True)
        os.replace(tmp, dest)
    raw.unlink(missing_ok=True)
    return dict(nrows=n, n_box=n_box, query_s=round(t_q, 1))


def cmd_run(a):
    pix = a.pixels or json.loads(PLAN.read_text())
    led = Ledger(LEDGER)
    locks = LEDGER.parent / "locks" / "virac2"
    for p in my_units(pix, a.task_id, a.ntasks):
        u = f"hp{p:05d}"
        if u in led.done() or not claim(locks, u):
            continue
        try:
            rec = do_pixel(int(p))
        except Exception as e:  # noqa: BLE001
            print(f"  FAILED {u}: {e!r}", flush=True)
            release(locks, u)
            continue
        led.add(u, **rec)
        release(locks, u)
        print(f"  {u}: {rec}", flush=True)


def cmd_manifest(a):
    done = Ledger(LEDGER).done()
    plan = json.loads(PLAN.read_text())
    cols = {c: {"src": c} for c in COLS}
    table = {
        "status": "partial", "format": "fits", "hdu": 1, "files": "virac2/hp32/{pix:05d}.fits", "partition": "position",
        "file_nside": 32, "nest": True, "pixels": sorted(int(u[2:]) for u, r in done.items() if r["nrows"] > 0),
        "sort_key": {"col": "_HPX", "shift": 0, "pad_arcsec": 0}, "ra_col": "ra", "dec_col": "de", "columns": cols,
        "nrows_total": int(sum(r["nrows"] for r in done.values())), "nrows_expected": 545_346_533,
        "pixels_done": len(done), "pixels_planned": len(plan),
        "provenance": {"source": f"{TAP} {TABLE}", "tool": "tools/mirror/virac2.py", "created": time.strftime("%Y-%m-%d")},
    }
    p = write_manifest(SURVEYS / "vvv", {"vvv.virac2": table}, "vvv")
    print(f"{p}: {table['pixels_done']}/{table['pixels_planned']} pixels, {table['nrows_total']:,} rows")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    r = sub.add_parser("run")
    r.add_argument("--task-id", type=int, default=None); r.add_argument("--ntasks", type=int, default=None)
    r.add_argument("--pixels", type=int, nargs="*", default=None)
    sub.add_parser("manifest")
    a = ap.parse_args()
    {"plan": cmd_plan, "run": cmd_run, "manifest": cmd_manifest}[a.cmd](a)


if __name__ == "__main__":
    main()

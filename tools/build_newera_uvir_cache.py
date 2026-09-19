"""Build the NewEra UV-to-mid-IR model cache (``newera_uvir_cache.npz``) from the
per-model HSR HDF5 files of the PHOENIX/1D NewEra grid (Hauschildt et al. 2025,
FDR record 16738), for GALEX FUV/NUV and WISE W1/W2 leverage and the metal-poor
extension of the fit grid.

Each HSR file (~110 MB) carries a ``PHOENIX_SPECTRUM_LSR`` group: the spectrum
sampled at 0.1 A from 10 A to 9.5 mm as log10 F_lambda (erg s^-1 cm^-2 cm^-1).
We keep that dataset only, rebinned to the cache grid

    900-25000 A at 2 A   (GALEX FUV/NUV, Gaia XP, optical, JHKs)
  25000-60000 A at 20 A  (WISE W1/W2),

converted to the cache convention (surface flux, W m^-2 nm^-1: 10^(log F - 10)),
and delete the HDF5.  The sub-grid: alpha = 0, Teff >= 3200 K, log g 0-6 in 0.5
steps, [M/H] -2.0 .. +0.5 (4,366 models, ~490 GB of downloads, ~200 MB out).

Resumable: one .npy per model in the shard directory; rerun to continue.
Meant for a NERSC shared-QOS job (8 download threads; see the sbatch header
in the docstring of main()).  Usage:

  python build_newera_uvir_cache.py --list list_of_available_NewEra_models.txt \
      --shards $SCRATCH/newera_uvir_shards --out newera_uvir_cache.npz [--threads 8]
  python build_newera_uvir_cache.py ... --assemble      # shards -> npz
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

W_OPT = np.arange(900.0, 25000.0, 2.0)        # bin lower edges
W_MIR = np.arange(25000.0, 60000.0, 20.0)
WAVE = np.concatenate([W_OPT + 1.0, W_MIR + 10.0])   # bin centres
TEFF_MIN = 3200.0
LOGG = tuple(x / 2 for x in range(0, 13))
MH = (-2.0, -1.5, -1.0, -0.5, 0.0, 0.5)
PAT = re.compile(r"lte(\d{5})([-+]\d\.\d\d)([-+]\d\.\d)(?:\.alpha=([-+]\d\.\d))?\.PHOENIX")


def select_models(list_path: str) -> list[tuple[float, float, float, str, str]]:
    """(Teff, logg, MH, filename, url) for the sub-grid, from the FDR model list."""
    out = []
    for line in open(list_path):
        f = line.split()
        if len(f) < 5 or not f[0].isdigit():
            continue
        m = PAT.match(f[1])
        if not m:
            continue
        teff, logg, mh = float(m.group(1)), -float(m.group(2)), float(m.group(3))
        alpha = float(m.group(4) or 0.0)
        if alpha != 0.0 or teff < TEFF_MIN or logg not in LOGG or mh not in MH:
            continue
        out.append((teff, logg, mh + 0.0, f[1], f[4]))
    return sorted(out)


def rebin_lsr(wl: np.ndarray, logf: np.ndarray) -> np.ndarray:
    """Mean of 10^(logF-10) in each cache bin (W m^-2 nm^-1); NaN-free by interpolation."""
    flux = 10.0 ** (logf - 10.0)
    out = np.empty(WAVE.size, np.float32)
    for k, (edges, step, sl) in enumerate(((W_OPT, 2.0, slice(0, W_OPT.size)),
                                           (W_MIR, 20.0, slice(W_OPT.size, WAVE.size)))):
        sel = (wl >= edges[0]) & (wl < edges[-1] + step)
        idx = np.floor((wl[sel] - edges[0]) / step).astype(int)
        ok = (idx >= 0) & (idx < edges.size)
        sums = np.bincount(idx[ok], weights=flux[sel][ok], minlength=edges.size)
        cnts = np.bincount(idx[ok], minlength=edges.size).astype(float)
        prof = np.where(cnts > 0, sums / np.maximum(cnts, 1), np.nan)
        gd = np.isfinite(prof)
        centres = WAVE[sl]
        out[sl] = np.interp(centres, centres[gd], prof[gd]) if gd.any() else 0.0
    return out


def shard_name(teff, logg, mh) -> str:
    return f"lte{teff:05.0f}-{logg:.2f}{mh:+.1f}.npy"


def process(model, shards: str, tmpdir: str, retries: int = 4) -> str:
    import h5py

    teff, logg, mh, fname, url = model
    dest = os.path.join(shards, shard_name(teff, logg, mh))
    if os.path.exists(dest):
        return "cached"
    fd, tmp = tempfile.mkstemp(dir=tmpdir, suffix=".h5")
    os.close(fd)
    try:
        for k in range(retries):
            try:
                req = urllib.request.Request(url + "?download=1", headers={"User-Agent": "dustline"})
                with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as f:
                    while True:
                        b = r.read(1 << 20)
                        if not b:
                            break
                        f.write(b)
                break
            except Exception as e:  # noqa: BLE001
                if k == retries - 1:
                    return f"FAILED download {fname}: {e}"
                time.sleep(10 * (k + 1))
        with h5py.File(tmp, "r") as h:
            g = h["PHOENIX_SPECTRUM_LSR"]
            wl, logf = g["wl"][:], g["fl"][:]
        prof = rebin_lsr(wl, logf)
        np.save(dest + ".part.npy", prof)
        os.replace(dest + ".part.npy", dest)
        return "ok"
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def assemble(models, shards: str, out: str) -> None:
    flux, meta, missing = [], [], 0
    for teff, logg, mh, _, _ in models:
        p = os.path.join(shards, shard_name(teff, logg, mh))
        if not os.path.exists(p):
            missing += 1
            continue
        flux.append(np.load(p))
        meta.append((teff, logg, mh))
    m = np.array(meta)
    np.savez_compressed(out, wave=WAVE, flux=np.array(flux, np.float32), meta=m)
    print(f"cached {len(flux)} models ({missing} missing) -> {out}")
    print(f"Teff {m[:, 0].min():.0f}-{m[:, 0].max():.0f}  logg {m[:, 1].min():.1f}-{m[:, 1].max():.1f}  "
          f"[M/H] {sorted(set(m[:, 2]))}  wave {WAVE.min():.0f}-{WAVE.max():.0f} A ({WAVE.size} pts)")


def main():
    """sbatch example (Perlmutter, shared QOS; downloads are the bottleneck):

        #SBATCH -q shared -C cpu -c 8 -t 08:00:00 -J newera_uvir
        module load python
        python build_newera_uvir_cache.py --list ... --shards $SCRATCH/newera_uvir_shards \\
               --out $SCRATCH/newera_uvir_cache.npz --threads 8 && \\
        python build_newera_uvir_cache.py --list ... --shards ... --out ... --assemble
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", required=True)
    ap.add_argument("--shards", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0, help="only the first N models (testing)")
    ap.add_argument("--assemble", action="store_true", help="only assemble the shards")
    a = ap.parse_args()
    models = select_models(a.list)
    if a.limit:
        models = models[:a.limit]
    print(f"{len(models)} models in the sub-grid", flush=True)
    if a.assemble:
        assemble(models, a.shards, a.out)
        return
    os.makedirs(a.shards, exist_ok=True)
    tmpdir = tempfile.mkdtemp(prefix="newera_", dir=os.environ.get("SCRATCH", None))
    t0, done, fails = time.time(), 0, 0
    with ThreadPoolExecutor(a.threads) as ex:
        futs = {ex.submit(process, m, a.shards, tmpdir): m for m in models}
        for k, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            if r.startswith("FAILED"):
                fails += 1
                print(r, flush=True)
            elif r == "ok":
                done += 1
            if k % 25 == 0 or k == len(models):
                print(f"  {k}/{len(models)} ({done} new, {fails} failed, {time.time() - t0:.0f} s)",
                      flush=True)
    os.rmdir(tmpdir)
    if fails:
        sys.exit(f"{fails} models failed; rerun to retry")
    assemble(models, a.shards, a.out)


if __name__ == "__main__":
    main()

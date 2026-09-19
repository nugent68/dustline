"""One-off: repack the MIST v1.2 basic isochrone .iso files ([Fe/H] -0.5, 0.0, +0.5)
into the single npz release asset that dustline.models.load_mist() reads, and copy
it plus the NewEra cache into the local asset cache for testing.

  python tools/build_mist_npz.py /path/to/MIST_v1.2_vvcrit0.4_basic_isos [out.npz]

Layout of the npz: for each feh tag in (m0.50, p0.00, p0.50):
  <tag>_cols  (ncol,) unicode column names
  <tag>_data  (N, ncol) float64
"""

import glob
import os
import sys

import numpy as np

TAGS = ["m0.50", "p0.00", "p0.50"]
KEEP = ["log10_isochrone_age_yr", "log_Teff", "log_g", "log_R", "phase",
        "initial_mass", "star_mass", "log_L", "EEP"]


def read_iso(path):
    cols, rows = None, []
    for line in open(path):
        if line.startswith("# EEP"):
            cols = line[1:].split()
        elif line.startswith("#") or not line.strip():
            continue
        else:
            rows.append([float(v) for v in line.split()])
    a = np.array(rows)
    keep_idx = [cols.index(c) for c in KEEP if c in cols]
    keep_cols = [cols[i] for i in keep_idx]
    return keep_cols, a[:, keep_idx]


def main():
    src = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else "mist_v1.2_basic.npz"
    payload = {}
    for tag in TAGS:
        ftag = f"feh_{tag[0]}{tag[1:]}"
        files = glob.glob(os.path.join(src, f"*{ftag}*.iso"))
        if not files:
            raise FileNotFoundError(f"{ftag} in {src}")
        cols, data = read_iso(files[0])
        payload[f"{tag}_cols"] = np.array(cols)
        payload[f"{tag}_data"] = data.astype(np.float64)
        print(f"{tag}: {files[0]} -> {data.shape[0]} rows, cols {cols}")
    np.savez_compressed(out, **payload)
    print(f"wrote {out} ({os.path.getsize(out) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()

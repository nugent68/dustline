"""Validate a mirrored table against its network service with the package's own query helpers:
the same box query with DUSTLINE_MIRROR=off (network) and forced local, on several fields.
A pass (same rows, values equal to float precision) marks the table complete in its manifest.

  python check_table.py decaps_dr2.stellar_inference [--fields ra,dec,r ...] [--mark]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
import types
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(os.path.dirname(os.path.abspath(__file__)))
SRC = Path(os.environ.get("DUSTLINE_SRC", HERE.parent.parent / "src"))
if not (SRC / "dustline").exists():
    SRC = HERE.parent / "pkgsrc"                       # the NERSC layout: _tools/{mirror,pkgsrc}


def _load():
    """dustline.catalogs.{_healpix,local,datalab,vizier} without importing the whole package."""
    sys.modules.setdefault("dustline", types.ModuleType("dustline")).__path__ = [str(SRC / "dustline")]
    pkg = types.ModuleType("dustline.catalogs"); pkg.__path__ = [str(SRC / "dustline" / "catalogs")]
    sys.modules["dustline.catalogs"] = pkg
    mods = {}
    for m in ("_healpix", "local", "datalab", "vizier"):
        spec = importlib.util.spec_from_file_location(f"dustline.catalogs.{m}", SRC / "dustline" / "catalogs" / f"{m}.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[f"dustline.catalogs.{m}"] = mod
        spec.loader.exec_module(mod)
        setattr(pkg, m, mod)
        mods[m] = mod
    return mods


# what each table is queried with in the package (columns, helper, coordinate columns, filters)
CASES = {
    "decaps_dr2.stellar_inference": dict(helper="datalab", ra_col="ra", dec_col="dec", where=None,
        columns=["ra", "dec", "gaia_id"] + [f"mag_{k}" for k in range(1, 14)] + [f"magerr_{k}" for k in range(1, 14)]
        + [f"decaps_fracflux_{k}" for k in range(1, 6)],
        fields=[(267.86642, -33.13517, 0.05), (201.0, -62.5, 0.04), (123.0, -35.0, 0.05), (272.0, -20.0, 0.04)]),
    "II/335/galex_ais": dict(helper="vizier", ra_col="RAJ2000", dec_col="DEJ2000", where=None,
        columns=["RAJ2000", "DEJ2000", "FUVmag", "e_FUVmag", "NUVmag", "e_NUVmag", "Fafl", "Nafl", "Fexf", "Nexf", '"E(B-V)"'],
        fields=[(10.0, -45.0, 0.3), (150.12, 2.21, 0.15), (230.0, 60.0, 0.2), (300.0, -60.0, 0.2)]),
}


def _decimals(v: np.ndarray) -> int:
    """Decimals the service printed (e.g. VizieR stores 4 for magnitudes): the most any value uses."""
    s = pd.Series(v[np.isfinite(v)][:2000]).map(lambda x: repr(float(x)))
    d = s.str.split(".").str[1].str.len().fillna(0) if len(s) else pd.Series([12])
    return int(min(d.max(), 12))


def compare(net: pd.DataFrame, loc: pd.DataFrame, key: list[str]) -> dict:
    """Rows paired by nearest position (0.2"), numbers equal within half a unit of the service's
    last printed decimal (the network side may be rounded), strings/flags exactly."""
    from scipy.spatial import cKDTree
    out = {"n_net": len(net), "n_loc": len(loc)}
    if not len(net) and not len(loc):
        return dict(out, ok=True)
    if len(net) != len(loc):
        return dict(out, ok=False, why="row count")
    rc, dc = key
    cd = np.cos(np.radians(float(np.nanmedian(net[dc]))))
    xa = np.c_[net[rc].values * cd, net[dc].values]
    xb = np.c_[loc[rc].values * cd, loc[dc].values]
    dist, j = cKDTree(xb).query(xa, distance_upper_bound=0.2 / 3600)
    if not np.all(np.isfinite(dist)) or len(set(j)) != len(j):
        return dict(out, ok=False, why=f"{int((~np.isfinite(dist)).sum())} rows unmatched within 0.2 arcsec")
    a = net.reset_index(drop=True)
    b = loc.iloc[j].reset_index(drop=True)
    bad, maxd = {}, {}
    for c in a.columns:
        x, y = a[c], b[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            xv, yv = x.astype(float).values, y.astype(float).values
            tol = 0.5 * 10.0 ** -_decimals(xv) + 1e-6 * np.abs(xv)
            with np.errstate(invalid="ignore"):
                ok = (np.isnan(xv) & np.isnan(yv)) | (np.abs(xv - yv) <= tol)
            if np.isfinite(xv - yv).any():
                maxd[c] = float(np.nanmax(np.abs(xv - yv)))
        else:
            ok = (x.isna() & y.isna()) | (x.astype(str).str.strip() == y.astype(str).str.strip())
        if not np.all(ok):
            bad[c] = int((~np.asarray(ok)).sum())
    return dict(out, ok=not bad, bad=bad, max_abs_diff=maxd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table"); ap.add_argument("--mark", action="store_true", help="set status complete on a pass")
    a = ap.parse_args()
    mods = _load()
    L = mods["local"]
    case = CASES[a.table]
    helper = mods[case["helper"]]
    reg = L._registry()
    if a.table not in reg:
        raise SystemExit(f"{a.table} has no manifest under {L.root()}")
    base, spec = reg[a.table]
    t = L.Table(a.table, base, spec)                       # regardless of status
    results = []
    for ra, dec, r in case["fields"]:
        os.environ["DUSTLINE_MIRROR"] = "off"
        t0 = time.time()
        net = helper.box_query(a.table, case["columns"], ra, dec, r, ra_col=case["ra_col"], dec_col=case["dec_col"],
                               where=case["where"])
        tn = time.time() - t0
        t0 = time.time()
        loc = t.box(case["columns"], ra, dec, r, case["ra_col"], case["dec_col"], case["where"])
        tl = time.time() - t0
        res = compare(net, loc, [case["ra_col"], case["dec_col"]])
        res.update(field=[ra, dec, r], t_network_s=round(tn, 1), t_local_s=round(tl, 2))
        print(res, flush=True)
        results.append(res)
    passed = all(r["ok"] for r in results)
    print("PASS" if passed else "FAIL")
    if a.mark:
        mf = Path(base) / "manifest.json"
        d = json.loads(mf.read_text())
        d["tables"][a.table]["validation"] = {"against": case["helper"], "fields": results, "result": "pass" if passed else "fail",
                                              "date": time.strftime("%Y-%m-%d")}
        if passed:
            d["tables"][a.table]["status"] = "complete"
        mf.write_text(json.dumps(d, indent=1, sort_keys=True))
        print(f"{mf}: status {d['tables'][a.table]['status']}")


if __name__ == "__main__":
    main()

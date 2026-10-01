"""Stage B of the microlensed-dwarf benchmark: full dustline runs at every Bensby+2017 event
(benchmarks/bensby2017/events.csv), made to run many fields at once on one node.

Three phases, each resumable (a finished field / grid is skipped):

  data   catalogs, XP spectra (fetch + calibrate) and the photometry of every field, in
         parallel; records each field's band list (stage_b/<name>.data.json)
  grids  every distinct model grid (band list, law g23, solar [M/H]) built ONCE - concurrent
         fields would otherwise race to write the same cached grid file
  fit    Sightline.run() per field in parallel (fits, law, run, clump bridge); result.json and
         the A_I(D) table (stage_b/<name>.extinction_I.csv), wall time per field

  python tools/bensby_stage_b.py data  --jobs 24
  python tools/bensby_stage_b.py grids --jobs 8
  python tools/bensby_stage_b.py fit   --jobs 24 --threads 5 [--only NAME ...] [--max-hours 3.6]

Environment: DUSTLINE_SURVEYS / DUSTLINE_CACHE_DIR / DUSTLINE_MIRROR (the NERSC mirror), set
before calling. Every field uses the official 0095 settings: 5' radius, parallax errors x1.7.
"""
import argparse
import json
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "..", "benchmarks", "bensby2017")
OUT = os.path.join(BENCH, "stage_b")
RADIUS, PLX_INFLATE = 5.0, 1.7


def events(only=None):
    ev = pd.read_csv(os.path.join(BENCH, "events.csv"))
    return ev[ev.name.isin(only)] if only else ev


def sightline(ra, dec, shallow=False):
    """shallow: let the plan pick PS1 over DECaPS (fields at the DECaPS edge where DECaPS returns
    nothing, e.g. OGLE-2014-BLG-1418 at l 5.4, b -4.35: 2MASS-only otherwise)."""
    from dustline.api import Sightline
    return Sightline(ra, dec, RADIUS, plx_inflate=PLX_INFLATE, prefer_deep=not shallow)


def _data(ev):
    from dustline import bands as bands_mod
    from dustline.catalogs import gaia
    t0 = time.time()
    out = dict(name=ev["name"])
    try:
        sl = sightline(ev["ra"], ev["dec"], ev.get("shallow", False))
        if not sl.ws.has("xp_sampled.npz"):
            sl.ws.seed_from_sibling()
        g = gaia.cone(sl.ws)
        if not sl.ws.has("xp_sampled.npz"):
            gaia.fetch_xp(sl.ws, g)
        gaia.calibrate_xp(sl.ws)
        stars, band_list = bands_mod.assemble(sl.ws, g, sl.plan)
        out.update(ok=True, workspace=str(sl.ws.dir), bands=band_list, spectro=sl.plan.spectro,
                   n_gaia=int(len(g)), n_xp=int((g.has_xp_continuous == True).sum()))   # noqa: E712
    except Exception as e:  # noqa: BLE001
        out.update(ok=False, error=f"{type(e).__name__}: {e}", tb=traceback.format_exc()[-2000:])
    out["seconds"] = round(time.time() - t0, 1)
    return out


def _grid(item):
    from dustline import fit as F
    bands, ra, dec = item
    t0 = time.time()
    # the grid lives in the shared asset cache (keyed by bands/law/[M/H]/corrections), not the workspace
    F.build_grid(sightline(ra, dec).ws, list(bands), "g23", mh=F.MH_DEFAULT)
    return dict(bands=list(bands), seconds=round(time.time() - t0, 1))


def _fit(ev, threads):
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(threads)
    t0 = time.time()
    out = dict(name=ev["name"])
    try:
        sl = sightline(ev["ra"], ev["dec"], ev.get("shallow", False))
        res = sl.run()
        res.save(os.path.join(OUT, f"{ev['name']}.extinction_I.csv"), band="I")
        law = res.law
        out.update(ok=True, workspace=str(sl.ws.dir), rv=law.get("rv"), rv_mad=law.get("rv_mad"),
                   n_stars=law.get("n_stars"), mode=law.get("mode", "law"), clump_anchor=law.get("clump_anchor"),
                   ratio_I=float(res.extinction("I").ratio_av.iloc[0]))
    except Exception as e:  # noqa: BLE001
        out.update(ok=False, error=f"{type(e).__name__}: {e}", tb=traceback.format_exc()[-2000:])
    out["seconds"] = round(time.time() - t0, 1)
    return out


def _pool(fn, items, jobs, _tag, *extra):
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(fn, it, *extra): it for it in items}
        for f in as_completed(futs):
            yield f.result()


def cmd_data(a):
    os.makedirs(OUT, exist_ok=True)
    todo = [dict(r._asdict(), shallow=r.name in (a.shallow or [])) for r in events(a.only).itertuples(index=False)
            if a.redo or not os.path.exists(os.path.join(OUT, f"{r.name}.data.json"))]
    print(f"data: {len(todo)} fields", flush=True)
    for r in _pool(_data, todo, a.jobs, "data"):
        json.dump(r, open(os.path.join(OUT, f"{r['name']}.data.json"), "w"), indent=1)
        print(f"  {r['name']:22s} {'ok' if r['ok'] else 'FAIL ' + r.get('error', '')}  {r['seconds']:.0f} s  "
              f"{r.get('n_xp', '')} XP  {r.get('bands', '')}", flush=True)


def cmd_grids(a):
    sets = {}
    for r in events().itertuples(index=False):
        p = os.path.join(OUT, f"{r.name}.data.json")
        if os.path.exists(p):
            d = json.load(open(p))
            if d.get("ok"):
                sets.setdefault(tuple(d["bands"]), []).append((r.name, r.ra, r.dec))
    print(f"grids: {len(sets)} distinct band lists", flush=True)
    for b, names in sets.items():
        print(f"  {len(names):3d} fields: {list(b)}", flush=True)
    for r in _pool(_grid, [(b, v[0][1], v[0][2]) for b, v in sets.items()], a.jobs, "grid"):
        print(f"  built {r['bands']} in {r['seconds']:.0f} s", flush=True)


def cmd_fit(a):
    os.makedirs(OUT, exist_ok=True)
    t_start = time.time()
    ev = events(a.only)
    ready = []
    for r in ev.itertuples(index=False):
        if not a.redo and os.path.exists(os.path.join(OUT, f"{r.name}.fit.json")):
            if json.load(open(os.path.join(OUT, f"{r.name}.fit.json"))).get("ok"):
                continue
        p = os.path.join(OUT, f"{r.name}.data.json")
        if os.path.exists(p) and json.load(open(p)).get("ok"):
            ready.append(dict(r._asdict(), shallow=r.name in (a.shallow or [])))
    # largest fields first (the wall time follows the XP count)
    ready.sort(key=lambda e: -e.get("n_xp_5arcmin", 0))
    print(f"fit: {len(ready)} fields, {a.jobs} at a time x {a.threads} threads", flush=True)
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = {}
        it = iter(ready)
        for _ in range(a.jobs):
            e = next(it, None)
            if e is not None:
                futs[ex.submit(_fit, e, a.threads)] = e
        while futs:
            done = next(as_completed(futs))
            futs.pop(done)
            r = done.result()
            json.dump(r, open(os.path.join(OUT, f"{r['name']}.fit.json"), "w"), indent=1, default=float)
            ca = r.get("clump_anchor") or {}
            print(f"  {r['name']:22s} {'ok' if r['ok'] else 'FAIL ' + r.get('error', '')}  {r['seconds'] / 60:.1f} min  "
                  + (f"R_V {r['rv']:.2f} +/- {r['rv_mad']:.2f} (N {r['n_stars']})  E(J-Ks)_RC {ca.get('E_JK', float('nan')):.3f}" if r["ok"] else ""),
                  flush=True)
            if (time.time() - t_start) / 3600 > a.max_hours:
                print("time budget reached: no new fields started", flush=True)
                continue
            e = next(it, None)
            if e is not None:
                futs[ex.submit(_fit, e, a.threads)] = e


def cmd_summary(a):
    import numpy as np
    ev = events()
    rows = []
    for r in ev.itertuples(index=False):
        pf, pe = os.path.join(OUT, f"{r.name}.fit.json"), os.path.join(OUT, f"{r.name}.extinction_I.csv")
        if not (os.path.exists(pf) and os.path.exists(pe)):
            continue
        f = json.load(open(pf))
        if not f.get("ok"):
            continue
        x = pd.read_csv(pe)
        d = json.load(open(os.path.join(OUT, f"{r.name}.data.json")))
        row = dict(name=r.name, rv=f["rv"], rv_mad=f["rv_mad"], n_law=f["n_stars"], bands="+".join(sorted({b.split("_")[0] for b in d["bands"]})),
                   D_last_measured=float(x.D_kpc[~x.bridged.astype(bool)].max()) if (~x.bridged.astype(bool)).any() else np.nan,
                   minutes=f["seconds"] / 60)
        for D, tag in ((r.D_src_kpc, "src"), (r.D_RC_kpc, "RC")):
            if np.isfinite(D):
                row[f"AI_run_{tag}"] = float(np.interp(D, x.D_kpc, x.A_med))
                row[f"AI16_run_{tag}"] = float(np.interp(D, x.D_kpc, x.A_16))
                row[f"AI84_run_{tag}"] = float(np.interp(D, x.D_kpc, x.A_84))
                row[f"bridged_{tag}"] = bool(np.interp(D, x.D_kpc, x.bridged.astype(float)) > 0.5)
        rows.append(row)
    s = ev.merge(pd.DataFrame(rows), on="name", how="inner")
    s.round(5).to_csv(os.path.join(BENCH, "stage_b_summary.csv"), index=False)
    rs = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))   # noqa: E731
    t = s.AI_src.notna() & s.AI_run_src.notna()
    L = [f"Stage B: {len(s)} fields fitted; with source truth {t.sum()}"]
    L.append(f"  R_V per field: median {s.rv.median():.2f}, 16-84 % {s.rv.quantile(.16):.2f}-{s.rv.quantile(.84):.2f} "
             f"(star-to-star MAD median {s.rv_mad.median():.2f}, N_law median {s.n_law.median():.0f}); last measured bin median "
             f"{s.D_last_measured.median():.1f} kpc; fit wall time median {s.minutes.median():.0f} min")
    q = s.AI_run_src[t] / s.AI_src[t]
    sig_run = 0.5 * (s.AI84_run_src - s.AI16_run_src)[t]
    z = (s.AI_run_src - s.AI_src)[t] / np.sqrt(sig_run ** 2 + s.e_AI_src[t] ** 2)
    L.append(f"  A_I,run(D_src) / A_I,src: median {np.median(q):.3f}, robust sd {rs(q):.3f} (relative {rs(q) / np.median(q):.3f}); "
             f"z = (run - truth) / sqrt(sig_run^2 + sig_truth^2): median {np.median(z):+.2f}, robust sd {rs(z):.2f}, |z|<1 {np.mean(np.abs(z) < 1):.2f}")
    for lab, m in (("source in the measured run (not bridged)", t & ~s.bridged_src.astype(bool)), ("source in the bridged part", t & s.bridged_src.astype(bool))):
        if m.sum():
            qq = (s.AI_run_src / s.AI_src)[m]
            L.append(f"    {lab:42s}: N {m.sum():3d}  ratio median {np.median(qq):.3f}, robust sd {rs(qq):.3f}")
    if "ebv_decaps_Dsrc" in s:
        m = t & np.isfinite(s.ebv_decaps_Dsrc) & (s.ebv_decaps_Dsrc > 0)
        qm = (s.AI_src / s.ebv_decaps_Dsrc)[m]; qr = (s.AI_run_src / s.AI_src)[m]
        L.append(f"  same {m.sum()} sources: DECaPS 3D map relative scatter {rs(qm) / np.median(qm):.3f} (after its median A_I/E {np.median(qm):.3f}); "
                 f"dustline run {rs(qr) / np.median(qr):.3f} (bias {np.median(qr):.3f}); truth error {np.median((s.e_AI_src / s.AI_src)[m]):.3f}")
    for b, g in s[t].groupby("bands"):
        qq = (g.AI_run_src / g.AI_src)
        L.append(f"  bands {b:16s}: N {len(g):3d}  A_I ratio median {np.median(qq):.3f}  R_V median {g.rv.median():.2f}")
    txt = "\n".join(L)
    open(os.path.join(BENCH, "stage_b_summary.txt"), "w").write(txt + "\n")
    print(txt)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("summary")
    for n in ("data", "grids", "fit"):
        p = sub.add_parser(n)
        p.add_argument("--jobs", type=int, default=8)
        p.add_argument("--only", nargs="*")
        p.add_argument("--redo", action="store_true")
        p.add_argument("--shallow", nargs="*", help="fields where the plan should prefer PS1 over DECaPS")
        if n == "fit":
            p.add_argument("--threads", type=int, default=4)
            p.add_argument("--max-hours", type=float, default=99.0, help="start no new field after this")
    a = ap.parse_args()
    {"data": cmd_data, "grids": cmd_grids, "fit": cmd_fit, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())

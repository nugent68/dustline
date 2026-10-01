"""Stage A of the microlensed-dwarf benchmark: the dustline red-clump anchor at every Bensby+2017
event, against the source truth and the OGLE-III (Nataf+2013) / Gonzalez+2012 clump values.

Per event: Gaia DR3 cone + the plan's optical/NIR photometry (bands.assemble, cached in the
event's workspace), clump.find_clump in (J-Ks, Ks) -> E(J-Ks), Ks_RC, D_RC. No XP spectra.
The clump's E(J-Ks) is turned into A_I (Cousins I) and E(V-I) through G23 at several R_V, with
the band ratios of the clump SED (4700 K, log g 2.5).

  run:        python tools/bensby_stage_a.py run  [--jobs 4] [--only NAME ...]
  summarise:  python tools/bensby_stage_a.py summary
"""
import argparse
import json
import os
import sys
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "..", "benchmarks", "bensby2017")
EVENTS = os.path.join(BENCH, "events.csv")
OUTDIR = os.path.join(BENCH, "stage_a")
RVS = (2.5, 2.8, 3.1, 3.4, 3.7)
SED = dict(teff=4700.0, logg=2.5)
PLX_INFLATE = 1.7


def law_ratios(j_band, ks_band):
    from dustline import ensemble
    out = {}
    for rv in RVS:
        r = {b: ensemble.band_ratio_at_rv(b, rv, **SED) for b in (j_band, ks_band, "Cousins_I", "Johnson_V")}
        e_jk = r[j_band] - r[ks_band]
        out[f"{rv:.1f}"] = dict(AI_per_EJK=r["Cousins_I"] / e_jk, EVI_per_EJK=(r["Johnson_V"] - r["Cousins_I"]) / e_jk,
                                AKs_per_EJK=r[ks_band] / e_jk, ratios_av=r)
    return out


def one(ev: dict) -> dict:
    from dustline import bands as bands_mod
    from dustline import clump
    from dustline.api import Sightline
    from dustline.catalogs import gaia as gaia_mod

    out = dict(name=ev["name"], ra=ev["ra"], dec=ev["dec"])
    try:
        sl = Sightline(ev["ra"], ev["dec"], 5.0, plx_inflate=PLX_INFLATE)
        g = gaia_mod.cone(sl.ws)
        stars, band_list = bands_mod.assemble(sl.ws, g, sl.plan)
        j, k = clump.nir_bands(stars, sl.plan.nir)
        lr = law_ratios(j, k)
        law = dict(ratios_av=lr["3.1"]["ratios_av"])
        a = clump.find_clump(stars, law, j, k)
        out.update(workspace=str(sl.ws.dir), optical=sl.plan.optical, nir=sl.plan.nir, nir_used=j.split("_")[0], bands=band_list,
                   n_gaia=int(len(g)), n_jk=int((np.isfinite(stars.get(f"mag_{j}", np.nan)) & np.isfinite(stars.get(f"mag_{k}", np.nan))).sum()),
                   anchor=a, law=lr, ok=a is not None)
    except Exception as e:  # noqa: BLE001
        out.update(ok=False, error=f"{type(e).__name__}: {e}", tb=traceback.format_exc()[-1500:])
    return out


def run(a):
    os.makedirs(OUTDIR, exist_ok=True)
    ev = pd.read_csv(EVENTS)
    if a.only:
        ev = ev[ev.name.isin(a.only)]
    todo = [r._asdict() for r in ev.itertuples(index=False)
            if a.redo or not os.path.exists(os.path.join(OUTDIR, f"{r.name}.json"))]
    print(f"stage A: {len(todo)} events to do ({len(ev) - len(todo)} done)", flush=True)
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(one, e): e["name"] for e in todo}
        for f in as_completed(futs):
            r = f.result()
            json.dump(r, open(os.path.join(OUTDIR, f"{r['name']}.json"), "w"), indent=1, default=float)
            aa = r.get("anchor") or {}
            print(f"  {r['name']:22s} {'ok ' if r.get('ok') else 'FAIL'} "
                  + (f"E(J-Ks) {aa['E_JK']:.3f} Ks_RC {aa['ks_rc']:.2f} N {aa['n_window']} D_RC {aa['D_RC']:.1f}" if aa else r.get("error", "no clump")),
                  flush=True)


def rs(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def summary(a):
    ev = pd.read_csv(EVENTS)
    rows = []
    for r in ev.itertuples(index=False):
        p = os.path.join(OUTDIR, f"{r.name}.json")
        if not os.path.exists(p):
            continue
        j = json.load(open(p))
        an = j.get("anchor") or {}
        row = dict(name=r.name, ok=bool(j.get("ok")), nir=j.get("nir"), nir_used=j.get("nir_used"), E_JK=an.get("E_JK", np.nan), ks_rc=an.get("ks_rc", np.nan),
                   n_window=an.get("n_window", np.nan), D_RC=an.get("D_RC", np.nan), mad_jk=an.get("mad_jk", np.nan))
        for rv, v in (j.get("law") or {}).items():
            row[f"AI_{rv}"] = row["E_JK"] * v["AI_per_EJK"]
            row[f"EVI_{rv}"] = row["E_JK"] * v["EVI_per_EJK"]
            row[f"AIperEJK_{rv}"] = v["AI_per_EJK"]; row[f"EVIperEJK_{rv}"] = v["EVI_per_EJK"]
        rows.append(row)
    s = ev.merge(pd.DataFrame(rows), on="name", how="inner")
    s["EJK_gonzalez"] = s.R_JKVI * s.EVI_RC
    s.round(5).to_csv(os.path.join(BENCH, "stage_a_summary.csv"), index=False)
    ok = s.ok & s.E_JK.notna()
    t = ok & s.AI_src.notna()
    L = [f"Stage A: {len(s)} events run, clump anchor found for {ok.sum()}; with source truth {t.sum()}"]
    L.append(f"  failures / no clump: {', '.join(s.name[~ok])}" if (~ok).any() else "  no failures")
    g = ok & s.EJK_gonzalez.notna()
    d = (s.E_JK - s.EJK_gonzalez)[g]
    L.append(f"\n1. clump colour: dustline E(J-Ks) - Gonzalez+12 E(J-Ks) (= R_JKVI x Nataf E(V-I)): median {np.median(d):+.3f}, "
             f"robust sd {rs(d):.3f}, ratio median {np.median((s.E_JK / s.EJK_gonzalez)[g]):.3f} (N {g.sum()})")
    dd = (s.D_RC - s.D_RC_kpc)[ok & s.D_RC_kpc.notna()]
    L.append(f"   D_RC dustline - Nataf: median {np.median(dd):+.2f} kpc, robust sd {rs(dd):.2f}")
    L.append("\n2. A_I of the clump / of the source from the dustline anchor under G23(R_V), vs truth:")
    L.append(f"   {'R_V':>5s} {'A_I/E(J-Ks)':>11s} {'A_I,anchor / A_I,src':>22s} {'A_I,anchor / A_I,RC(Nataf)':>28s} {'E(V-I) anchor / E(V-I)_src':>28s}")
    for rv in RVS:
        k = f"{rv:.1f}"
        r1 = (s[f"AI_{k}"] / s.AI_src)[t]; r2 = (s[f"AI_{k}"] / s.AI_RC)[ok & s.AI_RC.notna()]; r3 = (s[f"EVI_{k}"] / s.EVI_src)[t]
        L.append(f"   {rv:5.1f} {np.median(s[f'AIperEJK_{k}'][ok]):11.3f} {np.median(r1):14.3f} ± {rs(r1):.3f} "
                 f"{np.median(r2):20.3f} ± {rs(r2):.3f} {np.median(r3):20.3f} ± {rs(r3):.3f}")
    ai_ratio = (s.AI_src / s.E_JK)[t]; vi_ratio = (s.EVI_src / s.E_JK)[t]
    L.append(f"\n3. truth-implied ratios: A_I,src / E(J-Ks)_RC median {np.median(ai_ratio):.3f} (robust sd {rs(ai_ratio):.3f}, "
             f"error of median {1.2533 * rs(ai_ratio) / np.sqrt(t.sum()):.3f}); E(V-I)_src / E(J-Ks)_RC {np.median(vi_ratio):.3f} "
             f"(robust sd {rs(vi_ratio):.3f})")
    x = np.array([np.median(s[f"AIperEJK_{rv:.1f}"][ok]) for rv in RVS]); y = np.array(RVS)
    o = np.argsort(x)
    L.append(f"   -> G23 R_V that reproduces the median A_I/E(J-Ks): {np.interp(np.median(ai_ratio), x[o], y[o]):.2f}; "
             f"E(V-I)/E(J-Ks): {np.interp(np.median(vi_ratio), *[np.array([np.median(s[f'EVIperEJK_{rv:.1f}'][ok]) for rv in RVS])[np.argsort([np.median(s[f'EVIperEJK_{rv:.1f}'][ok]) for rv in RVS])], np.array(RVS)[np.argsort([np.median(s[f'EVIperEJK_{rv:.1f}'][ok]) for rv in RVS])]]):.2f}")
    if "ebv_decaps_Dsrc" in s:
        m = t & np.isfinite(s.ebv_decaps_Dsrc) & (s.ebv_decaps_Dsrc > 0)
        L.append(f"\n4. same sources, the DECaPS 3D map: E(V-I)_src / E(B-V)_map(D_src) median {np.median((s.EVI_src / s.ebv_decaps_Dsrc)[m]):.3f}, "
                 f"relative robust sd {rs((s.EVI_src / s.ebv_decaps_Dsrc)[m]) / np.median((s.EVI_src / s.ebv_decaps_Dsrc)[m]):.3f}; "
                 f"dustline anchor E(J-Ks) route: relative robust sd {rs(vi_ratio[m[t]]) / np.median(vi_ratio[m[t]]):.3f} (N {m.sum()})")
    txt = "\n".join(L)
    open(os.path.join(BENCH, "stage_a_summary.txt"), "w").write(txt + "\n")
    print(txt)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("--jobs", type=int, default=4); r.add_argument("--only", nargs="*")
    r.add_argument("--redo", action="store_true")
    sub.add_parser("summary")
    a = ap.parse_args()
    {"run": run, "summary": summary}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())

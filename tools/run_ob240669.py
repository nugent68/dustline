"""Seed the dustline workspace for OGLE-2024-BLG-0669 (RA 275.089125,
Dec -18.269389) from the research repo's cached Gaia table and calibrated XP
spectra (skipping the ~2 h archive fetch), then run the v0.2 pipeline:
PS1 + 2MASS photometry, per-star fits, law, dust run with the red-clump
bridge, and the two figures.

  .venv/bin/python tools/run_ob240669.py [--max-stars N]
"""

import argparse
import shutil
import sys
from pathlib import Path

RESEARCH = Path.home() / "claude" / "dustline"
RA, DEC = 275.089125, -18.269389
RADIUS = 9.0    # arcmin: the research Gaia pull is a 9' circle


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-stars", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    import dustline

    sl = dustline.Sightline(RA, DEC, radius_arcmin=RADIUS)
    print(f"workspace: {sl.ws.dir}")
    print(f"plan: optical={sl.plan.optical} nir={sl.plan.nir} "
          f"bulge_window={sl.plan.in_bulge_window}")

    # --- seed from the research repo ---
    seeds = {
        "gaia.csv": RESEARCH / "results" / "ob240669" / "gaia_dr3_5arcmin.csv",
        "xp_sampled.npz": RESEARCH / "data" / "ob240669" / "xp_sampled.npz",
    }
    for name, src in seeds.items():
        if not sl.ws.has(name):
            if not src.exists():
                sys.exit(f"seed missing: {src}")
            shutil.copy(src, sl.ws.path(name))
            print(f"seeded {name} <- {src}")

    if a.max_stars:
        import dustline.fit as fit_mod
        orig = fit_mod.fit_stars

        def limited(*args, **kw):
            kw["max_stars"] = a.max_stars
            return orig(*args, **kw)
        fit_mod.fit_stars = limited
        print(f"limiting to {a.max_stars} stars")

    res = sl.run(force=a.force)
    print()
    print(res)
    print("R_V:", res.rv)
    if res.clump_anchor:
        c = res.clump_anchor
        print(f"clump anchor: A_V {c['AV_column']:.2f} +/- {c['AV_column_err']:.2f} "
              f"at {c['D_RC']:.2f} kpc (E(J-Ks) {c['E_JK']:.2f}, n {c['n_window']})")
    tab = res.extinction("I")
    print("\nA_I(D):")
    for _, r in tab.iterrows():
        if abs(r.D_kpc * 2 - round(r.D_kpc * 2)) < 1e-6:
            b = " (bridged)" if r.bridged else ""
            print(f"  {r.D_kpc:5.1f} kpc  A_I = {r.A_med:5.2f} [{r.A_16:5.2f}, {r.A_84:5.2f}]"
                  f"  N {int(r.n_stars)}{b}")
    res.plots(band="I", directory="ob240669_plots")


if __name__ == "__main__":
    main()

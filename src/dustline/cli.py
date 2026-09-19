"""Command-line interface.

    dustline run RA DEC [--filter I] [--radius 5] [--phot my.csv --phot-bands g=DECam_g,Ks=VISTA_Ks]
                        [-o out.csv] [--force]
    dustline fetch-assets            # pre-download the model assets
    dustline clear-cache [--all]     # remove sightline caches (and assets with --all)
"""

from __future__ import annotations

import argparse
import shutil
import sys


def _cmd_run(a) -> int:
    from .api import Sightline
    from .userphot import UserPhotometry

    phot = None
    if a.phot:
        if not a.phot_bands:
            print("--phot needs --phot-bands key=BandName[,key=BandName...]", file=sys.stderr)
            return 2
        bands = dict(kv.split("=") for kv in a.phot_bands.split(","))
        phot = UserPhotometry(a.phot, bands)
    sl = Sightline(a.ra, a.dec, radius_arcmin=a.radius, photometry=phot,
                   prefer_deep=not a.simple_surveys)
    res = sl.run(force=a.force)
    tab = res.extinction(a.filter)
    rv = res.rv
    print(f"\nR_V = {rv['rv']:.2f} (MAD {rv['rv_mad']:.2f}, N {rv['n_stars']})")
    print(f"A_{a.filter}(D)  [A_{a.filter}/A_V = {tab.ratio_av.iloc[0]:.3f}]:")
    for _, r in tab.iterrows():
        if abs(r.D_kpc * 2 - round(r.D_kpc * 2)) < 1e-6:   # print every 0.5 kpc
            print(f"  {r.D_kpc:5.1f} kpc  A = {r.A_med:5.2f}  [{r.A_16:5.2f}, {r.A_84:5.2f}]"
                  f"  N = {int(r.n_stars)}")
    if a.out:
        res.save(a.out, band=a.filter)
        print(f"wrote {a.out}")
    if a.plots:
        res.plots(band=a.filter, directory=a.plots if a.plots != "." else None)
    return 0


def _cmd_fetch_assets(a) -> int:
    from . import assets

    for name in assets.ASSETS:
        assets.fetch(name, force=a.force)
    return 0


def _cmd_clear_cache(a) -> int:
    from . import assets

    root = assets.cache_dir()
    target = root if a.all else root / "sightlines"
    if target.exists():
        shutil.rmtree(target)
        print(f"removed {target}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="dustline",
                                description="Measured extinction law and 3D dust run per sightline")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run the pipeline for a sightline")
    r.add_argument("ra", type=float, help="RA (deg, ICRS)")
    r.add_argument("dec", type=float, help="Dec (deg, ICRS)")
    r.add_argument("--filter", default="I", help="output filter (default Cousins I)")
    r.add_argument("--radius", type=float, default=5.0, help="field radius (arcmin)")
    r.add_argument("--phot", default=None, help="user photometry csv (ra, dec, mag_<k>, magerr_<k>)")
    r.add_argument("--phot-bands", default=None, help="key=BandName[,key=BandName...] for --phot")
    r.add_argument("--simple-surveys", action="store_true",
                   help="use PS1+2MASS even where DECaPS/VVV are available")
    r.add_argument("-o", "--out", default=None, help="write the A(D) table to this csv")
    r.add_argument("--plots", nargs="?", const=".", default=None,
                   help="write the R_V and extinction-run figures (optionally to this directory)")
    r.add_argument("--force", action="store_true", help="ignore cached results")
    r.set_defaults(func=_cmd_run)

    f = sub.add_parser("fetch-assets", help="pre-download the model assets")
    f.add_argument("--force", action="store_true")
    f.set_defaults(func=_cmd_fetch_assets)

    c = sub.add_parser("clear-cache", help="remove cached sightline products")
    c.add_argument("--all", action="store_true", help="also remove the model assets")
    c.set_defaults(func=_cmd_clear_cache)

    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

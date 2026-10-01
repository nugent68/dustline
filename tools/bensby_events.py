"""Truth table for the microlensed-dwarf benchmark: the 91 spectroscopically analysed microlensed
bulge dwarfs/subgiants of Bensby et al. (2017, A&A 605, A89; VizieR J/A+A/605/A89).

For each source the table gives (V-I)_0 and M_V from SPECTROSCOPY (T_eff, log g, [Fe/H] + isochrones;
independent of dust) and (V-I)_0 and M_I from the MICROLENSING technique (Yoo et al. 2004: the source
dereddened relative to the red clump of its field, and placed at the clump distance). So

  dE(V-I) = (V-I)_0,mu - (V-I)_0,spec = E(V-I)_source - E(V-I)_clump
  dM_I    = M_I,mu - M_I,spec         = (mu_source - mu_clump) + (A_I,source - A_I,clump)

With the clump's A_I, E(V-I) and distance modulus of the field (Nataf et al. 2013, OGLE-III map, via
the OGLE extinction calculator, nearest grid point):

  E(V-I)_src = E(V-I)_RC + dE,   A_I,src = A_I,RC + (A_I/E(V-I))_RC dE,   mu_src = mu_RC + dM_I - (A_I/E(V-I))_RC dE

Errors: spectroscopic (V-I)_0 half-range (+) 0.05 mag for the microlensing colour (+) the map's
sigma E(V-I); mu from the spectroscopic M_V half-range. Caveat: the microlensing values used each
paper's own clump calibration, not necessarily Nataf's map, so dE/dM_I are relative to a clump
that may differ from the map's by a few hundredths.

Also: the dustline survey plan at each position (footprint.plan), the number of Gaia DR3 XP spectra
within 5' (pipeline cost ~1.8 s each), and (--map-only, any Python with dustmaps) the DECaPS 3D map
(Zucker et al. 2025) E(B-V) at the source and clump distances.

  python tools/bensby_events.py -o benchmarks/bensby2017/events.csv
  python tools/bensby_events.py --map-only --map-dir DUSTMAPS_DIR -o benchmarks/bensby2017/events.csv
"""
import argparse
import re
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

VIZ = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/605/A89/"
OGLE_EXT = "https://ogle.astrouw.edu.pl/cgi-ogle/getext.py"
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
UA = {"User-Agent": "dustline benchmark (bensby_events.py)"}


def get(url, data=None, timeout=120):
    return urllib.request.urlopen(urllib.request.Request(url, data=data, headers=UA), timeout=timeout).read().decode()


def bensby_table() -> pd.DataFrame:
    from astropy.io import ascii
    readme, dat = get(VIZ + "ReadMe"), get(VIZ + "tablea2.dat")
    import os
    import tempfile
    with tempfile.TemporaryDirectory() as tdir:
        open(os.path.join(tdir, "ReadMe"), "w").write(readme)
        open(os.path.join(tdir, "tablea2.dat"), "w").write(dat)
        t = ascii.read(os.path.join(tdir, "tablea2.dat"), readme=os.path.join(tdir, "ReadMe"), format="cds")
    f = lambda c: np.ma.filled(t[c].astype(float), np.nan)   # noqa: E731
    sgn = np.where(np.array(t["DE-"]) == "-", -1, 1)
    return pd.DataFrame(dict(
        name=[str(n).strip() for n in t["Name"]],
        ra=15 * (f("RAh") + f("RAm") / 60 + f("RAs") / 3600), dec=sgn * (f("DEd") + f("DEm") / 60 + f("DEs") / 3600),
        l_table=f("GLON"), b_table=f("GLAT"), teff=f("Teff"), e_teff=f("e_Teff"), logg=f("logg"), feh=f("[Fe/H]"),
        age=f("Age"), mass=f("M"), vi0_spec=f("(V-I)0"), vi0_spec_lo=f("b_(V-I)0"), vi0_spec_hi=f("B_(V-I)0"),
        MV_spec=f("VMAG"), MV_spec_lo=f("b_VMAG"), MV_spec_hi=f("B_VMAG"), MI_mu=f("IMAGmu"), vi0_mu=f("(V-I)0mu"),
        teff_mu=f("Teffmu"), tmax_hjd=f("Tmax"), amax=f("Amax")))


def nataf(ra, dec, pause=1.0):
    """OGLE-III clump extinction (Nataf+2013) at the nearest grid point: good (QF=0) points, else all."""
    cols = ["AI_RC", "EVI_RC", "sEVI_RC", "R_JKVI", "mu_RC", "smu_RC", "QF", "dist_ng_arcsec"]
    pat = re.compile(r"dist_NG\s*\n\s*" + r"\s+".join([r"([-\d.]+)"] * 8) + r"\s+(\d+)\s+NG\s+([-\d.]+)")
    for ds in ("good", "all"):
        q = urllib.parse.urlencode(dict(coord="dd", radec=f"{ra:.5f} {dec:.5f}", mode="ng", dataset=ds))
        m = pat.search(re.sub(r"<[^>]*>", " ", get(OGLE_EXT + "?" + q, timeout=60)))
        time.sleep(pause)
        if m:
            v = [float(x) for x in m.groups()]
            if v[2] > -9:
                return dict(zip(cols, v[2:]), nataf_set=ds)
    return dict({c: np.nan for c in cols}, nataf_set="none")


def n_xp(ra, dec, radius_arcmin=5.0):
    adql = (f"SELECT COUNT(*) AS nxp FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
            f"CIRCLE('ICRS', {ra:.5f}, {dec:.5f}, {radius_arcmin / 60:.6f})) AND has_xp_continuous = 'true'")
    q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql)).encode()
    return int(get(GAIA_TAP, data=q).splitlines()[1])


def build(out):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from dustline.catalogs import footprint

    d = bensby_table()
    g = SkyCoord(d.ra.values * u.deg, d.dec.values * u.deg).galactic
    d["l"] = (g.l.deg + 180) % 360 - 180
    d["b"] = g.b.deg
    d["lb_table_mismatch_deg"] = np.hypot(d.l - d.l_table, d.b - d.b_table).round(3)   # the table's l,b has typos
    d["MI_spec"] = d.MV_spec - d.vi0_spec
    d["dE_VI"] = d.vi0_mu - d.vi0_spec
    d["dM_I"] = d.MI_mu - d.MI_spec
    plans = [footprint.plan(r, de) for r, de in zip(d.ra, d.dec)]
    d["optical"] = [p.optical for p in plans]; d["nir"] = [p.nir for p in plans]
    d["mode"] = [p.mode for p in plans]; d["bulge_window"] = [bool(p.in_bulge_window) for p in plans]
    print(f"Bensby+2017: {len(d)} sources; querying the OGLE-III extinction calculator ...", flush=True)
    d = pd.concat([d, pd.DataFrame([nataf(r, de) for r, de in zip(d.ra, d.dec)])], axis=1)
    R = d.AI_RC / d.EVI_RC
    d["RI_RC"] = R
    d["EVI_src"] = d.EVI_RC + d.dE_VI
    d["AI_src"] = d.AI_RC + R * d.dE_VI
    d["mu_src"] = d.mu_RC + d.dM_I - R * d.dE_VI
    d["D_src_kpc"] = 10 ** (d.mu_src / 5 - 2)
    d["D_RC_kpc"] = 10 ** (d.mu_RC / 5 - 2)
    d["fAI_src_over_RC"] = d.AI_src / d.AI_RC
    e_vi = np.sqrt(((d.vi0_spec_hi - d.vi0_spec_lo) / 2) ** 2 + 0.05 ** 2)
    d["e_EVI_src"] = np.sqrt(e_vi ** 2 + d.sEVI_RC.fillna(0) ** 2)
    d["e_AI_src"] = R * d.e_EVI_src
    d["e_mu_src"] = np.sqrt(((d.MV_spec_hi - d.MV_spec_lo) / 2) ** 2 + (R * e_vi) ** 2)
    print("counting Gaia XP spectra within 5' ...", flush=True)
    counts = []
    for r, de in zip(d.ra, d.dec):
        counts.append(n_xp(r, de)); time.sleep(0.3)
    d["n_xp_5arcmin"] = counts
    d.round(5).to_csv(out, index=False)
    ok = d.AI_src.notna()
    print(f"wrote {out}: truth for {ok.sum()} sources; A_I,src/A_I,RC median {d.fAI_src_over_RC[ok].median():.3f}; "
          f"D_src median {d.D_src_kpc[ok].median():.2f} kpc; XP spectra per field median {np.median(counts):.0f}")


def add_map(out, map_dir):
    import os

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from dustmaps.config import config
    config["data_dir"] = os.path.expanduser(map_dir)
    from dustmaps.decaps import DECaPSQueryLite
    q = DECaPSQueryLite(mean_only=True)
    d = pd.read_csv(out)
    for col, D in (("ebv_decaps_Dsrc", d.D_src_kpc), ("ebv_decaps_DRC", d.D_RC_kpc)):
        ok = D.notna()
        c = SkyCoord(d.ra.values * u.deg, d.dec.values * u.deg, distance=D.fillna(8.0).values * u.kpc)
        d[col] = np.where(ok, np.asarray(q.query(c, mode="mean"), float), np.nan)
    for D in (2.0, 4.0, 8.0):
        c = SkyCoord(d.ra.values * u.deg, d.dec.values * u.deg, distance=np.full(len(d), D) * u.kpc)
        d[f"ebv_decaps_{D:g}kpc"] = np.asarray(q.query(c, mode="mean"), float)
    # first batch for full runs: 3 per quintile of A_I,src among the DECaPS-photometry, map-covered events
    # with truth, best-measured source distance first, then the cheapest field
    elig = d.AI_src.notna() & d.D_src_kpc.notna() & (d.optical == "decaps") & np.isfinite(d.ebv_decaps_Dsrc)
    e = d[elig].assign(AI_bin=pd.qcut(d.AI_src[elig], 5, labels=False))
    pick = e.sort_values(["AI_bin", "e_mu_src", "n_xp_5arcmin"]).groupby("AI_bin").head(3)
    d["stage_b"] = d.name.isin(pick.name)
    d.round(5).to_csv(out, index=False)
    print(f"stage B batch: {len(pick)} fields, {pick.n_xp_5arcmin.sum():.0f} XP spectra")
    print(f"DECaPS 3D map added: {np.isfinite(d.ebv_decaps_DRC).sum()} of {len(d)} positions covered")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--map-only", action="store_true", help="only add the DECaPS 3D map columns (needs dustmaps)")
    ap.add_argument("--map-dir", default="~/claude/declens/data/dustmaps")
    a = ap.parse_args()
    if a.map_only:
        add_map(a.out, a.map_dir)
    else:
        build(a.out)


if __name__ == "__main__":
    main()

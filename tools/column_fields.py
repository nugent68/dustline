"""Column-mode zero point against the DESI stellar-reddening map (Zhou+24) over many high-latitude
fields: is the COSMOS deficit (dustline 0.02-0.03 below DESI/SFD) an offset or a scale error?

select   candidate field centres (|b| > 30, Dec > -25: column mode with DESI priors) on an
         nside-64 grid; per candidate the DESI map averaged over the field (nside 512, stars
         weighted), its smoothness, SFD, and the DESI DR1 MWS density (rows of the mirror's hp32
         file); one field per log bin of DESI E(B-V) between 0.01 and 0.1 (most MWS stars, >= 5 deg
         apart), plus COSMOS -> fields.csv
compare  after `field_batch.py ... fit`: dustline column vs DESI per field, the straight-line fit
         dustline = a + b x DESI (bootstrap over fields), also vs SFD x 0.86

  python tools/column_fields.py select  --out DIR [--n 19]
  python tools/column_fields.py compare --out DIR
"""
import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
DESI_MAPS = os.environ.get("DESI_DUST_MAPS", "/global/cfs/cdirs/desi/public/papers/mws/desi_dust/y2/v1/maps")
MWS_HP32 = os.environ.get("MWS_HP32", "/global/cfs/projectdirs/newera/surveys/desi_dr1/mws/hp32")
RADIUS_DEG = 0.5
COSMOS = dict(name="COSMOS", ra=150.12, dec=2.21)


def a_v_per_e():
    """E(g-r)/A_V and E(r-z)/A_V under G23 R_V 3.1 for a 5800 K dwarf, DECam (as for COSMOS)."""
    from dustline import ensemble
    r = {b: ensemble.band_ratio_at_rv(b, 3.1, av=0.05, teff=5800.0, logg=4.3) for b in ("DECam_g", "DECam_r", "DECam_z")}
    return r["DECam_g"] - r["DECam_r"], r["DECam_r"] - r["DECam_z"]


def load_map(tag, nside=512):
    from astropy.io import fits
    t = fits.getdata(os.path.join(DESI_MAPS, f"desi_dust_{tag}_{nside}.fits"))
    full = {}
    hp = np.asarray(t["HPXPIXEL"]).astype(np.int64)
    col = "EGR" if tag == "gr" else "ERZ"
    for k, src in ((col, col), ("N", "N_STAR"), ("SFD", "EBV_SFD")):
        a = np.full(12 * nside ** 2, np.nan)
        a[hp] = np.asarray(t[src]).astype(float)
        full[k] = a
    return full


def field_stats(maps, ra, dec, k_gr, k_rz):
    import healpy as hp
    v = hp.ang2vec(ra, dec, lonlat=True)
    pix = hp.query_disc(512, v, np.deg2rad(RADIUS_DEG))
    out = {}
    for tag, k in (("gr", k_gr), ("rz", k_rz)):
        m = maps[tag]
        e, n = m["EGR" if tag == "gr" else "ERZ"][pix], m["N"][pix]
        ok = np.isfinite(e) & (n > 0)
        out[f"n_{tag}"] = int(n[ok].sum())
        out[f"AV_desi_{tag}"] = float(np.sum(e[ok] * n[ok]) / n[ok].sum() / k) if ok.any() else np.nan
        out[f"AV_desi_{tag}_err"] = float(np.std(e[ok]) / np.sqrt(ok.sum()) / k) if ok.sum() > 1 else np.nan
        out[f"AV_desi_{tag}_pixsd"] = float(np.std(e[ok]) / k) if ok.any() else np.nan
        if tag == "gr":
            s = m["SFD"][pix]
            out["AV_sfd086"] = float(2.742 * np.nanmean(s))
    return out


def mws_rows(ra, dec):
    import healpy as hp
    from astropy.io import fits
    p = hp.ang2pix(32, ra, dec, lonlat=True, nest=True)
    f = os.path.join(MWS_HP32, f"{p:05d}.fits")
    if not os.path.exists(f):
        return 0
    with fits.open(f, memmap=True) as h:
        return int(h[1].header.get("NAXIS2", 0))


def cmd_select(a):
    import healpy as hp
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    k_gr, k_rz = a_v_per_e()
    maps = {t: load_map(t) for t in ("gr", "rz")}
    ra, dec = hp.pix2ang(64, np.arange(12 * 64 ** 2), lonlat=True)
    gb = SkyCoord(ra * u.deg, dec * u.deg).galactic.b.deg
    cand = np.where((np.abs(gb) > 30.0) & (dec > -25.0))[0]
    # quick pre-screen on the nside-512 map at the centre (has DESI stars)
    cpix = hp.ang2pix(512, ra[cand], dec[cand], lonlat=True)
    cand = cand[np.isfinite(maps["gr"]["EGR"][cpix]) & (maps["gr"]["N"][cpix] > 0)]
    rows = []
    for i in cand:
        st = field_stats(maps, ra[i], dec[i], k_gr, k_rz)
        if st["n_gr"] < 150:
            continue
        rows.append(dict(name=f"F{i:05d}", ra=round(float(ra[i]), 4), dec=round(float(dec[i]), 4), b=round(float(gb[i]), 2), **st))
    c = pd.DataFrame(rows)
    c["mws_rows_hp32"] = [mws_rows(r, d) for r, d in zip(c.ra, c.dec)]
    c["EBV_desi_gr"] = c.AV_desi_gr / 3.1
    smooth = c.AV_desi_gr_pixsd < np.maximum(0.6 * c.AV_desi_gr, 0.03)
    c = c[smooth & (c.mws_rows_hp32 >= 1000) & c.AV_desi_gr.notna()].reset_index(drop=True)
    print(f"{len(c)} candidate fields after cuts (>= 150 DESI-map stars in 30', smooth, >= 1000 DR1 MWS rows in the hp32 pixel)")
    edges = np.logspace(np.log10(0.01), np.log10(0.10), a.n + 1)
    picked = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = c[(c.EBV_desi_gr >= lo) & (c.EBV_desi_gr < hi)].sort_values("mws_rows_hp32", ascending=False)
        for r in s.itertuples():
            if all(SkyCoord(r.ra * u.deg, r.dec * u.deg).separation(SkyCoord(p.ra * u.deg, p.dec * u.deg)).deg > 5 for p in picked):
                picked.append(r)
                break
    out = pd.DataFrame([p._asdict() for p in picked]).drop(columns=["Index"])
    cs = field_stats(maps, COSMOS["ra"], COSMOS["dec"], k_gr, k_rz)
    out = pd.concat([out, pd.DataFrame([dict(COSMOS, b=42.1, **cs, mws_rows_hp32=mws_rows(COSMOS["ra"], COSMOS["dec"]),
                                             EBV_desi_gr=cs["AV_desi_gr"] / 3.1)])], ignore_index=True)
    os.makedirs(a.out, exist_ok=True)
    out.round(5).to_csv(os.path.join(a.out, "fields.csv"), index=False)
    print(out[["name", "ra", "dec", "b", "EBV_desi_gr", "AV_desi_gr", "AV_desi_rz", "AV_sfd086", "n_gr", "mws_rows_hp32"]].round(4).to_string(index=False))


def cmd_compare(a):
    f = pd.read_csv(os.path.join(a.out, "fields.csv"))
    rows = []
    for r in f.itertuples():
        p = os.path.join(a.out, f"{r.name}.fit.json")
        if not os.path.exists(p):
            continue
        j = json.load(open(p))
        c = (j.get("column") or {}).get("clean") or {}
        if not j.get("ok") or not c:
            continue
        rows.append(dict(name=r.name, AV_dl=c["av"], AV_dl_err=c["av_err"], n_dl=c["n"], AV_dl_all=j["column"]["av"],
                         offsets=j.get("phot_offsets")))
    s = f.merge(pd.DataFrame(rows), on="name")
    s.drop(columns=["offsets"]).round(5).to_csv(os.path.join(a.out, "column_vs_desi.csv"), index=False)
    L = [f"{len(s)} fields with a dustline column", s[["name", "b", "AV_desi_gr", "AV_desi_rz", "AV_sfd086", "AV_dl", "AV_dl_err", "n_dl"]].round(4).to_string(index=False)]
    rng = np.random.default_rng(1)
    for ref in ("AV_desi_gr", "AV_desi_rz", "AV_sfd086"):
        x, y, e = s[ref].values, s.AV_dl.values, s.AV_dl_err.values
        ok = np.isfinite(x) & np.isfinite(y)
        x, y, e = x[ok], y[ok], e[ok]
        w = 1 / np.maximum(e, 0.003) ** 2
        b, a0 = np.polyfit(x, y, 1, w=np.sqrt(w))
        bs = []
        for _ in range(2000):
            i = rng.integers(0, len(x), len(x))
            if np.ptp(x[i]) > 0:
                bs.append(np.polyfit(x[i], y[i], 1, w=np.sqrt(w[i])))
        bs = np.array(bs)
        d = y - x
        L.append(f"dustline = a + b x {ref}: a {a0:+.4f} +/- {bs[:, 1].std():.4f}, b {b:.3f} +/- {bs[:, 0].std():.3f}; "
                 f"median difference {np.median(d):+.4f} (robust sd {1.4826 * np.median(np.abs(d - np.median(d))):.4f}), "
                 f"median ratio {np.median(y / x):.3f}  (N {ok.sum()})")
    txt = "\n".join(L)
    open(os.path.join(a.out, "column_vs_desi.txt"), "w").write(txt + "\n")
    print(txt)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select"); s.add_argument("--out", required=True); s.add_argument("--n", type=int, default=19)
    c = sub.add_parser("compare"); c.add_argument("--out", required=True)
    a = ap.parse_args()
    {"select": cmd_select, "compare": cmd_compare}[a.cmd](a)


if __name__ == "__main__":
    main()

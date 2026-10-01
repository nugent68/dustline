"""How much dust do the template-correction calibrators really sit behind?

The dwarf template corrections (tools/build_template_corrections.py, dustline.calib) are measured
on DESI DR1 MWS dwarfs within 250 pc at |b| > 40 deg, dereddened with the Edenhofer+2023 3D map
(calib.map_extinction, A_V = 2.8 E) and then fitted at A_V = 0. Dust the map does not count in
front of them ends up in the corrected templates and is subtracted from every column measured with
them (column mode: COSMOS reads A_V 0.019-0.027 against DESI/SFD 0.046-0.059).

Per calibrator (re-selected exactly as calib.select_calibrators does):
  A_front_E   Edenhofer A_V to the star (the dereddening that was applied)
  A_tot_E     Edenhofer A_V to 1.25 kpc (its total column at |b| > 40)
  A_tot_DESI  Zhou+24 DESI stellar-reddening total column, E(g-r) and E(r-z) -> A_V (G23 3.1, DECam)
  A_tot_SFD   SFD x 0.86 (A_V = 2.742 E_SFD, Schlafly & Finkbeiner 2011)
  f_front     fraction of the column in front of the star: Edenhofer's own profile, and an
              exponential layer (scale height 125 pc) as a cross-check
  missing     f_front x A_tot(DESI or SFD) - A_front_E

  python tools/calibrator_dust.py [--out DIR]          (needs dustmaps; see the NERSC surveys README)
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
EDEN_DIR = os.environ.get("EDENHOFER_DIR", "/global/cfs/projectdirs/newera/surveys")
DESI_MAPS = os.environ.get("DESI_DUST_MAPS", "/global/cfs/cdirs/desi/public/papers/mws/desi_dust/y2/v1/maps")
H_DUST_PC = 125.0
COSMOS = (150.12, 2.21, 0.5)        # deg, deg, radius deg


def edenhofer():
    from dustmaps.config import config
    config["data_dir"] = EDEN_DIR
    from dustmaps.edenhofer2023 import Edenhofer2023Query
    return Edenhofer2023Query(integrated=True)


def eden_av(q, ra, dec, d_pc):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg,
                 distance=np.clip(np.asarray(d_pc, float), 70.0, 1240.0) * u.pc, frame="icrs")
    return 2.8 * np.nan_to_num(np.asarray(q(c), float), nan=0.0)


def desi_lookup(ra, dec, nside=512):
    import healpy as hp
    from astropy.io import fits
    out = {}
    pix = hp.ang2pix(nside, np.asarray(ra), np.asarray(dec), lonlat=True)    # RING, as the maps
    for tag, col in (("gr", "EGR"), ("rz", "ERZ")):
        t = fits.getdata(os.path.join(DESI_MAPS, f"desi_dust_{tag}_{nside}.fits"))
        idx = pd.Series(np.arange(len(t)), index=t["HPXPIXEL"])
        j = idx.reindex(pix).values
        ok = np.isfinite(j)
        jj = np.where(ok, j, 0).astype(int)
        good = ok & (t["N_STAR"][jj] > 0)
        out[f"E_{tag}"] = np.where(good, t[col][jj], np.nan)
        out[f"n_{tag}"] = np.where(ok, t["N_STAR"][jj], 0)
        out[f"sfd_{tag}"] = np.where(ok, t["EBV_SFD"][jj], np.nan)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    from dustline import calib, ensemble
    sed = dict(teff=5800.0, logg=4.3)
    r = {b: ensemble.band_ratio_at_rv(b, 3.1, av=0.05, **sed) for b in ("DECam_g", "DECam_r", "DECam_z")}
    k_gr, k_rz = r["DECam_g"] - r["DECam_r"], r["DECam_r"] - r["DECam_z"]
    cal = calib.select_calibrators()
    cal = cal[np.isfinite(cal.parallax) & (cal.parallax > 0)].reset_index(drop=True)
    d = 1000.0 / cal.parallax.values
    print(f"calibrators: {len(cal)} DESI dwarfs, D {np.median(d):.0f} pc [{np.percentile(d, 16):.0f}, {np.percentile(d, 84):.0f}], "
          f"|b| {np.median(np.abs(cal.b)):.0f} deg", flush=True)
    q = edenhofer()
    cal["A_front_E"] = eden_av(q, cal.ra, cal.dec, d)
    cal["A_tot_E"] = eden_av(q, cal.ra, cal.dec, np.full(len(cal), 1240.0))
    m = desi_lookup(cal.ra.values, cal.dec.values)
    cal["A_tot_DESI_gr"] = m["E_gr"] / k_gr
    cal["A_tot_DESI_rz"] = m["E_rz"] / k_rz
    cal["A_tot_SFD"] = 2.742 * m["sfd_gr"]
    z = d * np.abs(np.sin(np.deg2rad(cal.b.values)))
    cal["f_front_E"] = np.where(cal.A_tot_E > 0.003, cal.A_front_E / cal.A_tot_E, np.nan)
    cal["f_front_exp"] = 1.0 - np.exp(-z / H_DUST_PC)
    for ref in ("DESI_gr", "DESI_rz", "SFD"):
        for f in ("E", "exp"):
            cal[f"missing_{ref}_{f}"] = cal[f"f_front_{f}"] * cal[f"A_tot_{ref}"] - cal.A_front_E
    cal.round(5).to_csv(os.path.join(a.out, "calibrator_dust.csv"), index=False)

    med = lambda x: float(np.nanmedian(x))                                   # noqa: E731
    L = [f"{len(cal)} calibrators; with a DESI map value: {np.isfinite(cal.A_tot_DESI_gr).sum()}",
         f"  A_V in front, Edenhofer (applied):       {med(cal.A_front_E):.4f}",
         f"  total column, Edenhofer (to 1.25 kpc):    {med(cal.A_tot_E):.4f}",
         f"  total column, DESI g-r / r-z:             {med(cal.A_tot_DESI_gr):.4f} / {med(cal.A_tot_DESI_rz):.4f}",
         f"  total column, SFD x 0.86:                 {med(cal.A_tot_SFD):.4f}",
         f"  fraction in front: Edenhofer profile {med(cal.f_front_E):.2f}, exponential layer (h {H_DUST_PC:.0f} pc) {med(cal.f_front_exp):.2f}",
         "  missing dust in front of the calibrators (median A_V; f x total - Edenhofer front):"]
    for ref in ("DESI_gr", "DESI_rz", "SFD"):
        L.append(f"    vs {ref:8s}: Edenhofer profile {med(cal[f'missing_{ref}_E']):+.4f}, exponential layer {med(cal[f'missing_{ref}_exp']):+.4f}")
    bins = [(3500, 4500), (4500, 5500), (5500, 6500), (6500, 7500)]
    L.append("  by T_eff (DESI label), missing vs DESI g-r, Edenhofer profile / exponential layer:")
    for lo, hi in bins:
        s = cal[(cal.teff >= lo) & (cal.teff < hi)]
        if len(s):
            L.append(f"    {lo}-{hi} K: N {len(s):4d}  {med(s.missing_DESI_gr_E):+.4f} / {med(s.missing_DESI_gr_exp):+.4f}"
                     f"   (front E {med(s.A_front_E):.4f}, D {med(1000 / s.parallax):.0f} pc)")
    # the COSMOS field itself: the three total columns
    rng = np.random.default_rng(1)
    rr = COSMOS[2] * np.sqrt(rng.random(400)); th = 2 * np.pi * rng.random(400)
    ra = COSMOS[0] + rr * np.cos(th) / np.cos(np.deg2rad(COSMOS[1])); dec = COSMOS[1] + rr * np.sin(th)
    cm = desi_lookup(ra, dec)
    L.append(f"COSMOS (30'): total column Edenhofer {np.mean(eden_av(q, ra, dec, np.full(400, 1240.0))):.4f}, "
             f"DESI g-r {np.nanmean(cm['E_gr']) / k_gr:.4f}, r-z {np.nanmean(cm['E_rz']) / k_rz:.4f}, SFD x 0.86 {2.742 * np.nanmean(cm['sfd_gr']):.4f}; "
             "dustline column 0.019-0.027")
    txt = "\n".join(L)
    open(os.path.join(a.out, "calibrator_dust.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()

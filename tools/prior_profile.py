"""Convert a dustline run into the inputs of the declens source-distance prior
(`sed_distance_prior.py --profile <xp_dust_run.csv> --profile-unit ai --laws measured`,
which reads the research prototype's formats):

  xp_dust_run.csv   D_kpc, AV_target, AV_nb_min, AV_nb_med, AV_nb_max, N   -- in A_i (DECam i)
  xp_law.json       rv, rv_mad, n_stars, ratios_source_g23 {g r i z Y J H Ks Ks2 Kp: A_X/A_i}
                    for the microlensed source (a NewEra SED at --teff/--logg under G23 at the
                    measured R_V), plus the per-star field ratios and the clump anchor
  xp_stars_phot.csv the per-star table in the prototype's column names (teff, av, rv, A_i,
                    D_med/D_lo/D_hi, plx_inflate, ...) for the declens figure scripts
  clump_anchor.json the clump anchor in the prototype's keys (Ai_colour, D_RC_Ks, ...)

  python tools/prior_profile.py RA DEC --radius 5 [--plx-inflate 1.7] [--teff 6000 --logg 4.0] -o OUTDIR
"""

from __future__ import annotations

import argparse
import json
import os

from dustline import ensemble
from dustline.api import Sightline

# prototype key -> package filter name (DECam + VISTA for the 0095 sightline; Kp = Keck NIRC2)
SOURCE_BANDS = {"g": "DECam_g", "r": "DECam_r", "i": "DECam_i", "z": "DECam_z", "Y": "DECam_Y",
                "J": "VISTA_J", "H": "VISTA_H", "Ks": "VISTA_Ks", "Ks2": "2MASS_Ks", "Kp": "Keck_Kp"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ra", type=float); ap.add_argument("dec", type=float)
    ap.add_argument("--radius", type=float, default=5.0)
    ap.add_argument("--teff", type=float, default=6000.0, help="source SED for the band ratios")
    ap.add_argument("--logg", type=float, default=4.0)
    ap.add_argument("--band", default="DECam_i", help="profile band (the prior works in DECam i)")
    ap.add_argument("--plx-inflate", type=float, default=1.0, help="select the workspace run with this parallax-error inflation")
    ap.add_argument("-o", "--outdir", default=".")
    a = ap.parse_args()
    res = Sightline(a.ra, a.dec, radius_arcmin=a.radius, plx_inflate=a.plx_inflate).run()   # cached: instant after a run
    law = res.law
    ext = res.extinction(a.band)
    prof = ext.rename(columns={"A_med": "AV_target", "A_16": "AV_nb_min", "A_84": "AV_nb_max", "n_stars": "N"})
    prof["AV_nb_med"] = prof.AV_target
    prof = prof[["D_kpc", "AV_target", "AV_nb_min", "AV_nb_med", "AV_nb_max", "N"]]
    os.makedirs(a.outdir, exist_ok=True)
    prof.round(4).to_csv(os.path.join(a.outdir, "xp_dust_run.csv"), index=False)
    rv = law["rv"]
    src = {k: ensemble.band_ratio_at_rv(b, rv, teff=a.teff, logg=a.logg) for k, b in SOURCE_BANDS.items()}
    ratios = {k: v / src["i"] for k, v in src.items()}
    lo = {k: ensemble.band_ratio_at_rv(b, max(rv - law["rv_mad"], 2.3), teff=a.teff, logg=a.logg) for k, b in SOURCE_BANDS.items()}
    hi = {k: ensemble.band_ratio_at_rv(b, min(rv + law["rv_mad"], 5.55), teff=a.teff, logg=a.logg) for k, b in SOURCE_BANDS.items()}
    out = dict(rv=rv, rv_mad=law["rv_mad"], rv_raw=law.get("rv_raw"), rv_mean=law.get("rv_mean"), rv_mean_err=law.get("rv_mean_err"),
               n_stars=law["n_stars"], source_model=[a.teff, a.logg], profile_band=a.band, profile_ratio_av=float(ext.ratio_av.iloc[0]),
               ratios_source_g23=ratios,
               ratios_source_g23_lo={k: min(lo[k] / lo["i"], hi[k] / hi["i"]) for k in ratios},
               ratios_source_g23_hi={k: max(lo[k] / lo["i"], hi[k] / hi["i"]) for k in ratios},
               ratios_av_field=law["ratios_av"], clump_anchor=law.get("clump_anchor"), apogee_check=law.get("apogee_check"),
               phot_offsets=law.get("phot_offsets"), template_corrections=law.get("template_corrections"), model_cache=law.get("model_cache"),
               source="dustline package run; ratios for the source SED under G23 at the measured R_V, normalised to DECam i")
    json.dump(out, open(os.path.join(a.outdir, "xp_law.json"), "w"), indent=1)
    # per-star table and clump anchor in the prototype's conventions (declens dust_run_figure.py)
    from dustline import ensemble as _ens
    st = res.stars.copy()
    _ens.distance_posterior(st)
    ri = law["ratios_av"].get("DECam_i", ext.ratio_av.iloc[0])
    st["A_i"] = st["A_DECam_i"] if "A_DECam_i" in st else st.av * ri
    st["A_i_err"] = st["A_DECam_i_err"] if "A_DECam_i_err" in st else st.av_err * ri
    if "plx_inflate" not in st:
        st["plx_inflate"] = 1.0
    st.round(5).to_csv(os.path.join(a.outdir, "xp_stars_phot.csv"), index=False)
    c = law.get("clump_anchor")
    if c:
        json.dump(dict(Ai_colour=c["AV_column"] * ri, Ai_conv_err=0.0, Ai_diff_spread=c["AV_column_err"] * ri,
                       D_RC_Ks=c["D_RC"], D_RC_i=c["D_RC"], n_window=c["n_window"], n_est=1, E_JK=c["E_JK"], AV_column=c["AV_column"],
                       ratio_i=ri, source="dustline clump.find_clump (J-Ks), converted to DECam i"),
                  open(os.path.join(a.outdir, "clump_anchor.json"), "w"), indent=1)
    print(f"R_V {rv:.2f} +/- {law['rv_mad']:.2f} ({law['n_stars']} stars; raw {law.get('rv_raw')}); A_X/A_i for the {a.teff:.0f} K source: "
          + " ".join(f"{k} {v:.3f}" for k, v in ratios.items()))
    for d in (1, 2, 3, 4, 5, 6, 8):
        row = prof.iloc[(prof.D_kpc - d).abs().argmin()]
        print(f"  D {d} kpc: A_i {row.AV_target:.2f} [{row.AV_nb_min:.2f}, {row.AV_nb_max:.2f}] N {int(row.N)}")
    print("wrote", os.path.join(a.outdir, "xp_dust_run.csv"), "and xp_law.json")


if __name__ == "__main__":
    main()

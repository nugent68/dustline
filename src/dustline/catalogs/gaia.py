"""Gaia DR3: cone query (astrometry + photometry + flags + 2MASS join) and the
XP continuous-coefficient fetch (resumable, chunked POST) + gaiaxpy calibration
to sampled spectra.

All products land in the sightline Workspace:
  gaia.csv                the cone query result
  xp_continuous_raw.csv   raw XP coefficients (append-mode, resumable)
  xp_sampled.npz          calibrated spectra: source_id, wave_nm, flux, flux_err
"""

from __future__ import annotations

import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

from ..cache import Workspace

GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
GAIA_DATA = "https://gea.esac.esa.int/data-server/data"
UA = {"User-Agent": "dustline"}

_ADQL = """
SELECT g.source_id, g.ra, g.dec, g.l, g.b,
       g.parallax, g.parallax_error, g.parallax_over_error,
       g.pmra, g.pmra_error, g.pmdec, g.pmdec_error,
       g.ruwe, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude,
       g.astrometric_params_solved, g.visibility_periods_used,
       g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
       g.phot_g_mean_flux_over_error, g.phot_bp_mean_flux_over_error,
       g.phot_rp_mean_flux_over_error, g.phot_bp_rp_excess_factor,
       g.phot_bp_n_obs, g.phot_rp_n_obs, g.bp_rp,
       g.has_xp_continuous, g.has_xp_sampled, g.has_rvs,
       g.teff_gspphot, g.ag_gspphot, g.distance_gspphot,
       g.non_single_star, g.in_qso_candidates, g.in_galaxy_candidates,
       t.j_m AS mag_J, t.j_msigcom AS magerr_J,
       t.h_m AS mag_H, t.h_msigcom AS magerr_H,
       t.ks_m AS mag_Ks, t.ks_msigcom AS magerr_Ks,
       t.ph_qual AS tmass_ph_qual, xm.angular_distance AS tmass_sep
FROM gaiadr3.gaia_source AS g
LEFT JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xm ON xm.source_id = g.source_id
LEFT JOIN gaiadr1.tmass_original_valid AS t ON t.designation = xm.original_ext_source_id
WHERE 1=CONTAINS(POINT('ICRS', g.ra, g.dec),
                 CIRCLE('ICRS', {ra}, {dec}, {radius_deg}))
"""

# the TAP service lower-cases aliases; restore the band-case names the fits expect
_CASE_FIX = [("mag_j", "mag_J"), ("magerr_j", "magerr_J"), ("mag_h", "mag_H"),
             ("magerr_h", "magerr_H"), ("mag_ks", "mag_Ks"), ("magerr_ks", "magerr_Ks")]


def by_ids(source_ids, chunk: int = 100) -> pd.DataFrame:
    """The cone-query columns (incl. the 2MASS join) for an explicit list of Gaia DR3
    source_ids, fetched in chunks (for calibration samples that are not a cone)."""
    ids = [int(i) for i in source_ids]
    frames = []
    for k in range(0, len(ids), chunk):
        part = ids[k:k + chunk]
        adql = _ADQL[:_ADQL.index("WHERE 1=CONTAINS")] + \
            "WHERE g.source_id IN (" + ",".join(map(str, part)) + ")"
        q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql))
        req = urllib.request.Request(GAIA_TAP, data=q.encode(), headers=UA)   # POST: long id lists
        text = urllib.request.urlopen(req, timeout=900).read().decode()
        head, body = text.split("\n", 1)
        for lo, hi in _CASE_FIX:
            head = head.replace(lo, hi)
        import io
        frames.append(pd.read_csv(io.StringIO(head + "\n" + body), low_memory=False))
        print(f"  gaia by_ids: {min(k + chunk, len(ids))}/{len(ids)}", flush=True)
    return pd.concat(frames, ignore_index=True)


def cone(ws: Workspace, force: bool = False) -> pd.DataFrame:
    """Gaia DR3 sources within the workspace radius; cached as gaia.csv."""
    out = ws.path("gaia.csv")
    if out.exists() and not force:
        return pd.read_csv(out, low_memory=False)
    adql = _ADQL.format(ra=ws.ra, dec=ws.dec, radius_deg=ws.radius_arcmin / 60.0)
    q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql))
    req = urllib.request.Request(GAIA_TAP + "?" + q, headers=UA)
    text = urllib.request.urlopen(req, timeout=900).read().decode()
    nrows = text.count("\n") - 1
    if nrows < 50:
        raise RuntimeError(f"Gaia cone query returned only {nrows} rows:\n{text[:500]}")
    head, body = text.split("\n", 1)
    for lo, hi in _CASE_FIX:
        head = head.replace(lo, hi)
    out.write_text(head + "\n" + body)
    g = pd.read_csv(out, low_memory=False)
    print(f"gaia cone: {len(g)} sources, {int((g.has_xp_continuous == True).sum())} with XP")  # noqa: E712
    return g


def _fetch_chunk(ids: list[int], retries: int = 4) -> str:
    data = urllib.parse.urlencode(dict(
        RETRIEVAL_TYPE="XP_CONTINUOUS", DATA_STRUCTURE="RAW", FORMAT="csv",
        ID=",".join(f"Gaia DR3 {i}" for i in ids))).encode()
    for k in range(retries):
        try:  # POST: a GET URL with >~100 ids is rejected (HTTP 414)
            req = urllib.request.Request(GAIA_DATA, data=data, headers=UA)
            return urllib.request.urlopen(req, timeout=900).read().decode()
        except Exception as e:  # noqa: BLE001
            print(f"  retry {k + 1}: {e}", flush=True)
            time.sleep(10 * (k + 1))
    raise RuntimeError("Gaia DataLink XP fetch failed after retries")


def fetch_xp(ws: Workspace, gaia: pd.DataFrame, chunk: int = 200) -> None:
    """XP continuous coefficients for every has_xp_continuous source; resumable."""
    out = ws.path("xp_continuous_raw.csv")
    ids = gaia.loc[gaia.has_xp_continuous == True, "source_id"].astype(int).tolist()  # noqa: E712
    done: set[int] = set()
    if out.exists():
        done = set(pd.read_csv(out, usecols=["source_id"]).source_id.astype(int))
    todo = [i for i in ids if i not in done]
    print(f"XP: {len(ids)} sources, {len(done)} fetched, {len(todo)} to go "
          f"(~{len(todo) * 1.8 / 60:.0f} min)", flush=True)
    header_written = out.exists()
    for k in range(0, len(todo), chunk):
        txt = _fetch_chunk(todo[k:k + chunk])
        lines = txt.splitlines()
        if not lines or not lines[0].startswith("source_id"):
            print(f"  chunk {k}: unexpected response ({txt[:200]!r})", flush=True)
            continue
        with open(out, "a") as fh:
            if not header_written:
                fh.write(lines[0] + "\n")
                header_written = True
            fh.write("\n".join(lines[1:]) + "\n")
        print(f"  {min(k + chunk, len(todo))}/{len(todo)}", flush=True)


def calibrate_xp(ws: Workspace, chunk: int = 500, truncation: bool = False,
                 force: bool = False):
    """gaiaxpy calibrate() of the raw coefficients -> xp_sampled.npz (336-1020 nm, 2 nm)."""
    out = ws.path("xp_sampled.npz")
    if out.exists() and not force:
        return np.load(out)
    from gaiaxpy import calibrate
    df = pd.read_csv(ws.path("xp_continuous_raw.csv")).drop_duplicates("source_id")
    print(f"calibrating {len(df)} XP spectra", flush=True)
    ids, flux, err = [], [], []
    sampling = None
    for k in range(0, len(df), chunk):
        cal, sampling = calibrate(df.iloc[k:k + chunk], truncation=truncation,
                                  save_file=False)
        ids += cal.source_id.astype(np.int64).tolist()
        flux += [np.asarray(f, dtype=np.float32) for f in cal.flux]
        err += [np.asarray(f, dtype=np.float32) for f in cal.flux_error]
        print(f"  {min(k + chunk, len(df))}/{len(df)}", flush=True)
    np.savez_compressed(out, source_id=np.array(ids, dtype=np.int64),
                        wave_nm=np.asarray(sampling), flux=np.array(flux),
                        flux_err=np.array(err))
    return np.load(out)

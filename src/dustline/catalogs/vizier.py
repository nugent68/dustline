"""CDS VizieR TAP sync helper (stdlib only), the same box-then-circle pattern as
catalogs.datalab.box_query."""

from __future__ import annotations

import io
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
UA = {"User-Agent": "dustline"}


def box_query(table: str, columns: list[str], ra: float, dec: float, radius_deg: float,
              ra_col: str = "RAJ2000", dec_col: str = "DEJ2000", extra: str = "") -> pd.DataFrame:
    cosd = max(np.cos(np.radians(dec)), 1e-3)
    dra = radius_deg / cosd
    adql = (f'SELECT {", ".join(columns)} FROM "{table}" '
            f"WHERE {ra_col} BETWEEN {ra - dra} AND {ra + dra} "
            f"AND {dec_col} BETWEEN {dec - radius_deg} AND {dec + radius_deg} {extra}")
    q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql))
    req = urllib.request.Request(TAP + "?" + q, headers=UA)
    txt = urllib.request.urlopen(req, timeout=900).read().decode()
    df = pd.read_csv(io.StringIO(txt))
    if not len(df):
        return df
    sep = np.hypot((df[ra_col] - ra) * cosd, df[dec_col] - dec)
    return df[sep <= radius_deg].reset_index(drop=True)


def positions_query(table: str, columns: list[str], ra, dec, radius_arcsec: float,
                    ra_col: str = "RAJ2000", dec_col: str = "DEJ2000",
                    chunk: int = 80, extra: str = "") -> pd.DataFrame:
    """All rows within radius_arcsec of any of the positions: OR'd boxes, chunked
    (for calibration samples scattered over the sky)."""
    ra, dec = np.asarray(ra, float), np.asarray(dec, float)
    r = radius_arcsec / 3600.0
    frames = []
    for k in range(0, len(ra), chunk):
        boxes = []
        for a, d in zip(ra[k:k + chunk], dec[k:k + chunk]):
            dra = r / max(np.cos(np.radians(d)), 1e-3)
            boxes.append(f"({ra_col} BETWEEN {a - dra} AND {a + dra} AND {dec_col} BETWEEN {d - r} AND {d + r})")
        adql = f'SELECT {", ".join(columns)} FROM "{table}" WHERE (' + " OR ".join(boxes) + f") {extra}"
        q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql))
        req = urllib.request.Request(TAP, data=q.encode(), headers=UA)     # POST: long OR lists
        txt = urllib.request.urlopen(req, timeout=900).read().decode()
        df = pd.read_csv(io.StringIO(txt))
        if len(df):
            frames.append(df)
        print(f"  positions_query {table}: {min(k + chunk, len(ra))}/{len(ra)}", flush=True)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=[c.strip('"') for c in columns])

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

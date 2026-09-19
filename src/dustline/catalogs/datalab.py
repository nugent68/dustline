"""NOIRLab Astro Data Lab TAP sync helper (stdlib only).

Note (learned the hard way): the Data Lab TAP rejects CONTAINS/q3c in some
configurations; use ``ra BETWEEN .. AND dec BETWEEN ..`` box queries and cut
to the circle client-side.
"""

from __future__ import annotations

import io
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

TAP = "https://datalab.noirlab.edu/tap/sync"
UA = {"User-Agent": "dustline"}


def box_query(table: str, columns: list[str], ra: float, dec: float, radius_deg: float,
              ra_col: str = "ra", dec_col: str = "dec", extra: str = "") -> pd.DataFrame:
    cosd = max(np.cos(np.radians(dec)), 1e-3)
    dra = radius_deg / cosd
    adql = (f"SELECT {', '.join(columns)} FROM {table} "
            f"WHERE {ra_col} BETWEEN {ra - dra} AND {ra + dra} "
            f"AND {dec_col} BETWEEN {dec - radius_deg} AND {dec + radius_deg} {extra}")
    q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql))
    req = urllib.request.Request(TAP + "?" + q, headers=UA)
    txt = urllib.request.urlopen(req, timeout=900).read().decode()
    df = pd.read_csv(io.StringIO(txt))
    # cut the box to the circle
    sep = np.hypot((df[ra_col] - ra) * cosd, df[dec_col] - dec)
    return df[sep <= radius_deg].reset_index(drop=True)

"""Local survey mirror: read catalogs from disk (e.g. at NERSC) instead of the network.

Set ``DUSTLINE_SURVEYS`` to the mirror root (e.g. /global/cfs/projectdirs/newera/surveys).  Each
survey directory holds a ``manifest.json`` (written by tools/mirror/) describing its tables,
keyed by the NETWORK table name the fetchers query (``gaiadr3.gaia_source``, ``decaps_dr2.object``,
``II/335/galex_ais`` ...).  The query helpers (datalab / vizier / gaia) ask ``table(name)`` first
and fall back to the network when it returns None, so the package runs anywhere and, with the
variable unset, behaves exactly as before.

``DUSTLINE_MIRROR``: ``auto`` (default: local when mirrored, else network), ``off`` (always the
network), ``only`` (never the network: a query for a table that is not mirrored raises
NetworkForbidden - for batch jobs that must not stall on a remote service).

Local frames reproduce the network ones: same column names and order, the network path's cone
geometry (Data Lab / VizieR: RA-Dec box then a flat-sky circle; ESA: great circle), and a final
CSV round trip so dtypes match what ``pd.read_csv`` makes of a TAP response.  Anything that
goes wrong reading a mirrored table raises MirrorError - it must never be cached as an empty
result the way a network failure is.
"""

from __future__ import annotations

import io
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import _healpix

ENV_ROOT, ENV_MODE = "DUSTLINE_SURVEYS", "DUSTLINE_MIRROR"
FINE_NSIDE = 4096


class MirrorError(RuntimeError):
    """A mirrored table could not be read (never swallowed into an empty cache)."""


class NetworkForbidden(MirrorError):
    """DUSTLINE_MIRROR=only and the table is not in the mirror."""


# ---------------------------------------------------------------------------------------- setup

def root() -> Path | None:
    r = os.environ.get(ENV_ROOT, "").strip()
    return Path(r) if r else None


def mode() -> str:
    m = os.environ.get(ENV_MODE, "auto").strip().lower() or "auto"
    if m not in ("auto", "off", "only"):
        raise ValueError(f"{ENV_MODE}={m!r}: expected auto, off or only")
    return m


_MANIFESTS: dict = {}


def _registry() -> dict:
    r = root()
    if r is None:
        return {}
    key = str(r)
    if key not in _MANIFESTS:
        reg = {}
        if r.is_dir():
            for mf in sorted(r.glob("*/manifest.json")):
                try:
                    d = json.loads(mf.read_text())
                except (OSError, json.JSONDecodeError) as e:
                    raise MirrorError(f"unreadable manifest {mf}: {e}") from e
                for name, spec in d.get("tables", {}).items():
                    reg[name] = (mf.parent, spec)
        _MANIFESTS[key] = reg
    return _MANIFESTS[key]


def table(name: str) -> "Table | None":
    """The mirrored table, or None (no mirror, mode off, not mirrored, or not complete)."""
    if root() is None or mode() == "off":
        return None
    hit = _registry().get(name)
    if hit is None or hit[1].get("status") != "complete":
        return None
    return Table(name, hit[0], hit[1])


def require_network(what: str) -> None:
    if root() is not None and mode() == "only":
        raise NetworkForbidden(f"{what}: not in the local mirror ({ENV_ROOT}={root()}) and {ENV_MODE}=only")


def render_adql(where) -> str:
    """Structured filters [(col, op, value)] -> the ADQL fragment the network path appends."""
    if not where:
        return ""
    parts = []
    for col, op, val in where:
        if isinstance(val, bool):
            v = "'t'" if val else "'f'"
        elif isinstance(val, str):
            v = f"'{val}'"
        else:
            v = repr(val)
        parts.append(f"{col} {'=' if op == '==' else op} {v}")
    return "AND " + " AND ".join(parts)


# ---------------------------------------------------------------------------------------- FITS access

class _Fits:
    """Memory-mapped FITS binary table: rows are read only where touched (the files are row-major,
    so reading one column of a whole bulge file would read gigabytes)."""

    def __init__(self, path: Path, hdu: int = 1):
        from astropy.io import fits
        self.path = Path(path)
        if not self.path.exists():
            raise MirrorError(f"mirror file missing: {self.path}")
        with fits.open(self.path, memmap=True, lazy_load_hdus=True) as hl:
            h = hl[hdu]                      # header only: touching h.data builds the whole record array
            hdr = h.header
            self.offset = h.fileinfo()["datLoc"]
            self.nrows = hdr["NAXIS2"]
            nf = hdr["TFIELDS"]
            self.tforms = {hdr[f"TTYPE{i}"]: hdr[f"TFORM{i}"].strip() for i in range(1, nf + 1)}
            # scaled integer columns (e.g. signed bytes stored as 'B' with TZERO = -128)
            self.scale = {hdr[f"TTYPE{i}"]: (hdr.get(f"TSCAL{i}", 1), hdr.get(f"TZERO{i}", 0)) for i in range(1, nf + 1)
                          if hdr.get(f"TSCAL{i}", 1) != 1 or hdr.get(f"TZERO{i}", 0) != 0}
            names = [hdr[f"TTYPE{i}"] for i in range(1, nf + 1)]
            self.raw_dtype = np.dtype([(n, _tform_dtype(self.tforms[n])) for n in names])
            if self.raw_dtype.itemsize != hdr["NAXIS1"]:
                raise MirrorError(f"{self.path.name}: row size {self.raw_dtype.itemsize} != NAXIS1 {hdr['NAXIS1']}")
        self.mm = np.memmap(self.path, dtype=self.raw_dtype, mode="r", offset=self.offset, shape=(self.nrows,)) \
            if self.nrows else np.zeros(0, self.raw_dtype)

    def _bisect(self, col: str, value, right: bool) -> int:
        """Row index where `value` would insert into the sorted column (reads ~log2(n) rows; np.searchsorted
        on a big-endian memmap byte-swaps the whole column first, i.e. reads the whole file)."""
        keys = self.mm[col]
        lo, hi = 0, self.nrows
        while lo < hi:
            mid = (lo + hi) // 2
            k = keys[mid]
            if k < value or (right and k == value):
                lo = mid + 1
            else:
                hi = mid
        return lo

    def column_ranges(self, col: str, lo_vals, hi_vals) -> list[tuple[int, int]]:
        """Row ranges whose sorted `col` lies in [lo, hi] for each pair."""
        out = []
        for a, b in zip(np.atleast_1d(lo_vals), np.atleast_1d(hi_vals)):
            i, j = self._bisect(col, a, False), self._bisect(col, b, True)
            if j > i:
                out.append((i, j))
        return out

    def read(self, cols: list[str], ranges=None) -> dict[str, np.ndarray]:
        if ranges is None:
            ranges = [(0, self.nrows)]
        out = {}
        for c in cols:
            if c not in self.raw_dtype.names:
                raise MirrorError(f"{self.path.name}: no column {c}")
            parts = [np.asarray(self.mm[c][a:b]) for a, b in ranges]
            v = np.concatenate(parts) if parts else np.zeros(0, self.raw_dtype[c])
            if v.dtype.kind in "iuf":
                v = v.astype(v.dtype.newbyteorder("="))
            if c in self.scale:
                sc, zero = self.scale[c]
                if v.dtype.kind in "iu" and float(sc) == 1 and float(zero) == int(zero):
                    v = v.astype(np.int64) + int(zero)
                else:
                    v = v * sc + zero
            if self.tforms.get(c, "").endswith("L"):
                v = v == ord("T")
            elif v.dtype.kind == "S":
                v = np.char.strip(v.astype("U"))
            out[c] = v
        return out


_TFORM = {"L": "i1", "B": "u1", "I": ">i2", "J": ">i4", "K": ">i8", "E": ">f4", "D": ">f8"}


def _tform_dtype(tform: str) -> np.dtype | tuple:
    m = re.match(r"^(\d*)([LXBIJKAEDCM])", tform)
    if not m:
        raise MirrorError(f"unsupported TFORM {tform}")
    rep, code = int(m.group(1) or 1), m.group(2)
    if code == "A":
        return np.dtype(f"S{rep}")
    if code not in _TFORM:
        raise MirrorError(f"unsupported TFORM {tform}")
    return (np.dtype(_TFORM[code]), (rep,)) if rep > 1 else np.dtype(_TFORM[code])


def _merge_runs(fine: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sorted pixel ids -> (starts, ends) of consecutive runs."""
    if not len(fine):
        return fine, fine
    brk = np.flatnonzero(np.diff(fine) != 1)
    starts = np.r_[fine[0], fine[brk + 1]]
    ends = np.r_[fine[brk], fine[-1]]
    return starts, ends


# ---------------------------------------------------------------------------------------- tables

_ALIAS = re.compile(r"^\s*(.+?)\s+AS\s+(\w+)\s*$", re.IGNORECASE)


def _parse_columns(columns) -> list[tuple[str, str]]:
    """['ra', 'ra2000 AS ra', '"E(B-V)"'] -> [(network name, output name)]."""
    out = []
    for c in columns:
        m = _ALIAS.match(c)
        src, dst = (m.group(1), m.group(2)) if m else (c, c)
        src, dst = src.strip().strip('"'), dst.strip().strip('"')
        out.append((src, dst))
    return out


_OPS = {"==": np.equal, "!=": np.not_equal, "<": np.less, "<=": np.less_equal, ">": np.greater, ">=": np.greater_equal}


class Table:
    """One mirrored table (see tools/mirror/ for the manifest schema)."""

    def __init__(self, name: str, base: Path, spec: dict):
        self.name, self.base, self.spec = name, Path(base), spec
        self.cols = spec.get("columns", {})
        self._single = None

    # -- helpers ------------------------------------------------------------------------
    def _path(self, **kw) -> Path:
        p = Path(self.spec.get("root") or self.base) / self.spec["files"].format(**kw)
        return p

    def _src(self, net: str) -> dict:
        c = self.cols.get(net)
        if c is None:
            raise MirrorError(f"{self.name}: column {net!r} is not in the mirror")
        return c

    def _finish(self, raw: dict[str, np.ndarray], wanted: list[tuple[str, str]], n: int) -> pd.DataFrame:
        data = {}
        for net, out in wanted:
            spec = self._src(net)
            if "chars" in spec:                 # a string built from per-band code columns (2MASS ph_qual)
                cmap = spec["map"]
                # codes stored as single raw bytes (e.g. b'\x04'); the manifest maps their ordinals
                parts = [np.vectorize(lambda x: cmap.get(str(ord(x[0])) if len(x) else "0", "?"), otypes=["U1"])(raw[c])
                         if len(raw[c]) else np.zeros(0, "U1") for c in spec["chars"]]
                data[out] = np.char.add(np.char.add(parts[0], parts[1]), parts[2]) if len(parts) == 3 else \
                    np.array(["".join(t) for t in zip(*parts)])
                continue
            v = raw[spec["src"]]
            if "index" in spec:
                v = v[:, spec["index"]]
            if spec.get("null_zero"):
                v = np.where(v == 0, np.nan, v.astype(float))
            if spec.get("null_value") is not None:
                v = np.where(v == spec["null_value"], np.nan, v.astype(float))
            if spec.get("bool_text"):
                v = np.asarray(v).astype(str) == "True"
            data[out] = v
        df = pd.DataFrame(data, index=range(n)) if n else pd.DataFrame({o: [] for _, o in wanted})
        # the network frames come from pd.read_csv: reproduce its dtypes exactly
        return pd.read_csv(io.StringIO(df.to_csv(index=False)), low_memory=False)

    def _raw_cols(self, wanted, where, extra=()) -> list[str]:
        need = set()
        for net, _ in wanted:
            s = self._src(net)
            need |= set(s["chars"]) if "chars" in s else {s["src"]}
        need |= {self._src(c)["src"] for c, _, _ in (where or [])}
        need |= set(extra)
        return sorted(need)

    def _apply_where(self, raw: dict, where) -> np.ndarray:
        n = len(next(iter(raw.values()))) if raw else 0
        keep = np.ones(n, bool)
        for col, op, val in where or []:
            spec = self._src(col)
            v = raw[spec["src"]]
            if spec.get("bool_text"):
                v = np.asarray(v).astype(str) == "True"
            if spec.get("null_zero"):
                v = np.where(v == 0, np.nan, v.astype(float))
            with np.errstate(invalid="ignore"):
                ok = _OPS[op](v, val)
            if np.asarray(v).dtype.kind == "f":
                ok &= np.isfinite(v)                    # SQL: NULL compares false
            keep &= ok
        return keep

    def _positional(self, ra: float, dec: float, radius_deg: float, raw_cols: list[str]) -> dict:
        """All rows of the files/pixels that may hold sources within radius_deg (conservative)."""
        kind = self.spec["partition"]
        out = {c: [] for c in raw_cols}
        pad = self.spec.get("sort_key", {}).get("pad_arcsec", 0) / 3600.0 if self.spec.get("sort_key") else 0.0
        if kind == "position":
            nside = self.spec["file_nside"]
            fine = _healpix.query_disc_fine(ra, dec, radius_deg + pad, nside_file=nside, nside_fine=FINE_NSIDE)
            have = set(self.spec.get("pixels", []))
            sk = self.spec.get("sort_key")
            for pix, fpix in fine.items():
                if have and pix not in have:
                    continue                        # outside the survey: legitimately empty
                f = _Fits(self._path(pix=pix), self.spec.get("hdu", 1))
                if sk:
                    s, e = _merge_runs(fpix)
                    sh = sk.get("shift", 0)
                    ranges = f.column_ranges(sk["col"], s << sh, ((e + 1) << sh) - 1)
                else:
                    ranges = None
                part = f.read(raw_cols, ranges)
                for c in raw_cols:
                    out[c].append(part[c])
        elif kind == "bbox":
            cosd = max(np.cos(np.radians(dec)), 1e-3)
            dra = min(radius_deg / cosd, 180.0)
            for fn, ra0, ra1, de0, de1 in self.spec["bboxes"]:
                if de1 < dec - radius_deg or de0 > dec + radius_deg:
                    continue
                # file RA interval [ra0, ra1] (ra0 < 0 when it straddles RA 0); query [ra-dra, ra+dra] mod 360
                if not any(ra1 >= ra + s - dra and ra0 <= ra + s + dra for s in (-360.0, 0.0, 360.0)):
                    continue
                part = _Fits(Path(self.spec.get("root") or self.base) / fn, self.spec.get("hdu", 1)).read(raw_cols)
                for c in raw_cols:
                    out[c].append(part[c])
        elif kind == "single":
            if self._single is None:
                self._single = _Fits(self._path(), self.spec.get("hdu", 1))
            part = self._single.read(raw_cols)
            for c in raw_cols:
                out[c].append(part[c])
        else:
            raise MirrorError(f"{self.name}: partition {kind!r} is not positional")
        return {c: (np.concatenate(v) if v else np.zeros(0)) for c, v in out.items()}

    # -- queries --------------------------------------------------------------------------
    def box(self, columns, ra, dec, radius_deg, ra_col="ra", dec_col="dec", where=None) -> pd.DataFrame:
        """Data Lab / VizieR semantics: RA-Dec box, then the flat-sky circle cut."""
        wanted = _parse_columns(columns)
        rac, decc = self._src(ra_col)["src"], self._src(dec_col)["src"]
        raw = self._positional(ra, dec, radius_deg, self._raw_cols(wanted, where, (rac, decc)))
        cosd = max(np.cos(np.radians(dec)), 1e-3)
        dra = radius_deg / cosd
        r_, d_ = raw[rac].astype(float), raw[decc].astype(float)
        keep = (r_ >= ra - dra) & (r_ <= ra + dra) & (d_ >= dec - radius_deg) & (d_ <= dec + radius_deg)
        keep &= self._apply_where(raw, where)
        keep &= np.hypot((r_ - ra) * cosd, d_ - dec) <= radius_deg
        return self._finish({c: v[keep] for c, v in raw.items()}, wanted, int(keep.sum()))

    def cone(self, columns, ra, dec, radius_deg, ra_col="ra", dec_col="dec", where=None) -> pd.DataFrame:
        """Great-circle cone (ESA CONTAINS, MAST cone)."""
        wanted = _parse_columns(columns)
        rac, decc = self._src(ra_col)["src"], self._src(dec_col)["src"]
        raw = self._positional(ra, dec, radius_deg, self._raw_cols(wanted, where, (rac, decc)))
        keep = _healpix._sep_deg(ra, dec, raw[rac].astype(float), raw[decc].astype(float)) <= radius_deg
        keep &= self._apply_where(raw, where)
        return self._finish({c: v[keep] for c, v in raw.items()}, wanted, int(keep.sum()))

    def positions(self, columns, ra, dec, radius_arcsec, ra_col="ra", dec_col="dec", where=None) -> pd.DataFrame:
        """Rows within radius of any position (the OR'd-box positions_query), per position."""
        frames = [self.box(columns, a, d, radius_arcsec / 3600.0 * 1.5, ra_col, dec_col, where)
                  for a, d in zip(np.atleast_1d(ra), np.atleast_1d(dec))]
        frames = [f for f in frames if len(f)]
        if not frames:
            return self._finish({}, _parse_columns(columns), 0)
        return pd.concat(frames, ignore_index=True).drop_duplicates().reset_index(drop=True)

    def by_ids(self, ids, raw_cols: list[str]) -> dict[str, np.ndarray]:
        """Raw columns for the given ids (id-partitioned tables); rows in the order of `ids` present."""
        ids = np.unique(np.asarray(ids, dtype=np.int64))
        kind = self.spec["partition"]
        idc = self.spec["id_col"]
        out = {c: [] for c in [idc] + [c for c in raw_cols if c != idc]}
        if kind == "source_id":
            groups = {}
            for i in ids:
                groups.setdefault(int(i) >> self.spec["id_shift"], []).append(i)
            files = {p: self._path(pix=p) for p in groups}
            have = set(self.spec.get("pixels", []))
            groups = {p: v for p, v in groups.items() if not have or p in have}
        elif kind == "source_id_range":
            sh = self.spec["range_shift"]
            rng = np.array([(a, b) for a, b, _ in self.spec["ranges"]])
            names = [f for _, _, f in self.spec["ranges"]]
            k = np.searchsorted(rng[:, 0], ids >> sh, "right") - 1
            groups, files = {}, {}
            for i, kk in zip(ids, k):
                if kk < 0 or (int(i) >> sh) > rng[kk, 1]:
                    raise MirrorError(f"{self.name}: source_id {i} outside the mirrored ranges")
                groups.setdefault(int(kk), []).append(i)
                files[int(kk)] = self._path(file=names[kk])
        else:
            raise MirrorError(f"{self.name}: partition {kind!r} has no id lookup")
        for g, members in groups.items():
            f = _Fits(files[g], self.spec.get("hdu", 1))
            m = np.asarray(members, np.int64)
            if f.nrows * f.raw_dtype.itemsize < (64 << 20) or len(m) > f.nrows // 20:
                part = f.read(list(out))                  # small file or dense request: one read + isin
                keep = np.isin(part[idc], m)
                part = {c: v[keep] for c, v in part.items()}
            else:
                part = f.read(list(out), f.column_ranges(idc, m, m))
            for c in out:
                out[c].append(part[c])
        return {c: (np.concatenate(v) if v else np.zeros(0)) for c, v in out.items()}


# ---------------------------------------------------------------------------------------- Gaia

def gaia_cone(ra: float, dec: float, radius_deg: float, columns: list[str], tmass_cols: dict) -> pd.DataFrame | None:
    """gaia_source in a great-circle cone + the 2MASS best-neighbour join, as the ESA query returns
    it (columns: `columns` then the 2MASS aliases in `tmass_cols` order).  None if gaia_source is
    not mirrored; the 2MASS columns are required whenever gaia_source is."""
    g = table("gaiadr3.gaia_source")
    if g is None:
        return None
    bn, tm = table("gaiadr3.tmass_psc_xsc_best_neighbour"), table("gaiadr1.tmass_original_valid")
    if bn is None or tm is None:
        raise MirrorError("gaia_source is mirrored but the 2MASS best-neighbour / PSC tables are not")
    df = g.cone(columns, ra, dec, radius_deg)
    if not len(df):
        raise MirrorError(f"local gaia_source cone ({ra}, {dec}, {radius_deg}) is empty")
    x = bn.by_ids(df.source_id.values, ["original_ext_source_id", "angular_distance"])
    xm = pd.DataFrame({"source_id": x["source_id"], "designation": x["original_ext_source_id"],
                       "tmass_sep": x["angular_distance"]})
    # the 2MASS rows: a padded box over the field, joined by designation
    pad = radius_deg + 10.0 / 3600.0
    t = tm.box(["designation"] + [c for c in tmass_cols.values() if c not in ("designation", "tmass_sep")],
               ra, dec, pad, ra_col="ra", dec_col="dec")
    t = t.drop_duplicates("designation")
    m = xm.merge(t, on="designation", how="left")
    out = df.merge(m.drop(columns="designation"), on="source_id", how="left")
    out = out.rename(columns={v: k for k, v in tmass_cols.items()})
    order = list(df.columns) + [k for k in tmass_cols if k != "tmass_sep"] + ["tmass_sep"]
    return pd.read_csv(io.StringIO(out[order].to_csv(index=False)), low_memory=False)


XP_COLS = ["source_id", "solution_id"] + [
    f"{b}_{c}" for b in ("bp", "rp") for c in (
        "basis_function_id", "degrees_of_freedom", "n_parameters", "n_measurements", "n_rejected_measurements",
        "standard_deviation", "chi_squared", "coefficients", "coefficient_errors", "coefficient_correlations",
        "n_relevant_bases", "relative_shrinking")]


def _fmt_array(v: np.ndarray) -> str:
    if v.dtype == np.float64:
        return "(" + ", ".join(repr(float(x)) for x in v) + ")"
    return "(" + ", ".join(str(np.float32(x)) for x in v) + ")"


def xp_csv_rows(ids) -> str | None:
    """The DataLink XP_CONTINUOUS RAW csv text (header + rows) for `ids` from the mirror; None if
    XP is not mirrored.  Every id must be present."""
    t = table("gaiadr3.xp_continuous_mean_spectrum")
    if t is None:
        return None
    raw = t.by_ids(ids, XP_COLS[1:] + ["has_bp", "has_rp"])
    got = set(int(i) for i in raw["source_id"])
    missing = [int(i) for i in ids if int(i) not in got]
    if missing:
        raise MirrorError(f"XP mirror lacks {len(missing)} of {len(ids)} requested sources (e.g. {missing[:3]})")
    lines = [",".join(XP_COLS)]
    for k in range(len(raw["source_id"])):
        row = []
        for c in XP_COLS:
            band = c[:2] if c[:3] in ("bp_", "rp_") else None
            if band and not raw[f"has_{band}"][k]:
                row.append("")
                continue
            v = raw[c][k]
            if isinstance(v, np.ndarray) and v.ndim == 1:
                row.append('"' + _fmt_array(v) + '"')
            elif np.issubdtype(type(v), np.floating):
                row.append(str(np.float32(v)) if raw[c].dtype == np.float32 else repr(float(v)))
            else:
                row.append(str(int(v)))
        lines.append(",".join(row))
    return "\n".join(lines) + "\n"

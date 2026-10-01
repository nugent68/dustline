"""Local survey mirror (dustline.catalogs.local): synthetic mirrors in a temp DUSTLINE_SURVEYS,
checked against the network code path (urlopen mocked to serve the same synthetic catalog)."""

from __future__ import annotations

import io
import json
import os
import urllib.parse
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from dustline.catalogs import _healpix as H
from dustline.catalogs import datalab, local

RA0, DEC0 = 267.86642, -33.13517


# ------------------------------------------------------------------------------------ healpix

def test_healpix_known_pixels_and_roundtrip():
    assert H.ang2pix_nest(16, RA0, DEC0)[0] == 1795            # the MMU Gaia file holding 0095
    assert H.ang2pix_nest(16, 150.12, 2.21)[0] == 1703          # COSMOS
    rng = np.random.default_rng(1)
    ra = rng.uniform(0, 360, 20000); dec = np.degrees(np.arcsin(rng.uniform(-1, 1, 20000)))
    for ns in (1, 8, 32, 4096):
        p = H.ang2pix_nest(ns, ra, dec)
        assert p.min() >= 0 and p.max() < 12 * ns * ns
        cra, cdec = H.pix2ang_nest(ns, p)
        assert np.all(H.ang2pix_nest(ns, cra, cdec) == p)
        assert np.all(H._sep_deg(ra, dec, cra, cdec) <= H.max_pixrad_deg(ns))


@pytest.mark.parametrize("ra,dec,r", [(RA0, DEC0, 0.5), (0.05, 1.0, 0.8), (359.95, -20, 0.3), (40, 89.9, 0.4)])
def test_query_disc_is_complete(ra, dec, r):
    rng = np.random.default_rng(2)
    n = 5000
    th = np.radians(r) * np.sqrt(rng.uniform(0, 1, n)); ph = rng.uniform(0, 2 * np.pi, n)
    dd = np.clip(dec + np.degrees(th) * np.sin(ph), -90, 90)
    rr = (ra + np.degrees(th) * np.cos(ph) / np.maximum(np.cos(np.radians(dd)), 1e-3)) % 360
    inside = H._sep_deg(ra, dec, rr, dd) <= r
    fine = np.concatenate(list(H.query_disc_fine(ra, dec, r).values()))
    assert np.isin(H.ang2pix_nest(4096, rr[inside], dd[inside]), fine).all()


# ------------------------------------------------------------------------------------ fixtures

def _catalog(n=4000, seed=3, ra=RA0, dec=DEC0, r=0.25):
    rng = np.random.default_rng(seed)
    cosd = np.cos(np.radians(dec))
    c = pd.DataFrame({"ra": ra + rng.uniform(-r, r, n) / cosd, "dec": dec + rng.uniform(-r, r, n)})
    c["mag_g"] = rng.uniform(14, 22, n).astype("f4")
    c["err_g"] = np.where(rng.uniform(0, 1, n) < 0.1, np.nan, rng.uniform(0.001, 0.3, n)).astype("f4")
    c["type"] = rng.choice(["PSF", "REX", "DEV"], n)
    c["snr"] = rng.uniform(0, 100, n).astype("f4")
    c["id"] = np.arange(n, dtype=np.int64) + 10**12
    return c


def _write_position_table(root: Path, survey: str, name: str, cat: pd.DataFrame, status="complete"):
    d = root / survey
    fine = H.ang2pix_nest(4096, cat.ra.values, cat.dec.values)
    cat = cat.assign(_hpx=fine).sort_values("_hpx")
    fpix = cat._hpx.values >> (2 * 7)
    pixels = []
    for p in np.unique(fpix):
        part = cat[fpix == p]
        cols = [fits.Column(name="RA", format="D", array=part.ra.values), fits.Column(name="DEC", format="D", array=part.dec.values),
                fits.Column(name="MAG_G", format="E", array=part.mag_g.values), fits.Column(name="ERR_G", format="E", array=part.err_g.values),
                fits.Column(name="TYPE", format="3A", array=part.type.values.astype("S3")),
                fits.Column(name="SNR", format="E", array=part.snr.values), fits.Column(name="ID", format="K", array=part.id.values),
                fits.Column(name="_HPX", format="K", array=part._hpx.values)]
        (d / "hp32").mkdir(parents=True, exist_ok=True)
        fits.BinTableHDU.from_columns(cols).writeto(d / "hp32" / f"{p:05d}.fits", overwrite=True)
        pixels.append(int(p))
    spec = {"status": status, "format": "fits", "hdu": 1, "files": "hp32/{pix:05d}.fits", "partition": "position",
            "file_nside": 32, "nest": True, "pixels": pixels, "sort_key": {"col": "_HPX", "shift": 0, "pad_arcsec": 0},
            "columns": {"ra": {"src": "RA"}, "dec": {"src": "DEC"}, "mag_g": {"src": "MAG_G"}, "err_g": {"src": "ERR_G"},
                        "type": {"src": "TYPE"}, "snr": {"src": "SNR"}, "id": {"src": "ID"}}}
    mf = d / "manifest.json"
    old = json.loads(mf.read_text()) if mf.exists() else {"schema": 1, "survey": survey, "tables": {}}
    old["tables"][name] = spec
    mf.write_text(json.dumps(old))
    return pixels


@pytest.fixture
def mirror(tmp_path, monkeypatch):
    monkeypatch.setenv("DUSTLINE_SURVEYS", str(tmp_path))
    monkeypatch.delenv("DUSTLINE_MIRROR", raising=False)
    local._MANIFESTS.clear()
    cat = _catalog()
    _write_position_table(tmp_path, "synth", "synth.object", cat)
    yield tmp_path, cat
    local._MANIFESTS.clear()


def _fake_network(cat: pd.DataFrame, calls: list):
    """urlopen stand-in serving `cat` with the Data Lab box semantics the ADQL asks for."""
    def urlopen(req, timeout=None):
        url = req.full_url if hasattr(req, "full_url") else req
        q = urllib.parse.parse_qs(urllib.parse.urlparse(url).query or (req.data or b"").decode())["QUERY"][0]
        calls.append(q)
        import re
        cols = [c.strip() for c in re.match(r"SELECT (.*?) FROM", q).group(1).split(",")]
        (ra_lo, ra_hi), (de_lo, de_hi) = [tuple(map(float, m)) for m in re.findall(r"BETWEEN (\S+) AND (\S+)", q)]
        d = cat[(cat.ra >= ra_lo) & (cat.ra <= ra_hi) & (cat.dec >= de_lo) & (cat.dec <= de_hi)]
        if "type = 'PSF'" in q:
            d = d[d.type == "PSF"]
        if "snr > 30" in q:
            d = d[d.snr > 30]
        out = pd.DataFrame({c.split(" AS ")[-1]: d[c.split(" AS ")[0]].values for c in cols})

        class R:
            def read(self):
                return out.to_csv(index=False).encode()
        return R()
    return urlopen


# ------------------------------------------------------------------------------------ behaviour

def test_unset_means_network_and_no_filesystem(monkeypatch, tmp_path):
    monkeypatch.delenv("DUSTLINE_SURVEYS", raising=False)
    local._MANIFESTS.clear()
    assert local.table("synth.object") is None
    local.require_network("anything")                     # no raise without a mirror


def test_box_parity_with_network(mirror, monkeypatch):
    root, cat = mirror
    cols = ["ra", "dec", "mag_g", "err_g", "id"]
    loc = datalab.box_query("synth.object", cols, RA0, DEC0, 0.2)
    calls = []
    monkeypatch.setattr(datalab.urllib.request, "urlopen", _fake_network(cat, calls))
    monkeypatch.setenv("DUSTLINE_MIRROR", "off")
    net = datalab.box_query("synth.object", cols, RA0, DEC0, 0.2)
    assert calls, "network path was not taken with DUSTLINE_MIRROR=off"
    pd.testing.assert_frame_equal(loc.sort_values("id").reset_index(drop=True),
                                  net.sort_values("id").reset_index(drop=True), check_exact=False, rtol=1e-6)
    assert 0 < len(loc) < len(cat)


def test_where_filters_and_aliases(mirror, monkeypatch):
    root, cat = mirror
    where = [("type", "==", "PSF"), ("snr", ">", 30)]
    loc = datalab.box_query("synth.object", ["ra", "dec", "type", "snr AS s", "id"], RA0, DEC0, 0.2, where=where)
    assert list(loc.columns) == ["ra", "dec", "type", "s", "id"]
    assert set(loc.type) == {"PSF"} and (loc.s > 30).all()
    calls = []
    monkeypatch.setattr(datalab.urllib.request, "urlopen", _fake_network(cat, calls))
    monkeypatch.setenv("DUSTLINE_MIRROR", "off")
    net = datalab.box_query("synth.object", ["ra", "dec", "type", "snr AS s", "id"], RA0, DEC0, 0.2, where=where)
    assert calls[0].rstrip().endswith("AND type = 'PSF' AND snr > 30")
    assert sorted(net.id) == sorted(loc.id)


def test_null_compares_false(mirror):
    root, cat = mirror
    d = datalab.box_query("synth.object", ["err_g", "id"], RA0, DEC0, 0.2, where=[("err_g", "<", 1.0)])
    assert d.err_g.notna().all()


def test_render_adql_matches_fetcher_strings():
    from dustline.catalogs import apogee, desi
    assert local.render_adql(desi._WHERE) == desi._EXTRA
    assert local.render_adql(apogee._WHERE) == apogee._EXTRA
    assert local.render_adql([("type", "==", "PSF")]) == "AND type = 'PSF'"


def test_raw_extra_refused_locally(mirror):
    with pytest.raises(local.MirrorError):
        datalab.box_query("synth.object", ["ra"], RA0, DEC0, 0.1, extra="AND snr > 3")


def test_missing_file_raises_and_is_not_cached(mirror, tmp_path):
    root, cat = mirror
    for f in (root / "synth" / "hp32").glob("*.fits"):
        f.unlink()
    with pytest.raises(local.MirrorError):
        datalab.box_query("synth.object", ["ra"], RA0, DEC0, 0.1)


def test_partial_table_falls_back_and_only_mode_forbids(mirror, monkeypatch):
    root, cat = mirror
    _write_position_table(root, "synth", "synth.object", cat, status="partial")
    local._MANIFESTS.clear()
    assert local.table("synth.object") is None
    calls = []
    monkeypatch.setattr(datalab.urllib.request, "urlopen", _fake_network(cat, calls))
    datalab.box_query("synth.object", ["ra", "dec"], RA0, DEC0, 0.1)
    assert calls
    monkeypatch.setenv("DUSTLINE_MIRROR", "only")
    with pytest.raises(local.NetworkForbidden):
        datalab.box_query("synth.object", ["ra", "dec"], RA0, DEC0, 0.1)


def test_fetcher_does_not_cache_empty_on_mirror_error(mirror, tmp_path, monkeypatch):
    from dustline.catalogs import decaps
    monkeypatch.setattr(decaps, "box_query", lambda *a, **k: (_ for _ in ()).throw(local.MirrorError("boom")))

    class WS:
        ra, dec, radius_arcmin = RA0, DEC0, 5.0
        def path(self, n):
            return tmp_path / n
    with pytest.raises(local.MirrorError):
        decaps.fetch(WS(), pd.DataFrame(dict(source_id=[], ra=[], dec=[])))
    assert not (tmp_path / "decaps_gaia.csv").exists()


# ------------------------------------------------------------------------------------ XP

def _xp_rows_from_cache():
    for p in Path(os.path.expanduser("~/.cache/dustline/sightlines")).glob("*/xp_continuous_raw.csv"):
        return p
    return None


def test_xp_csv_round_trip(tmp_path, monkeypatch):
    src = _xp_rows_from_cache()
    if src is None:
        pytest.skip("no cached xp_continuous_raw.csv")
    raw = pd.read_csv(src, nrows=5)
    arr = lambda s, dt: np.array([float(x) for x in s.strip("()[]").split(",")], dtype=dt)  # noqa: E731
    cols = [fits.Column(name="source_id", format="K", array=raw.source_id.values),
            fits.Column(name="solution_id", format="K", array=raw.solution_id.values)]
    for c in local.XP_COLS[2:]:
        if c.endswith(("coefficients",)):
            cols.append(fits.Column(name=c, format="55D", array=np.stack([arr(s, "f8") for s in raw[c]])))
        elif c.endswith(("coefficient_errors",)):
            cols.append(fits.Column(name=c, format="55E", array=np.stack([arr(s, "f4") for s in raw[c]])))
        elif c.endswith("coefficient_correlations"):
            cols.append(fits.Column(name=c, format="1485E", array=np.stack([arr(s, "f4") for s in raw[c]])))
        elif raw[c].dtype.kind == "f":
            cols.append(fits.Column(name=c, format="E", array=raw[c].values.astype("f4")))
        else:
            cols.append(fits.Column(name=c, format="I", array=raw[c].values.astype("i2")))
    cols += [fits.Column(name="has_bp", format="L", array=np.ones(5, bool)),
             fits.Column(name="has_rp", format="L", array=np.ones(5, bool))]
    order = np.argsort(raw.source_id.values)
    d = tmp_path / "gaia_dr3"
    (d / "xp").mkdir(parents=True)
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.data = hdu.data[order]
    hdu.writeto(d / "xp" / "part.fits")
    sid = np.sort(raw.source_id.values)
    lo, hi = int(sid[0] >> 43), int(sid[-1] >> 43)
    (d / "manifest.json").write_text(json.dumps({"schema": 1, "survey": "gaia_dr3", "tables": {
        "gaiadr3.xp_continuous_mean_spectrum": {"status": "complete", "format": "fits", "hdu": 1,
                                                "partition": "source_id_range", "files": "xp/{file}", "range_shift": 43,
                                                "ranges": [[lo, hi, "part.fits"]], "id_col": "source_id"}}}))
    monkeypatch.setenv("DUSTLINE_SURVEYS", str(tmp_path))
    local._MANIFESTS.clear()
    txt = local.xp_csv_rows(raw.source_id.values)
    back = pd.read_csv(io.StringIO(txt)).set_index("source_id").loc[raw.source_id.values].reset_index()
    for c in local.XP_COLS:
        if not pd.api.types.is_numeric_dtype(back[c]):
            dt = "f8" if c.endswith("coefficients") else "f4"
            assert all(np.array_equal(arr(a, dt), arr(b, dt)) for a, b in zip(back[c], raw[c])), c
        else:
            assert np.allclose(back[c].values, raw[c].values, rtol=1e-7), c
    gaiaxpy = pytest.importorskip("gaiaxpy")
    f1, _ = gaiaxpy.calibrate(raw, save_file=False, truncation=False)
    f2, _ = gaiaxpy.calibrate(back, save_file=False, truncation=False)
    assert all(np.array_equal(np.asarray(a), np.asarray(b)) for a, b in zip(f1.flux, f2.flux))
    with pytest.raises(local.MirrorError):
        local.xp_csv_rows([int(sid[-1]) + 1])
    local._MANIFESTS.clear()


def test_postgres_nan_semantics_and_replace(mirror):
    """Data Lab (PostgreSQL) keeps NaN rows under `x > c`; SQL NULL logic drops them."""
    root, cat = mirror
    mf = root / "synth" / "manifest.json"
    d = json.loads(mf.read_text())
    d["tables"]["synth.object"]["nan_semantics"] = "postgres"
    d["tables"]["synth.object"]["columns"]["id"]["replace"] = [int(cat.id.min()), -1]
    mf.write_text(json.dumps(d))
    local._MANIFESTS.clear()
    pg = datalab.box_query("synth.object", ["err_g", "id"], RA0, DEC0, 0.2, where=[("err_g", ">", 0.0)])
    assert pg.err_g.isna().any()
    assert (pg.err_g.dropna() > 0).all()
    lt = datalab.box_query("synth.object", ["err_g", "id"], RA0, DEC0, 0.2, where=[("err_g", "<", 1.0)])
    assert lt.err_g.notna().all()
    allrows = datalab.box_query("synth.object", ["id"], RA0, DEC0, 0.5)
    assert int(cat.id.min()) not in set(allrows.id) and (-1 in set(allrows.id) or len(allrows) < len(cat))

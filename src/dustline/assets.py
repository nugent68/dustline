"""Model-asset management: download-on-first-use of the large model files.

The PHOENIX NewEra spectral cache (~91 MB) and the MIST v1.2 isochrone tables
are too large to ship inside the wheel; they are attached to a GitHub release
of the dustline repository and fetched on first use into the user cache
directory (default ``~/.cache/dustline/``, override with $DUSTLINE_CACHE_DIR).

Filter transmission curves are small and ship as package data
(``dustline/data/filters/*.dat``).

Usage:
    from dustline import assets
    path = assets.fetch("newera_full_cache.npz")   # downloads once, verifies sha256
"""

from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Registry of release assets: name -> (release tag, sha256, size bytes).
# A sha256 of None disables verification for that asset (with a warning).
# ---------------------------------------------------------------------------
RELEASE_BASE = "https://github.com/nugent68/dustline/releases/download"

ASSETS: dict[str, dict] = {
    # 2,517 PHOENIX NewEra models, [M/H] -0.5..+0.5, log g 0-6, 2500-25000 A at 2 A,
    # surface flux in W m^-2 nm^-1 (see dustline.models for the convention checks).
    "newera_full_cache.npz": dict(
        tag="v0.1.0",
        sha256="c4cc1bf53863a375a568b7c53408539555aedbeca001c291c753d9949169bce3",
        size=91_425_303),
    # MIST v1.2 vvcrit0.4 basic isochrones, [Fe/H] -0.5..+0.5, repacked npz used by
    # the radius prior (see dustline.models.radius_prior; tools/build_mist_npz.py).
    "mist_v1.2_basic.npz": dict(
        tag="v0.1.0",
        sha256="9ac0d3be95fda3373ad9aa3bcab54c2f4ee59159e2e85b2a42cec9aaac0f7062",
        size=13_773_902),
    # 4,366 NewEra models from the HSR files' LSR spectra: [M/H] -2..+0.5, log g 0-6,
    # Teff >= 3200, 900-25000 A at 2 A + 25000-60000 A at 20 A (GALEX + WISE coverage;
    # tools/build_newera_uvir_cache.py).
    "newera_uvir_cache.npz": dict(
        tag="v0.4.0",
        sha256="2e03d4f0c101e848faf9738132ea8db3af2f0d05afbf81dff1e2c38f2c6c8918",
        size=214_319_883),
    # Empirical dwarf-template corrections per DESI T_eff / [Fe/H] bin from 2,451
    # unreddened DESI x XP dwarfs (dustline.calib; tools/build_template_corrections.py).
    "template_corrections.npz": dict(
        tag="v0.5.0",
        sha256="b7956d09e743edad0e68ebb91b1bd8ad926eb94eaa17c80db31faa29968a071f",
        size=120060),
}


def cache_dir() -> Path:
    """The local asset cache directory (created on demand)."""
    d = os.environ.get("DUSTLINE_CACHE_DIR")
    if d:
        p = Path(d)
    else:
        p = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "dustline"
    p.mkdir(parents=True, exist_ok=True)
    return p


def _sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _download(url: str, dest: Path, expected_size: int | None = None) -> None:
    """Stream url to dest atomically, with a progress bar when tqdm is available."""
    try:
        from tqdm import tqdm
    except ImportError:  # pragma: no cover
        tqdm = None
    r = requests.get(url, stream=True, timeout=60)
    r.raise_for_status()
    total = int(r.headers.get("content-length", 0)) or expected_size or 0
    fd, tmp = tempfile.mkstemp(dir=dest.parent, suffix=".part")
    bar = tqdm(total=total, unit="B", unit_scale=True, desc=dest.name) if tqdm and total else None
    try:
        with os.fdopen(fd, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                if bar:
                    bar.update(len(chunk))
        shutil.move(tmp, dest)
    finally:
        if bar:
            bar.close()
        if os.path.exists(tmp):
            os.unlink(tmp)


def fetch(name: str, force: bool = False) -> Path:
    """Return the local path of a release asset, downloading it on first use.

    Verifies the sha256 recorded in ASSETS (skipped, with a warning, while the
    registry hash is still None). ``force=True`` re-downloads.
    """
    if name not in ASSETS:
        raise KeyError(f"unknown asset {name!r}; known: {sorted(ASSETS)}")
    info = ASSETS[name]
    dest = cache_dir() / name
    if dest.exists() and not force:
        if info["sha256"] and _sha256(dest) != info["sha256"]:
            raise RuntimeError(
                f"{dest} exists but its sha256 does not match the registry; "
                f"delete it or run fetch({name!r}, force=True)"
            )
        return dest
    url = f"{RELEASE_BASE}/{info['tag']}/{name}"
    size = f"~{info['size'] / 1e6:.0f} MB" if info.get("size") else "size unknown"
    print(f"dustline: downloading {name} ({size}) from {url}")
    _download(url, dest, info["size"])
    if info["sha256"]:
        got = _sha256(dest)
        if got != info["sha256"]:
            dest.unlink()
            raise RuntimeError(f"sha256 mismatch for {name}: got {got}, expected {info['sha256']}")
    else:
        import warnings

        warnings.warn(f"asset {name} has no registry sha256 yet; downloaded without verification")
    return dest


def filter_dir() -> Path:
    """Directory of the packaged filter transmission curves."""
    return Path(__file__).parent / "data" / "filters"

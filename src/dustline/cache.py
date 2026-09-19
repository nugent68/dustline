"""Per-sightline workspace: where a sightline's downloaded catalogs, XP spectra
and fit products live so reruns are instant.

A sightline is keyed by its position rounded to 0.1" plus the radius and a hash
of the pipeline configuration; the directory sits under the same user cache root
as the model assets (~/.cache/dustline/sightlines/<key>/).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from . import assets


class Workspace:
    """Directory manager for one sightline run."""

    def __init__(self, ra: float, dec: float, radius_arcmin: float, config: dict | None = None):
        self.ra, self.dec, self.radius_arcmin = float(ra), float(dec), float(radius_arcmin)
        self.config = dict(config or {})
        blob = json.dumps(self.config, sort_keys=True).encode()
        chash = hashlib.sha256(blob).hexdigest()[:8]
        key = f"ra{self.ra:+011.5f}_dec{self.dec:+010.5f}_r{self.radius_arcmin:g}_{chash}"
        self.dir = assets.cache_dir() / "sightlines" / key
        self.dir.mkdir(parents=True, exist_ok=True)
        cfg_path = self.dir / "config.json"
        if not cfg_path.exists():
            cfg_path.write_text(json.dumps(
                dict(ra=self.ra, dec=self.dec, radius_arcmin=self.radius_arcmin, **self.config),
                indent=2))

    def path(self, name: str) -> Path:
        return self.dir / name

    def has(self, name: str) -> bool:
        return (self.dir / name).exists()

    def __repr__(self) -> str:  # pragma: no cover
        return f"Workspace({self.dir})"

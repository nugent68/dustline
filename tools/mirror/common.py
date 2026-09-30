"""Shared helpers for the survey-mirror ingest scripts (tools/mirror/*.py).

Every ingest script is idempotent and resumable: each unit of work (one input file, one sky
pixel) appends one JSON line to a ledger when it is complete, and a rerun skips units already in
the ledger.  Outputs are written to a temporary name and renamed into place, so a killed job
never leaves a truncated product that looks finished.  Nothing here imports dustline, so the
scripts run in any environment with numpy / astropy.
"""

from __future__ import annotations

import hashlib
import json
import os
import socket
import time
from pathlib import Path

SURVEYS = Path(os.environ.get("DUSTLINE_SURVEYS", "/global/cfs/projectdirs/newera/surveys"))
STAGING = Path(os.environ.get("MIRROR_STAGING", os.path.join(os.environ.get("PSCRATCH", "/tmp"), "dustline_mirror_staging")))


def md5sum(path: Path, chunk: int = 1 << 24) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


class Ledger:
    """Append-only JSONL record of finished work units; safe for concurrent array tasks
    (each line is written with a single O_APPEND write)."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def done(self) -> dict[str, dict]:
        out = {}
        if self.path.exists():
            for line in self.path.read_text().splitlines():
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue            # a torn last line from a killed task: that unit is redone
                out[rec["unit"]] = rec
        return out

    def add(self, unit: str, **rec) -> None:
        rec = dict(unit=unit, host=socket.gethostname(), time=time.strftime("%Y-%m-%dT%H:%M:%S"), **rec)
        line = (json.dumps(rec, sort_keys=True) + "\n").encode()
        fd = os.open(self.path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o664)
        try:
            os.write(fd, line)
        finally:
            os.close(fd)


def claim(lockdir: Path, unit: str, stale_s: float = 3 * 3600) -> bool:
    """Atomically claim a work unit (O_EXCL lock file) so concurrent runners never process the
    same unit; a lock older than stale_s (a killed runner) is taken over."""
    lockdir = Path(lockdir)
    lockdir.mkdir(parents=True, exist_ok=True)
    p = lockdir / f"{unit}.lock"
    try:
        fd = os.open(p, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o664)
    except FileExistsError:
        try:
            if time.time() - p.stat().st_mtime < stale_s:
                return False
            p.unlink()
            return claim(lockdir, unit, stale_s)
        except FileNotFoundError:
            return claim(lockdir, unit, stale_s)
    os.write(fd, f"{socket.gethostname()} {os.getpid()} {time.time():.0f}\n".encode())
    os.close(fd)
    return True


def release(lockdir: Path, unit: str) -> None:
    (Path(lockdir) / f"{unit}.lock").unlink(missing_ok=True)


def atomic_path(final: Path) -> Path:
    """A temporary sibling of `final`; os.replace(tmp, final) publishes it."""
    final = Path(final)
    final.parent.mkdir(parents=True, exist_ok=True)
    return final.with_name(f".{final.name}.tmp{os.getpid()}")


def my_units(units: list, task_id: int | None, ntasks: int | None) -> list:
    """Round-robin share of the work for a Slurm array task (SLURM_ARRAY_TASK_ID / _COUNT by default)."""
    if task_id is None:
        task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    if ntasks is None:
        ntasks = int(os.environ.get("SLURM_ARRAY_TASK_COUNT", 1))
    return units[task_id::ntasks]


def write_manifest(survey_dir: Path, tables: dict, survey: str) -> Path:
    """surveys/<survey>/manifest.json (schema 1): {"schema", "survey", "tables": {network_name: {...}}}."""
    path = Path(survey_dir) / "manifest.json"
    old = json.loads(path.read_text()) if path.exists() else {"schema": 1, "survey": survey, "tables": {}}
    old["tables"].update(tables)
    tmp = atomic_path(path)
    tmp.write_text(json.dumps(old, indent=1, sort_keys=True))
    os.replace(tmp, path)
    return path

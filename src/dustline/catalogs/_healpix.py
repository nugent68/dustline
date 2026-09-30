"""Minimal HEALPix (NESTED scheme) in numpy: enough for the local survey mirror.

The mirror partitions every positional table into HEALPix nside-32 files and sorts rows within a
file by their nside-4096 index (= Gaia's level-12 index, source_id >> 35), so a cone read is:
file pixels overlapping the disc -> fine pixels overlapping the disc -> row ranges by
searchsorted -> exact angular cut.  Only ang2pix / pix2ang (nested) and a conservative
hierarchical disc query are needed, so healpy / astropy-healpix are not dependencies.
"""

from __future__ import annotations

import numpy as np

_JRLL = np.array([2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4])
_JPLL = np.array([1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7])


def _spread_bits(v: np.ndarray) -> np.ndarray:
    """Put the bits of v into the even bit positions (v < 2**29)."""
    v = v.astype(np.int64)
    out = np.zeros_like(v)
    for b in range(30):
        out |= ((v >> b) & 1) << (2 * b)
    return out


def _compress_bits(v: np.ndarray) -> np.ndarray:
    """Inverse of _spread_bits: take the even bits of v."""
    v = v.astype(np.int64)
    out = np.zeros_like(v)
    for b in range(30):
        out |= ((v >> (2 * b)) & 1) << b
    return out


def ang2pix_nest(nside: int, ra_deg, dec_deg) -> np.ndarray:
    """Nested pixel index of (ra, dec) in degrees."""
    ra = np.atleast_1d(np.asarray(ra_deg, dtype=np.float64))
    dec = np.atleast_1d(np.asarray(dec_deg, dtype=np.float64))
    z = np.sin(np.radians(dec))
    tt = np.mod(np.radians(ra), 2 * np.pi) / (np.pi / 2)          # [0, 4)
    tt = np.where(tt >= 4.0, 0.0, tt)
    za = np.abs(z)
    face = np.empty(ra.shape, np.int64); ix = np.empty(ra.shape, np.int64); iy = np.empty(ra.shape, np.int64)
    eq = za <= 2.0 / 3.0
    if eq.any():
        t1 = nside * (0.5 + tt[eq]); t2 = nside * z[eq] * 0.75
        jp = np.floor(t1 - t2).astype(np.int64); jm = np.floor(t1 + t2).astype(np.int64)
        ifp = jp // nside; ifm = jm // nside
        face[eq] = np.where(ifp == ifm, ifp | 4, np.where(ifp < ifm, ifp, ifm + 8))
        ix[eq] = jm % nside
        iy[eq] = nside - (jp % nside) - 1
    po = ~eq
    if po.any():
        ntt = np.minimum(np.floor(tt[po]).astype(np.int64), 3)
        tp = tt[po] - ntt
        tmp = nside * np.sqrt(3.0 * (1.0 - za[po]))
        jp = np.minimum(np.floor(tp * tmp).astype(np.int64), nside - 1)
        jm = np.minimum(np.floor((1.0 - tp) * tmp).astype(np.int64), nside - 1)
        north = z[po] >= 0
        face[po] = np.where(north, ntt, ntt + 8)
        ix[po] = np.where(north, nside - jm - 1, jp)
        iy[po] = np.where(north, nside - jp - 1, jm)
    return face * nside * nside + (_spread_bits(ix) | (_spread_bits(iy) << 1))


def pix2ang_nest(nside: int, pix) -> tuple[np.ndarray, np.ndarray]:
    """(ra, dec) in degrees of nested pixel centres."""
    pix = np.atleast_1d(np.asarray(pix, dtype=np.int64))
    npface = nside * nside
    face = pix // npface
    ipf = pix % npface
    ix = _compress_bits(ipf); iy = _compress_bits(ipf >> 1)
    jr = _JRLL[face] * nside - ix - iy - 1
    nl4 = 4 * nside
    fact2 = 4.0 / (12.0 * npface)
    nr = np.where(jr < nside, jr, np.where(jr > 3 * nside, nl4 - jr, nside))
    z = np.where(jr < nside, 1.0 - nr * nr * fact2,
                 np.where(jr > 3 * nside, -1.0 + nr * nr * fact2, (2 * nside - jr) * (2.0 / (3.0 * nside))))
    kshift = np.where((jr >= nside) & (jr <= 3 * nside), (jr - nside) & 1, 0)
    jp = (_JPLL[face] * nr + ix - iy + 1 + kshift) // 2
    jp = np.where(jp > nl4, jp - nl4, np.where(jp < 1, jp + nl4, jp))
    phi = (jp - (kshift + 1) * 0.5) * (np.pi / 2 / nr)
    return np.degrees(phi) % 360.0, np.degrees(np.arcsin(np.clip(z, -1, 1)))


def max_pixrad_deg(nside: int) -> float:
    """A safe upper bound on the centre-to-corner angle of any pixel (healpy's max_pixrad
    is ~0.82 x sqrt(pixel area) at nside 1 and smaller at higher nside)."""
    return float(np.degrees(np.sqrt(4 * np.pi / (12.0 * nside * nside)))) * 1.2


def _sep_deg(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    s = np.sin((d2 - d1) / 2) ** 2 + np.cos(d1) * np.cos(d2) * np.sin((r2 - r1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(s, 0, 1))))


def query_disc(nside: int, ra: float, dec: float, radius_deg: float) -> np.ndarray:
    """Nested pixels at `nside` (a power of 2) that may overlap the disc (conservative)."""
    pix = np.arange(12, dtype=np.int64)
    ns = 1
    while True:
        cra, cdec = pix2ang_nest(ns, pix)
        pix = pix[_sep_deg(ra, dec, cra, cdec) <= radius_deg + max_pixrad_deg(ns)]
        if ns == nside:
            return np.sort(pix)
        pix = (pix[:, None] * 4 + np.arange(4)[None, :]).ravel()
        ns *= 2


def query_disc_fine(ra: float, dec: float, radius_deg: float, nside_file: int = 32,
                    nside_fine: int = 4096) -> dict[int, np.ndarray]:
    """{file pixel: sorted fine (nside_fine) pixels inside it that may overlap the disc}."""
    fine = query_disc(nside_fine, ra, dec, radius_deg)
    shift = 2 * int(np.log2(nside_fine // nside_file))
    fp = fine >> shift
    return {int(p): fine[fp == p] for p in np.unique(fp)}

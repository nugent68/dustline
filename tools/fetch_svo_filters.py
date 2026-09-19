"""One-off: fetch missing filter curves from the SVO Filter Profile Service in the
same .dat format as the existing tools/filters files (header: ZeroPoint <Jy>, then
two columns wave_A T).

  python tools/fetch_svo_filters.py
"""

import os
import re
import sys

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "src", "dustline", "data", "filters")

FILTERS = {
    "PS1_g": "PAN-STARRS/PS1.g",
    "PS1_r": "PAN-STARRS/PS1.r",
    "PS1_i": "PAN-STARRS/PS1.i",
    "Cousins_I": "Generic/Cousins.I",
    "Johnson_V": "Generic/Johnson.V",
    "GALEX_FUV": "GALEX/GALEX.FUV",
    "GALEX_NUV": "GALEX/GALEX.NUV",
}

FPS = "http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID={}"
DATA = "http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id={}"


def meta(fid):
    xml = requests.get(FPS.format(fid), timeout=30).text
    zp = float(re.search(r'name="ZeroPoint"[^/]*?value="([\d.eE+-]+)"', xml).group(1))
    det = re.search(r'name="DetectorType"[^/]*?value="(\d)"', xml)
    leff = float(re.search(r'name="WavelengthEff"[^/]*?value="([\d.eE+-]+)"', xml).group(1))
    return zp, (det.group(1) if det else "?"), leff


def main():
    os.makedirs(OUT, exist_ok=True)
    for name, fid in FILTERS.items():
        path = os.path.join(OUT, f"{name}.dat")
        if os.path.exists(path):
            print(f"{name}: exists")
            continue
        zp, det, leff = meta(fid)
        body = requests.get(DATA.format(fid), timeout=30).text
        rows = [ln for ln in body.splitlines() if ln.strip() and not ln.startswith("#")]
        with open(path, "w") as f:
            f.write(f"# SVO FPS {fid}; DetectorType {det} (0 energy, 1 photon counter); "
                    f"ZeroPoint {zp} Jy (Pogson); lambda_eff {leff} A\n")
            f.write("\n".join(rows) + "\n")
        print(f"{name}: {fid} zp {zp:.1f} Jy, {len(rows)} points -> {path}")


if __name__ == "__main__":
    sys.exit(main())

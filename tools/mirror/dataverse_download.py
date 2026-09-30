"""Download every file of a Harvard Dataverse dataset (md5-verified, resumable)."""
import hashlib, json, os, subprocess, sys, urllib.request
doi, dest = sys.argv[1], sys.argv[2]
os.makedirs(dest, exist_ok=True)
meta = json.load(urllib.request.urlopen(urllib.request.Request(f"https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId={doi}", headers={"User-Agent": "dustline-mirror (curl-compatible)"}), timeout=120))
files = [f["dataFile"] for f in meta["data"]["latestVersion"]["files"]]
json.dump(files, open(os.path.join(dest, "_files.json"), "w"), indent=1)
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()
for f in sorted(files, key=lambda f: f["filename"]):
    p = os.path.join(dest, f["filename"])
    want = f["checksum"]["value"]
    for k in range(3):
        if os.path.exists(p) and os.path.getsize(p) == f["filesize"] and md5(p) == want:
            print("ok", f["filename"], flush=True)
            break
        subprocess.run(["curl", "-sSfL", "-A", "dustline-mirror", "--retry", "5", "-C", "-", "-o", p,
                        f"https://dataverse.harvard.edu/api/access/datafile/{f['id']}"], check=False)
        if os.path.exists(p) and md5(p) != want and os.path.getsize(p) >= f["filesize"]:
            os.remove(p)
    else:
        print("FAILED", f["filename"], flush=True)

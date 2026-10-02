# -*- coding: utf-8 -*-
"""Download TARGET-OS masked MAF files from GDC (parallel, robust)."""
import os, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

BASE = r"D:/projects/OS_matrix_dualKO/data/raw"
IDS  = os.path.join(BASE, "TARGET-OS_maf_ids.txt")
OUTD = os.path.join(BASE, "TARGET-OS_maf")
URL   = "https://api.gdc.cancer.gov/data/{}"
HDRS  = {"User-Agent": "Mozilla/5.0 (research script)"}

def valid(path):
    try:
        with open(path, "rb") as f:
            magic = f.read(2)
        return magic == b"\x1f\x8b"
    except OSError:
        return False

def fetch(fid):
    out = os.path.join(OUTD, fid + ".maf.gz")
    if os.path.exists(out) and valid(out):
        return (fid, "skip")
    for attempt in range(3):
        try:
            req = urllib.request.Request(URL.format(fid), headers=HDRS)
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            with open(out, "wb") as f:
                f.write(data)
            return (fid, "ok")
        except Exception as e:
            if attempt == 2:
                return (fid, f"fail:{type(e).__name__}")
            time.sleep(1.5)
    return (fid, "fail:unknown")

def main():
    with open(IDS) as f:
        ids = [x.strip() for x in f if x.strip()]
    print("total ids:", len(ids))
    ok = skip = fail = 0
    fails = []
    with ThreadPoolExecutor(max_workers=8) as ex:
        for fid, st in ex.map(fetch, ids):
            if st == "ok":
                ok += 1
            elif st == "skip":
                skip += 1
            else:
                fail += 1
                fails.append((fid, st))
            done = ok + skip + fail
            if done % 40 == 0:
                print(f"progress {done}/{len(ids)}", flush=True)
    print(f"RESULT ok={ok} skip={skip} fail={fail}")
    for fid, st in fails[:10]:
        print("  ", fid, st)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
organize_seeds.py - group AlphaFold Server seed folders into one folder per condition,
matching the "folder" column of conditions.csv.

    top2b_6mg2atp_nodna_rndm10a_sd1  \
    top2b_6mg2atp_nodna_rndm10a_sd2   >  data/TOP2B_randomA_noDNA/
    top2b_6mg2atp_nodna_rndm10a_sd3  /

Run from the folder that holds conditions.csv (e.g. C:\\AF3 or ~/AF3):
    python organize_seeds.py                     # shows what it would do, moves nothing
    python organize_seeds.py --move              # does it
    python organize_seeds.py --data "PATH"       # if the seed folders are not in .\\data

It searches up to three folder levels below --data, ignores upper/lower case and extra
endings on the names (e.g. "..._sd1 2", "..._sd1-1"), and unzips .zip downloads with --move.
"""
import argparse
import re
import shutil
import zipfile
from pathlib import Path

DNA = [(re.compile(r"nodna", re.I), "noDNA"),
       (re.compile(r"75bps?1?_?50bps?2?|50.?75", re.I), "50_75bp")]
SET = re.compile(r"rndm10([abc])", re.I)
SEED = re.compile(r"sd(\d+)", re.I)


def classify(name):
    s = SET.search(name)
    if not s or not SEED.search(name):
        return None
    for pat, tag in DNA:
        if pat.search(name):
            return f"TOP2B_random{s.group(1).upper()}_{tag}"
    return None


def candidates(root, depth=3):
    out, level = [], [root]
    for _ in range(depth):
        nxt = []
        for d in level:
            for p in sorted(d.iterdir()):
                if p.name.startswith("TOP2B_random"):      # already-grouped condition folders
                    continue
                if p.is_dir():
                    (out if classify(p.name) else nxt).append(p)
                elif p.suffix.lower() == ".zip" and classify(p.stem):
                    out.append(p)
        level = nxt
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--move", action="store_true", help="actually move (and unzip) the seed folders")
    a = ap.parse_args()
    data = Path(a.data).expanduser()
    if not data.is_dir():
        raise SystemExit(f"Folder not found: {data.resolve()}\n"
                         "Run this from the folder that contains 'data', or give --data \"PATH\".")
    found = candidates(data)
    if not found:
        print(f"No seed folders with 'rndm10a/b/c' and 'sd1/2/3' in their names under {data.resolve()}")
        print("What is there (first two levels):")
        for p in sorted(data.iterdir())[:40]:
            print(f"    {p.name}{'/' if p.is_dir() else ''}")
            if p.is_dir():
                for q in sorted(p.iterdir())[:8]:
                    print(f"        {q.name}{'/' if q.is_dir() else ''}")
        raise SystemExit("\nSend me this list and I will adjust the patterns.")
    plan = {}
    for p in found:
        plan.setdefault(classify(p.stem if p.is_file() else p.name), []).append(p)
    for dest in sorted(plan):
        items = plan[dest]
        flag = "" if len(items) == 3 else f"   <-- {len(items)} seed folders, expected 3"
        print(f"{dest}{flag}")
        for p in items:
            print(f"    {p.relative_to(data)}{'  (zip: will be unzipped)' if p.is_file() else ''}")
            if a.move:
                target = data / dest
                target.mkdir(exist_ok=True)
                if p.is_file():
                    with zipfile.ZipFile(p) as z:
                        z.extractall(target / p.stem)
                else:
                    shutil.move(str(p), str(target / p.name))
    print("\nDone." if a.move else "\nNothing moved yet. If this looks right, run again with --move.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
run_all_pae.py — run the PAE analysis for every condition listed in conditions.csv,
then make the comparison (difference) maps and one summary table across conditions.

Needs af3_pae_analysis.py in the same folder (it reuses its functions).

Folder layout (Windows example):
    C:\\AF3\\
        af3_pae_analysis.py
        run_all_pae.py
        conditions.csv
        data\\
            TOP2B_noPTM_noDNA\\      <- unzipped AF3 output for all 3 seeds of this condition
            TOP2B_noPTM_50_75bp\\
            ...
Run:
    cd C:\\AF3
    python run_all_pae.py

Options:
    --only NAME[,NAME...]   run just these conditions (comparisons still use saved results)
    --force                 redo conditions that already have results (default: skip them)
    --per-model             also save a PDF map for each of the 15 models (slow, ~100 MB/condition)
    --no-compare            skip the difference maps
    --tables-only           just rebuild the summary + per-seed tables from saved results (fast)
    --data DIR / --results DIR / --conditions FILE   change the default locations
    Font sizes (points, at the figure sizes below):
    --tick-size 12 --label-size 14 --title-size 13 --cbar-size 13   (old script: 7 / 10 / 10 / 10)
    --font-scale 1.25        make every font 25% bigger (or 0.8 smaller)
    --zoom-scale 0.8         fonts on the small CTD zoom maps, relative to the full maps

What it makes, per condition (results\\<name>\\):
    full_map\\PAE_mean.svg/.pdf/.png, PAE_sd.*      whole complex, 7 x 6 in
    CTD_zoom\\PAE_mean.*, PAE_sd.*                  CTDs of both chains, 3.5 x 3 in, no title
    CTD_DNA_zoom\\PAE_mean.*                        CTDs + DNA (DNA conditions only)
    per_model_metrics.csv                          pTM/ipTM/ranking score + region PAE, per model
    PAE_mean.npy, residue_index.csv                saved for the comparison step
Comparisons (results\\comparisons\\<A>__minus__<B>\\):
    PAE_difference_full.* and PAE_difference_CTD.*  (blue = more confident in A)
    made automatically for: every PTM condition vs the matching no-PTM condition,
    every DNA condition vs the matching no-DNA condition, plus anything in compare_to.
Summary:
    results\\all_conditions_summary.csv   mean and SD of every metric, one row per condition
    results\\<name>\\per_seed_metrics.csv   mean and SD of every metric for each seed (5 models each)
    results\\all_conditions_per_seed.csv  the same, all conditions in one table (one row per seed)
    results\\all_conditions_seed_spread.csv  per condition: mean of the 3 seed means, SD of the
                                         seed means (n = 3) and the range (lowest-highest seed mean)
    The per-seed tables are rebuilt on every run from the saved per_model_metrics.csv files,
    so existing results do not need --force.
    results\\run_log.txt                  what ran, what was skipped, any errors
"""
import argparse
import csv
import json
import re
import sys
import time
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import af3_pae_analysis as pae  # noqa: E402

METRICS = ["ptm", "iptm", "ranking_score", "fraction_disordered",
           "PAE_CTD_to_CTD", "PAE_CTD_aligned_core_scored", "PAE_core_aligned_CTD_scored",
           "PAE_CTD_aligned_DNA_scored", "PAE_DNA_aligned_CTD_scored", "mean_PAE_all"]


def log(msg, fh):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


# ----------------------------------------------------------------------------- chain typing
def chain_types(conf_path):
    """Classify chains as protein / dna / ligand from one AF3 confidence file.

    Ligands/ions: < 5 residues. Nucleic acids have ~20 heavy atoms per residue,
    proteins ~8, so atoms-per-residue > 14 => DNA/RNA.
    """
    d = json.loads(Path(conf_path).read_text())
    res = {}
    for c, r in zip(d["token_chain_ids"], d["token_res_ids"]):
        res.setdefault(c, set()).add(r)
    atoms = {}
    for c in d.get("atom_chain_ids", []):
        atoms[c] = atoms.get(c, 0) + 1
    types = {}
    for c, rs in res.items():
        n = len(rs)
        if n < 5:
            types[c] = "ligand"
        elif atoms and atoms.get(c, 0) / n > 14:
            types[c] = "dna"
        else:
            types[c] = "protein"
    return types, {c: (min(rs), max(rs)) for c, rs in res.items()}


# ----------------------------------------------------------------------------- one condition
def run_condition(cond, data_dir, out_root, args, fh):
    name = cond["name"]
    folder = data_dir / cond["folder"]
    out = out_root / name
    if not folder.exists():
        log(f"SKIP {name}: folder not found ({folder})", fh)
        return None
    if (out / "PAE_mean.npy").exists() and not args.force:
        log(f"SKIP {name}: results already exist (use --force to redo)", fh)
        return out
    out.mkdir(parents=True, exist_ok=True)

    try:
        models = pae.find_models([folder])
    except SystemExit as e:
        log(f"ERROR {name}: {e}", fh)
        return None
    log(f"{name}: {len(models)} models found", fh)
    if len(models) != 15:
        log(f"  WARNING {name}: expected 15 models (3 seeds x 5), found {len(models)}", fh)

    types, ranges = chain_types(models[0]["conf_path"])
    prot = [c for c, t in types.items() if t == "protein"]
    dna = [c for c, t in types.items() if t == "dna"]
    if cond.get("dna_chains"):
        dna = [c.strip() for c in cond["dna_chains"].split(",") if c.strip()]
        prot = [c for c in prot if c not in dna]
    lo, hi = (int(x) for x in cond["ctd"].split("-"))
    ctd_regions = [(c, lo, hi) for c in prot]
    log(f"  protein chains {prot}, DNA chains {dna or 'none'}, CTD {lo}-{hi}", fh)

    keys0, rows = None, []
    total = sumsq = None
    for m in models:
        keys, mat = pae.load_residue_pae(m["conf_path"])
        mat = mat.astype(np.float64)
        if keys0 is None:
            keys0 = keys
            total = np.zeros_like(mat)
            sumsq = np.zeros_like(mat)
            is_ctd = pae.mask_for(keys, ctd_regions)
            is_dna = np.array([c in dna for c, _ in keys])
            is_prot = np.array([c in prot for c, _ in keys])
            is_core = is_prot & ~is_ctd
        elif keys != keys0:
            log(f"ERROR {name}: {m['name']} has a different residue layout — "
                "this folder mixes inputs from different conditions", fh)
            return None
        total += mat
        sumsq += mat * mat
        row = {"model": m["name"], **pae.load_summary(m["summary_path"])}
        row["mean_PAE_all"] = float(mat.mean())
        row["PAE_CTD_to_CTD"] = pae.block_mean(mat, is_ctd, is_ctd)
        row["PAE_CTD_aligned_core_scored"] = pae.block_mean(mat, is_ctd, is_core)
        row["PAE_core_aligned_CTD_scored"] = pae.block_mean(mat, is_core, is_ctd)
        if dna:
            row["PAE_CTD_aligned_DNA_scored"] = pae.block_mean(mat, is_ctd, is_dna)
            row["PAE_DNA_aligned_CTD_scored"] = pae.block_mean(mat, is_dna, is_ctd)
        rows.append(row)
        if args.per_model:
            pae.FORMATS, pae.FIGSIZE, pae.TITLES = ["pdf"], (7, 6), True
            (out / "per_model").mkdir(exist_ok=True)
            pae.FONT_SCALE = args.font_scale
            pae.plot_map(mat, keys, f"{name}\n{m['name']}", out / "per_model" / f"PAE_{m['name']}.pdf")
        del mat

    n = len(rows)
    mean = total / n
    sd = np.sqrt(np.maximum(sumsq / n - mean ** 2, 0) * n / (n - 1)) if n > 1 else np.zeros_like(mean)
    np.save(out / "PAE_mean.npy", mean.astype(np.float32))
    with open(out / "residue_index.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["index", "chain", "residue"])
        for i, (c, r) in enumerate(keys0):
            w.writerow([i, c, r])

    # figures
    sd_max = max(1.0, float(np.percentile(sd, 99)))
    views = [("full_map", None, (7, 6), True)]
    views.append(("CTD_zoom", ctd_regions, (3.5, 3.0), False))
    if dna:
        views.append(("CTD_DNA_zoom",
                      ctd_regions + [(c, *ranges[c]) for c in dna], (3.5, 3.0), False))
    for vname, zoom, size, titles in views:
        d = out / vname
        d.mkdir(exist_ok=True)
        pae.FORMATS, pae.FIGSIZE, pae.TITLES = args.formats, size, titles
        pae.FONT_SCALE = args.font_scale * (1.0 if size[0] >= 6 else args.zoom_scale)
        pae.plot_map(mean, keys0, f"{name}: mean PAE (n = {n} models)", d / "PAE_mean.png", zoom=zoom)
        pae.plot_map(sd, keys0, f"{name}: SD of PAE (n = {n})", d / "PAE_sd.png",
                     vmax=sd_max, cmap="magma", cbar="SD (Å)", zoom=zoom)

    # per-model table + mean/SD rows
    fields = ["model"] + [k for k in METRICS if any(k in r for r in rows)]
    summary = {}
    with open(out / "per_model_metrics.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
        mrow, srow = {"model": "MEAN"}, {"model": "SD"}
        for k in fields[1:]:
            vals = [r[k] for r in rows if isinstance(r.get(k), (int, float))
                    and not isinstance(r.get(k), bool)]
            if len(vals) > 1:
                mrow[k] = round(float(np.mean(vals)), 3)
                srow[k] = round(float(np.std(vals, ddof=1)), 3)
                summary[k] = (mrow[k], srow[k])
        w.writerow(mrow)
        w.writerow(srow)
    (out / "chains.json").write_text(json.dumps({"protein": prot, "dna": dna,
                                                 "ctd": [lo, hi], "n_models": n}))
    log(f"  done {name}", fh)
    return out


# ----------------------------------------------------------------------------- comparisons
def load_saved(out):
    keys = [(r["chain"], int(r["residue"])) for r in csv.DictReader(open(out / "residue_index.csv"))]
    info = json.loads((out / "chains.json").read_text())
    return keys, np.load(out / "PAE_mean.npy"), info


def compare(a_name, b_name, out_root, args, fh):
    a_out, b_out = out_root / a_name, out_root / b_name
    if not (a_out / "PAE_mean.npy").exists() or not (b_out / "PAE_mean.npy").exists():
        log(f"SKIP comparison {a_name} - {b_name}: one of them has no results yet", fh)
        return
    dest = out_root / "comparisons" / f"{a_name}__minus__{b_name}"
    if (dest / "PAE_difference_CTD.png").exists() and not args.force:
        return
    ka, ma, ia = load_saved(a_out)
    kb, mb, ib = load_saved(b_out)
    lookup_b = {k: i for i, k in enumerate(kb)}
    shared = [k for k in ka if k in lookup_b]
    prot_shared = [k for k in shared if k[0] in ia["protein"]]
    if len(prot_shared) < 100:
        log(f"SKIP comparison {a_name} - {b_name}: proteins differ (no shared residue numbering)", fh)
        return
    lookup_a = {k: i for i, k in enumerate(ka)}
    ia_idx = [lookup_a[k] for k in shared]
    ib_idx = [lookup_b[k] for k in shared]
    diff = ma[np.ix_(ia_idx, ia_idx)] - mb[np.ix_(ib_idx, ib_idx)]
    dest.mkdir(parents=True, exist_ok=True)
    lo, hi = ia["ctd"]
    ctd = [(c, lo, hi) for c in ia["protein"]]
    lim = max(1.0, float(np.percentile(np.abs(diff), 99)))
    title = f"Mean PAE: {a_name} minus {b_name}\n(blue = more confident in {a_name})"
    for fname, zoom, size in (("PAE_difference_full.png", None, (7, 6)),
                              ("PAE_difference_CTD.png", ctd, (3.5, 3.0))):
        pae.FORMATS, pae.FIGSIZE, pae.TITLES = args.formats, size, True
        pae.FONT_SCALE = args.font_scale * (1.0 if size[0] >= 6 else args.zoom_scale)
        pae.plot_map(diff, shared, title, dest / fname, vmin=-lim, vmax=lim,
                     cmap="RdBu_r", cbar="ΔPAE (Å)", zoom=zoom)
    log(f"  comparison {a_name} - {b_name} done", fh)


def auto_pairs(conds):
    """PTM vs matching no-PTM; DNA vs matching no-DNA; plus explicit compare_to."""
    by_key = {(c["protein"], c["ptm"], c["dna"]): c["name"] for c in conds}
    pairs = []
    for c in conds:
        if c["ptm"].lower() != "none":
            b = by_key.get((c["protein"], "none", c["dna"]))
            if b:
                pairs.append((c["name"], b))
        if c["dna"].lower() != "none":
            b = by_key.get((c["protein"], c["ptm"], "none"))
            if b:
                pairs.append((c["name"], b))
        for b in (c.get("compare_to") or "").split(";"):
            if b.strip():
                pairs.append((c["name"], b.strip()))
    return list(dict.fromkeys(pairs))



# ----------------------------------------------------------------------------- per-seed tables
SEED_RE = re.compile(r"(?:^|[_\-/])(?:sd|seed)[_\-]?(\d+)(?:$|[_\-/])", re.IGNORECASE)


def seed_of(model_name):
    """Seed for one model, from its name ('<job folder>/seed<N>_model<k>').
    Uses the seed read from the AF3 job_request file when there was one; otherwise looks
    for 'sd1', 'seed1', 'seed-1', 'seed_1' in the job folder name; otherwise uses the job
    folder name itself (each AlphaFold Server job is one seed)."""
    folder, _, tail = model_name.rpartition("/")
    m = re.match(r"seed([^_]+)_model", tail)
    if m and m.group(1) != "?":
        return m.group(1)
    m = SEED_RE.search(folder + "_") or SEED_RE.search(model_name)
    if m:
        return m.group(1)
    return folder or "?"


def _num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x else x


def seed_tables(conds, out_root, fh):
    long_rows, spread_rows = [], []
    for c in conds:
        f = out_root / c["name"] / "per_model_metrics.csv"
        if not f.exists():
            continue
        recs = [r for r in csv.DictReader(open(f)) if r["model"] not in ("MEAN", "SD")]
        metrics = [k for k in METRICS if any(_num(r.get(k)) is not None for r in recs)]
        groups = {}
        for r in recs:
            groups.setdefault(seed_of(r["model"]), []).append(r)
        seeds = sorted(groups, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else 0, x))
        if len(seeds) != 3:
            log(f"  NOTE {c['name']}: found {len(seeds)} seed group(s) {seeds} (expected 3)", fh)
        per_seed = []
        for sd_ in seeds:
            row = {"condition": c["name"], "protein": c["protein"], "ptm": c["ptm"], "dna": c["dna"],
                   "seed": sd_, "n_models": len(groups[sd_])}
            for k in metrics:
                vals = [x for x in (_num(r.get(k)) for r in groups[sd_]) if x is not None]
                if vals:
                    row[f"{k}_mean"] = round(float(np.mean(vals)), 3)
                    row[f"{k}_sd"] = round(float(np.std(vals, ddof=1)), 3) if len(vals) > 1 else ""
            per_seed.append(row)
        long_rows += per_seed
        fields = list(dict.fromkeys(k for r in per_seed for k in r))
        with open(out_root / c["name"] / "per_seed_metrics.csv", "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=fields)
            w.writeheader()
            w.writerows(per_seed)
        sp = {"condition": c["name"], "protein": c["protein"], "ptm": c["ptm"], "dna": c["dna"],
              "n_seeds": len(per_seed), "models_per_seed": "/".join(str(r["n_models"]) for r in per_seed)}
        for k in metrics:
            means = [r[f"{k}_mean"] for r in per_seed if f"{k}_mean" in r]
            if not means:
                continue
            sp[f"{k}_mean_of_seeds"] = round(float(np.mean(means)), 3)
            sp[f"{k}_sd_of_seed_means"] = round(float(np.std(means, ddof=1)), 3) if len(means) > 1 else ""
            sp[f"{k}_seed_range"] = f"{min(means):.3f}-{max(means):.3f}"
        spread_rows.append(sp)
    for name, rows in (("all_conditions_per_seed.csv", long_rows),
                       ("all_conditions_seed_spread.csv", spread_rows)):
        if rows:
            fields = list(dict.fromkeys(k for r in rows for k in r))
            with open(out_root / name, "w", newline="") as fo:
                w = csv.DictWriter(fo, fieldnames=fields)
                w.writeheader()
                w.writerows(rows)
    log(f"=== per-seed tables written for {len(spread_rows)} conditions: "
        "all_conditions_per_seed.csv, all_conditions_seed_spread.csv", fh)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conditions", default=str(HERE / "conditions.csv"))
    ap.add_argument("--data", default=str(HERE / "data"))
    ap.add_argument("--results", default=str(HERE / "results"))
    ap.add_argument("--only")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--per-model", action="store_true")
    ap.add_argument("--no-compare", action="store_true")
    ap.add_argument("--tables-only", action="store_true",
                    help="only rebuild the summary and per-seed tables from saved results")
    ap.add_argument("--format", default="svg,pdf,png")
    ap.add_argument("--tick-size", type=float, default=pae.FONT["tick"])
    ap.add_argument("--label-size", type=float, default=pae.FONT["label"])
    ap.add_argument("--title-size", type=float, default=pae.FONT["title"])
    ap.add_argument("--cbar-size", type=float, default=pae.FONT["cbar"])
    ap.add_argument("--font-scale", type=float, default=1.0)
    ap.add_argument("--zoom-scale", type=float, default=0.8)
    args = ap.parse_args()
    args.formats = [f.strip() for f in args.format.split(",")]
    pae.FONT.update(tick=args.tick_size, label=args.label_size, title=args.title_size,
                    cbar=args.cbar_size, cbar_tick=args.tick_size)

    conds = [r for r in csv.DictReader(open(args.conditions, newline="", encoding="utf-8-sig"))
             if r.get("name") and not r["name"].startswith("#")]
    for c in conds:
        for k in ("protein", "ptm", "dna"):
            c[k] = (c.get(k) or "none").strip()
    out_root = Path(args.results)
    out_root.mkdir(parents=True, exist_ok=True)
    data_dir = Path(args.data)
    only = set(args.only.split(",")) if args.only else None

    with open(out_root / "run_log.txt", "a", encoding="utf-8") as fh:
        log(f"=== run started: {len(conds)} conditions in {args.conditions}", fh)
        for i, c in enumerate(conds, 1):
            if args.tables_only:
                break
            if only and c["name"] not in only:
                continue
            log(f"[{i}/{len(conds)}] {c['name']}", fh)
            try:
                run_condition(c, data_dir, out_root, args, fh)
            except Exception:
                log(f"ERROR {c['name']}:\n{traceback.format_exc()}", fh)

        if not args.no_compare and not args.tables_only:
            log("=== comparisons", fh)
            for a, b in auto_pairs(conds):
                try:
                    compare(a, b, out_root, args, fh)
                except Exception:
                    log(f"ERROR comparison {a} - {b}:\n{traceback.format_exc()}", fh)

        # summary across all conditions
        rows = []
        for c in conds:
            f = out_root / c["name"] / "per_model_metrics.csv"
            if not f.exists():
                continue
            recs = list(csv.DictReader(open(f)))
            mean = next((r for r in recs if r["model"] == "MEAN"), {})
            sd = next((r for r in recs if r["model"] == "SD"), {})
            row = {"condition": c["name"], "protein": c["protein"], "ptm": c["ptm"], "dna": c["dna"],
                   "n_models": len(recs) - 2}
            for k in METRICS:
                if mean.get(k):
                    row[f"{k}_mean"] = mean[k]
                    row[f"{k}_sd"] = sd.get(k, "")
            rows.append(row)
        if rows:
            fields = list(dict.fromkeys(k for r in rows for k in r))
            with open(out_root / "all_conditions_summary.csv", "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields)
                w.writeheader()
                w.writerows(rows)
            log(f"=== summary written for {len(rows)} conditions: {out_root / 'all_conditions_summary.csv'}", fh)
        seed_tables(conds, out_root, fh)
        log("=== run finished", fh)


if __name__ == "__main__":
    main()

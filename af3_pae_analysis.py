#!/usr/bin/env python3
"""
af3_pae_analysis.py — PAE maps and confidence summaries for AlphaFold 3 runs.

Requirements: Python 3.9+, numpy, matplotlib.
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42


def find_models(folders):
    models = []
    for folder in folders:
        folder = Path(folder)
        if not folder.exists():
            sys.exit(f"Folder not found: {folder}")
        for p in sorted(folder.rglob("*_full_data_*.json")):
            m = re.search(r"(.*)_full_data_(\d+)(?:_\d+)?\.json$", p.name)
            if not m:
                continue
            cands = sorted(p.parent.glob(f"{m.group(1)}_summary_confidences_{m.group(2)}*.json"))
            summ = cands[0] if cands else p.with_name("__none__")
            seed = server_seed(p.parent)
            models.append(dict(name=f"{p.parent.name}/seed{seed}_model{m.group(2)}",
                               conf_path=p, summary_path=summ if summ.exists() else None))
        for p in sorted(folder.rglob("*confidences.json")):
            if "summary" in p.name or "_full_data_" in p.name:
                continue
            if not re.search(r"seed-\d+_sample-\d+", str(p.parent.name)):
                continue
            summ = p.with_name(p.name.replace("confidences.json", "summary_confidences.json"))
            models.append(dict(name=f"{p.parent.parent.name}/{p.parent.name}",
                               conf_path=p, summary_path=summ if summ.exists() else None))
    if not models:
        sys.exit("No AF3 confidence files found. Point the script at the folder(s) containing the unzipped AF3 output.")
    return models


def server_seed(job_dir):
    for req in job_dir.glob("*job_request*.json"):
        try:
            data = json.loads(req.read_text())
            data = data[0] if isinstance(data, list) else data
            seeds = data.get("modelSeeds") or []
            if seeds:
                return seeds[0]
        except Exception:
            pass
    return "?"


def load_residue_pae(conf_path):
    d = json.loads(Path(conf_path).read_text())
    pae = np.asarray(d["pae"], dtype=np.float32)
    chains = d["token_chain_ids"]
    resids = d["token_res_ids"]
    keys, index = [], {}
    tok2res = np.empty(len(chains), dtype=int)
    for t, (c, r) in enumerate(zip(chains, resids)):
        k = (c, int(r))
        if k not in index:
            index[k] = len(keys)
            keys.append(k)
        tok2res[t] = index[k]
    n_res, n_tok = len(keys), len(chains)
    M = np.zeros((n_res, n_tok), dtype=np.float32)
    M[tok2res, np.arange(n_tok)] = 1.0
    M /= M.sum(axis=1, keepdims=True)
    return keys, M @ pae @ M.T


def load_summary(path):
    if path is None:
        return {}
    d = json.loads(Path(path).read_text())
    return {k: d.get(k) for k in ("ptm", "iptm", "ranking_score", "fraction_disordered", "has_clash")}


def parse_regions(spec):
    out = []
    if not spec:
        return out
    for part in spec.split(","):
        chain, rng = part.strip().split(":")
        lo, hi = rng.split("-")
        out.append((chain.strip(), int(lo), int(hi)))
    return out


def mask_for(keys, regions):
    return np.array([any(c == rc and lo <= r <= hi for rc, lo, hi in regions) for c, r in keys])


def block_mean(pae, rows, cols):
    if rows.sum() == 0 or cols.sum() == 0:
        return float("nan")
    return float(pae[np.ix_(rows, cols)].mean())


FORMATS = ["png"]
FIGSIZE = (7, 6)
TITLES = True
FONT = dict(tick=12, label=14, title=13, cbar=13, cbar_tick=12)
FONT_SCALE = 1.0
N_TICKS = 6


def _pt(key):
    return FONT[key] * FONT_SCALE


def chain_bounds(keys):
    bounds, start = [], 0
    for i in range(1, len(keys) + 1):
        if i == len(keys) or keys[i][0] != keys[i - 1][0]:
            bounds.append(i)
            start = i
    return bounds[:-1]


def residue_ticks(keys, target=N_TICKS):
    n = len(keys)
    ticks = []
    for step in (10, 20, 25, 50, 100, 200, 250, 500, 1000):
        cand = [i for i, (c, r) in enumerate(keys) if r % step == 0]
        if 2 <= len(cand) <= target + 2:
            ticks = cand
            break
    if len(ticks) < 2:
        ticks = list(range(0, n, max(1, n // target)))
    min_gap = 0.06 * n
    kept = []
    for t in ticks:
        if not kept or t - kept[-1] >= min_gap:
            kept.append(t)
    return kept


def plot_map(mat, keys, title, path, vmax=31.75, cmap="Greens_r", cbar="Expected position error (Å)", vmin=0, zoom=None):
    if zoom:
        sel = mask_for(keys, zoom)
        mat = mat[np.ix_(sel, sel)]
        keys = [k for k, s in zip(keys, sel) if s]
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=300)
    im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    for x in chain_bounds(keys):
        ax.axhline(x - 0.5, color="k", lw=0.8)
        ax.axvline(x - 0.5, color="k", lw=0.8)
    ticks = residue_ticks(keys)
    names = [f"{keys[t][0]}{keys[t][1]}" for t in ticks]
    ax.set_xticks(ticks); ax.set_xticklabels(names, rotation=90, fontsize=_pt("tick"))
    ax.set_yticks(ticks); ax.set_yticklabels(names, fontsize=_pt("tick"))
    ax.set_xlabel("Scored residue", fontsize=_pt("label"))
    ax.set_ylabel("Aligned residue", fontsize=_pt("label"))
    if TITLES:
        ax.set_title(title, fontsize=_pt("title"))
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label(cbar, fontsize=_pt("cbar"))
    cb.ax.tick_params(labelsize=_pt("cbar_tick"))
    fig.tight_layout()
    path = Path(path)
    for fmt in FORMATS:
        fig.savefig(path.with_suffix("." + fmt.lstrip(".")))
    plt.close(fig)


def analyse_condition(folders, label, args, outdir):
    models = find_models(folders)
    print(f"[{label}] {len(models)} models found")
    keys0, maps, rows = None, [], []
    ctd = parse_regions(args.ctd)
    dna_chains = set(c.strip() for c in args.dna_chains.split(",")) if args.dna_chains else set()
    for m in models:
        keys, pae = load_residue_pae(m["conf_path"])
        if keys0 is None:
            keys0 = keys
        elif keys != keys0:
            sys.exit(f"Residue layout differs in {m['name']} — all models in one condition must come from the same input.")
        maps.append(pae)
        row = {"model": m["name"], **load_summary(m["summary_path"])}
        row["mean_PAE_all"] = float(pae.mean())
        if ctd:
            is_ctd = mask_for(keys, ctd)
            is_dna = np.array([c in dna_chains for c, _ in keys])
            is_core = ~is_ctd & ~is_dna
            row["PAE_CTD_to_CTD"] = block_mean(pae, is_ctd, is_ctd)
            row["PAE_CTD_aligned_core_scored"] = block_mean(pae, is_ctd, is_core)
            row["PAE_core_aligned_CTD_scored"] = block_mean(pae, is_core, is_ctd)
            if dna_chains:
                row["PAE_CTD_aligned_DNA_scored"] = block_mean(pae, is_ctd, is_dna)
                row["PAE_DNA_aligned_CTD_scored"] = block_mean(pae, is_dna, is_ctd)
        rows.append(row)
        if args.per_model:
            safe = re.sub(r"[^\w.-]+", "_", m["name"])
            plot_map(pae, keys, f"{label}\n{m['name']}", outdir / f"PAE_{safe}.png", zoom=parse_regions(args.zoom))
    stack = np.stack(maps)
    mean = stack.mean(0)
    sd = stack.std(0, ddof=1) if len(maps) > 1 else np.zeros_like(stack[0])
    zoom = parse_regions(args.zoom)
    plot_map(mean, keys0, f"{label}: mean PAE (n = {len(maps)} models)", outdir / "PAE_mean.png", zoom=zoom)
    plot_map(sd, keys0, f"{label}: SD of PAE across models (n = {len(maps)})", outdir / "PAE_sd.png",
             vmax=max(1.0, float(np.percentile(sd, 99))), cmap="magma", cbar="SD (Å)", zoom=zoom)
    np.save(outdir / "PAE_mean.npy", mean)
    np.save(outdir / "PAE_sd.npy", sd)
    with open(outdir / "residue_index.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["index", "chain", "residue"])
        for i, (c, r) in enumerate(keys0): w.writerow([i, c, r])
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(outdir / "per_model_metrics.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
        for stat, fn in (("MEAN", np.nanmean), ("SD", lambda a: np.nanstd(a, ddof=1))):
            agg = {"model": stat}
            for k in fields[1:]:
                vals = [r.get(k) for r in rows if isinstance(r.get(k), (int, float)) and not isinstance(r.get(k), bool)]
                agg[k] = round(float(fn(vals)), 3) if len(vals) > 1 else ""
            w.writerow(agg)
    print(f"[{label}] wrote {outdir}")
    return keys0, mean


def add_figure_args(ap):
    g = ap.add_argument_group("figure options")
    g.add_argument("--format", default=None)
    g.add_argument("--size")
    g.add_argument("--no-titles", action="store_true")
    g.add_argument("--tick-size", type=float, default=FONT["tick"])
    g.add_argument("--label-size", type=float, default=FONT["label"])
    g.add_argument("--title-size", type=float, default=FONT["title"])
    g.add_argument("--cbar-size", type=float, default=FONT["cbar"])
    g.add_argument("--font-scale", type=float, default=1.0)


def apply_figure_args(args):
    global FORMATS, FIGSIZE, TITLES, FONT_SCALE
    if args.format: FORMATS = [f.strip() for f in args.format.split(",") if f.strip()]
    if args.size:
        w, h = (float(x) for x in args.size.lower().split("x")); FIGSIZE = (w, h)
    if getattr(args, "no_titles", False): TITLES = False
    FONT.update(tick=args.tick_size, label=args.label_size, title=args.title_size,
                cbar=args.cbar_size, cbar_tick=args.tick_size)
    FONT_SCALE = args.font_scale


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folders", nargs="+")
    ap.add_argument("--label", default="condition")
    ap.add_argument("--out", default="pae_results")
    ap.add_argument("--ctd")
    ap.add_argument("--dna-chains")
    ap.add_argument("--zoom")
    ap.add_argument("--per-model", action="store_true")
    ap.add_argument("--compare", nargs="+")
    ap.add_argument("--compare-label", default="comparison")
    add_figure_args(ap)
    args = ap.parse_args(); apply_figure_args(args)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    keys_a, mean_a = analyse_condition(args.folders, args.label, args, out)
    if args.compare:
        out_b = out / "comparison_condition"; out_b.mkdir(exist_ok=True)
        keys_b, mean_b = analyse_condition(args.compare, args.compare_label, args, out_b)
        shared = [k for k in keys_a if k in set(keys_b)]
        ia = [keys_a.index(k) for k in shared]; ib_lookup = {k: i for i, k in enumerate(keys_b)}
        ib = [ib_lookup[k] for k in shared]
        diff = mean_a[np.ix_(ia, ia)] - mean_b[np.ix_(ib, ib)]
        lim = max(1.0, float(np.percentile(np.abs(diff), 99)))
        plot_map(diff, shared, f"Mean PAE: {args.label} minus {args.compare_label}\n(blue = lower error / more confident)",
                 out / "PAE_difference.png", vmin=-lim, vmax=lim, cmap="RdBu_r", cbar="ΔPAE (Å)", zoom=parse_regions(args.zoom))
        np.save(out / "PAE_difference.npy", diff)


if __name__ == "__main__":
    main()

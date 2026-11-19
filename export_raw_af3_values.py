#!/usr/bin/env python3
"""
export_raw_af3_values.py - every confidence value AlphaFold 3 reported, model by model, in one
Excel workbook (one tab per condition), read straight from the AF3 output files.

Nothing is recalculated except where stated: each number in the per-model rows is copied from
the AF3 JSON file named in the "source_file" column. Means and SDs are Excel formulas over
those rows, so anyone can check them in the workbook.

Per model (from fold_<job>_summary_confidences_<n>.json):
    ptm, iptm, ranking_score, fraction_disordered, has_clash, num_recycles,
    chain_ptm and chain_iptm for every chain, and the chain-pair ipTM / min PAE for
    protein chain A with every other chain
Optional, per model (from fold_<job>_full_data_<n>.json, slower: large files):
    mean pLDDT of all atoms and of each polymer chain (simple average of AF3's atom_plddts),
    mean of the full token-level PAE matrix (simple average of AF3's pae)

Tabs:
    README          what every column is and which file/key it comes from
    Index           one row per condition: folder, models found, seeds found
    <condition>     15 rows (3 seeds x 5 models), then MEAN/SD per seed, over all 15 models,
                    and the mean and SD of the three seed means (all Excel formulas)
    All models      every model of every condition in one long table

Run from the folder with conditions.csv (same layout as run_all_pae.py):
    python export_raw_af3_values.py                    # summary_confidences only (fast)
    python export_raw_af3_values.py --full-data        # also pLDDT / PAE means from full_data
Options: --conditions conditions.csv  --data data  --out results\\AF3_raw_values_by_model.xlsx
Requires: python -m pip install openpyxl
"""
import argparse
import csv
import json
import re
import string
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCALARS = ["ptm", "iptm", "ranking_score", "fraction_disordered", "has_clash", "num_recycles"]
SEED_IN_NAME = re.compile(r"(?:^|[_\-\s])(?:sd|seed)[_\-]?(\d+)(?=$|[_\-\s])", re.I)


# ----------------------------------------------------------------------------- finding models
def job_seed(job_dir):
    """Seed from the AF3 job_request file if present, else 'sd1'/'seed1' in the folder name."""
    for req in sorted(job_dir.glob("*job_request*.json")):
        try:
            d = json.loads(req.read_text())
            d = d[0] if isinstance(d, list) else d
            seeds = d.get("modelSeeds") or []
            if seeds:
                return str(seeds[0]), f"{req.name} (modelSeeds)"
        except Exception:
            pass
    m = SEED_IN_NAME.search(job_dir.name)
    if m:
        return m.group(1), "folder name"
    return job_dir.name, "folder (seed unknown)"


def find_models(folder):
    models = []
    for p in sorted(folder.rglob("*_summary_confidences_*.json")):
        m = re.search(r"(.*)_summary_confidences_(\d+)(?:_\d+)?\.json$", p.name)
        if not m:
            continue
        job, idx = m.group(1), int(m.group(2))
        full = sorted(p.parent.glob(f"{job}_full_data_{idx}*.json"))
        seed, seed_src = job_seed(p.parent)
        models.append(dict(job=job, job_dir=p.parent.name, model=idx, seed=seed, seed_source=seed_src,
                           summary=p, full=full[0] if full else None))
    # local AlphaFold 3 layout: seed-<s>_sample-<k>/summary_confidences.json
    for p in sorted(folder.rglob("*summary_confidences.json")):
        m = re.search(r"seed-(\d+)_sample-(\d+)", p.parent.name)
        if not m:
            continue
        full = p.with_name(p.name.replace("summary_confidences.json", "confidences.json"))
        models.append(dict(job=p.parent.parent.name, job_dir=p.parent.name, model=int(m.group(2)),
                           seed=m.group(1), seed_source="folder name", summary=p,
                           full=full if full.exists() else None))
    seen, out = set(), []
    for m in models:                                   # drop duplicate downloads of one model
        key = (m["job_dir"], m["model"])
        if key not in seen:
            seen.add(key)
            out.append(m)
    def seed_key(m):
        s = m["seed"]
        return (0, int(s), "") if s.isdigit() else (1, 0, s)
    return sorted(out, key=lambda m: (seed_key(m), m["model"]))


def chain_letters(n):
    letters = list(string.ascii_uppercase)
    return (letters + [a + b for a in letters for b in letters])[:n]


# ----------------------------------------------------------------------------- reading values
def read_model(m, full_data):
    s = json.loads(Path(m["summary"]).read_text())
    row = {"seed": m["seed"], "model": m["model"], "job": m["job"]}
    for k in SCALARS:
        v = s.get(k)
        row[k] = float(v) if isinstance(v, (int, float)) else v
    cp = s.get("chain_ptm") or []
    ci = s.get("chain_iptm") or []
    ids = chain_letters(max(len(cp), len(ci)))
    for c, v in zip(ids, cp):
        row[f"chain_ptm_{c}"] = v
    for c, v in zip(ids, ci):
        row[f"chain_iptm_{c}"] = v
    pi = s.get("chain_pair_iptm") or []
    pm = s.get("chain_pair_pae_min") or []
    for j, c in enumerate(ids[1:], start=1):
        if pi and len(pi[0]) > j:
            row[f"pair_iptm_A_{c}"] = pi[0][j]
        if pm and len(pm[0]) > j:
            row[f"pair_pae_min_A_{c}"] = pm[0][j]
    if full_data and m["full"]:
        f = json.loads(Path(m["full"]).read_text())
        pl, ach = f.get("atom_plddts"), f.get("atom_chain_ids")
        if pl:
            row["mean_pLDDT_all_atoms"] = round(sum(pl) / len(pl), 3)
            if ach:
                per = {}
                for c, v in zip(ach, pl):
                    per.setdefault(c, []).append(v)
                tok = f.get("token_chain_ids") or []
                ntok = {c: tok.count(c) for c in per}
                for c, vals in per.items():
                    if ntok.get(c, 0) > 1:               # polymer chains only (skip ions/ATP)
                        row[f"mean_pLDDT_chain_{c}"] = round(sum(vals) / len(vals), 3)
        pae = f.get("pae")
        if pae:
            tot = sum(sum(r) for r in pae)
            row["mean_PAE_all_tokens"] = round(tot / (len(pae) * len(pae[0])), 3)
        row["full_data_file"] = Path(m["full"]).name
    row["source_file"] = f"{m['job_dir']}/{Path(m['summary']).name}"
    row["seed_from"] = m["seed_source"]
    return row


# ----------------------------------------------------------------------------- workbook
def write_workbook(per_cond, conds_info, out, full_data):
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.utils import get_column_letter
    F = "Arial"
    bold, reg, ital = Font(name=F, bold=True), Font(name=F), Font(name=F, italic=True, size=9)
    head = PatternFill("solid", start_color="D9E1F2")
    stat = PatternFill("solid", start_color="F2F2F2")
    wb = Workbook()
    ws = wb.active
    ws.title = "README"
    lines = [
        ("AlphaFold 3 confidence values, one row per model, read directly from the AF3 output files", bold),
        ("", reg),
        ("Per-model rows: every value is copied unchanged from the AF3 JSON file named in 'source_file'.", reg),
        ("Rows below each table (grey) are Excel formulas over those rows: MEAN/SD for each seed (5 models),", reg),
        ("over all models of the condition, and the mean and SD of the seed means.", reg),
        ("", reg),
        ("Columns from fold_<job>_summary_confidences_<n>.json", bold),
        ("ptm, iptm, ranking_score, fraction_disordered, has_clash, num_recycles: AF3 keys of the same name.", reg),
        ("chain_ptm_X / chain_iptm_X: AF3 'chain_ptm' / 'chain_iptm' lists; X = chain ID in AF3 input order "
         "(A, B = protein chains; then ligands, ions and DNA in the order entered).", reg),
        ("pair_iptm_A_X / pair_pae_min_A_X: row A of AF3 'chain_pair_iptm' / 'chain_pair_pae_min'.", reg),
        ("seed: from the job_request file ('modelSeeds') when present, otherwise from 'sd1/sd2/sd3' in the "
         "folder name; 'seed_from' says which.", reg),
    ]
    if full_data:
        lines += [("", reg), ("Columns from fold_<job>_full_data_<n>.json", bold),
                  ("mean_pLDDT_all_atoms / mean_pLDDT_chain_X: arithmetic mean of AF3 'atom_plddts' "
                   "(all atoms, or atoms of chain X from 'atom_chain_ids'; polymer chains only).", reg),
                  ("mean_PAE_all_tokens: arithmetic mean of the full AF3 'pae' matrix (token level, no "
                   "collapsing of modified residues).", reg)]
    for i, (t, f) in enumerate(lines, 1):
        ws.cell(i, 1, t).font = f
    ws.column_dimensions["A"].width = 120

    wi = wb.create_sheet("Index")
    wi.append(["condition", "folder", "models found", "seeds found", "models per seed"])
    for c in wi[1]:
        c.font, c.fill = bold, head
    for info in conds_info:
        wi.append(list(info))
    for col in "ABCDE":
        wi.column_dimensions[col].width = 26

    all_rows = []
    for cname, rows in per_cond.items():
        sheet = wb.create_sheet(cname[:31])
        cols = list(dict.fromkeys(k for r in rows for k in r))
        front = ["seed", "model", "job"]
        tail = [k for k in ("source_file", "full_data_file", "seed_from") if k in cols]
        cols = front + [k for k in cols if k not in front + tail] + tail
        numeric = [k for k in cols if k not in front + tail]
        sheet.append(cols)
        for c in sheet[1]:
            c.font, c.fill = bold, head
            c.alignment = Alignment(wrap_text=True, vertical="center")
        for r in rows:
            sheet.append([r.get(k, "") for k in cols])
            all_rows.append({"condition": cname, **r})
        for row in sheet.iter_rows(min_row=2, max_row=len(rows) + 1):
            for c in row:
                c.font = reg
        # formula rows
        seeds = list(dict.fromkeys(r["seed"] for r in rows))
        blocks = []
        start = 2
        for sd in seeds:
            n = sum(1 for r in rows if r["seed"] == sd)
            blocks.append((sd, start, start + n - 1))
            start += n
        last = len(rows) + 1
        r0 = last + 2
        out_rows = []
        for sd, a, b in blocks:
            out_rows.append((f"seed {sd} MEAN", "AVERAGE", a, b))
            out_rows.append((f"seed {sd} SD", "STDEV", a, b))
        out_rows.append((f"all {len(rows)} models MEAN", "AVERAGE", 2, last))
        out_rows.append((f"all {len(rows)} models SD", "STDEV", 2, last))
        seed_mean_rows = []
        for i, (label, fn, a, b) in enumerate(out_rows):
            rr = r0 + i
            sheet.cell(rr, 1, label).font = bold
            for j, k in enumerate(cols, 1):
                if k in numeric:
                    L = get_column_letter(j)
                    cell = sheet.cell(rr, j, f"=IFERROR({fn}({L}{a}:{L}{b}),\"\")")
                    cell.number_format = "0.000"
                    cell.font, cell.fill = reg, stat
            if label.endswith("MEAN") and label.startswith("seed"):
                seed_mean_rows.append(rr)
        rr = r0 + len(out_rows)
        if len(seed_mean_rows) > 1:
            for label, fn in (("mean of seed means", "AVERAGE"), ("SD of seed means (n = seeds)", "STDEV")):
                sheet.cell(rr, 1, label).font = bold
                for j, k in enumerate(cols, 1):
                    if k in numeric:
                        L = get_column_letter(j)
                        refs = ",".join(f"{L}{x}" for x in seed_mean_rows)
                        cell = sheet.cell(rr, j, f"=IFERROR({fn}({refs}),\"\")")
                        cell.number_format = "0.000"
                        cell.font, cell.fill = reg, stat
                rr += 1
        sheet.cell(rr + 1, 1, "Per-model values copied from the AF3 file in 'source_file'; grey rows are "
                              "Excel formulas.").font = ital
        sheet.freeze_panes = "D2"
        sheet.column_dimensions["A"].width = 22
        sheet.column_dimensions["C"].width = 34
        for j, k in enumerate(cols, 1):
            if k in numeric:
                sheet.column_dimensions[get_column_letter(j)].width = 11
            elif k in tail:
                sheet.column_dimensions[get_column_letter(j)].width = 60

    if all_rows:
        sa = wb.create_sheet("All models")
        cols = list(dict.fromkeys(k for r in all_rows for k in r))
        sa.append(cols)
        for c in sa[1]:
            c.font, c.fill = bold, head
        for r in all_rows:
            sa.append([r.get(k, "") for k in cols])
        sa.freeze_panes = "B2"
        sa.auto_filter.ref = sa.dimensions
    out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conditions", default=str(HERE / "conditions.csv"))
    ap.add_argument("--data", default=str(HERE / "data"))
    ap.add_argument("--out", default=str(HERE / "results" / "AF3_raw_values_by_model.xlsx"))
    ap.add_argument("--full-data", action="store_true", help="also read full_data files (pLDDT, PAE means)")
    a = ap.parse_args()
    conds = [r for r in csv.DictReader(open(a.conditions, newline="", encoding="utf-8-sig"))
             if r.get("name") and not r["name"].startswith("#")]
    per_cond, info = {}, []
    for c in conds:
        folder = Path(a.data) / c["folder"]
        if not folder.exists():
            print(f"skip {c['name']}: folder not found ({folder})")
            continue
        models = find_models(folder)
        if not models:
            print(f"skip {c['name']}: no summary_confidences files in {folder}")
            continue
        rows = []
        for i, m in enumerate(models, 1):
            print(f"\r{c['name']}: model {i}/{len(models)}", end="", flush=True)
            rows.append(read_model(m, a.full_data))
        print()
        seeds = list(dict.fromkeys(r["seed"] for r in rows))
        per_seed = "/".join(str(sum(1 for r in rows if r["seed"] == s)) for s in seeds)
        if len(rows) != 15 or len(seeds) != 3:
            print(f"  NOTE {c['name']}: {len(rows)} models, seeds {seeds}")
        per_cond[c["name"]] = rows
        info.append((c["name"], c["folder"], len(rows), ", ".join(seeds), per_seed))
    if not per_cond:
        raise SystemExit("Nothing found. Check --data and conditions.csv.")
    write_workbook(per_cond, info, Path(a.out), a.full_data)
    print(f"wrote {a.out}: {len(per_cond)} conditions")


if __name__ == "__main__":
    main()

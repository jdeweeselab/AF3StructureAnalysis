#!/usr/bin/env python3
"""
verify_with_published_tools.py - check our analysis code against established, published tools,
on the same AlphaFold 3 model files.

    Check                     Our code                              Independent tool
    ------------------------  ------------------------------------  --------------------------------------
    1 Secondary structure     run_ss_plddt.py (assign_ss)           MDTraj compute_dssp (DSSP 2.2.0
                                                                    implementation; Kabsch & Sander 1983),
                                                                    3-state: H = H/G/I, E = E/B, C = rest
    2 pLDDT                   CA B-factors read by run_ss_plddt.py  Biopython MMCIFParser B-factors, and
                                                                    AF3's own atom_plddts (full_data JSON)
    3 PAE (CTD block)         run_all_pae.py per_model_metrics.csv  mean of AF3's raw token-level 'pae'
                                                                    matrix over CTD tokens (numpy only)
    4 Pairwise CTD RMSD       rmsd_ensemble.py summary              Biopython SVDSuperimposer: every pair
                                                                    of models fitted on the core CA atoms,
                                                                    CTD RMSD in that frame

Run from the folder with conditions.csv, after the other scripts have produced results:
    python -m pip install mdtraj biopython openpyxl
    python verify_with_published_tools.py                       # all conditions
    python verify_with_published_tools.py --only TOP2B_noPTM_noDNA,TOP2B_mitotic_50_75bp
Options: --conditions conditions.csv --data data --results results
         --rmsd-summary PATH   our RMSD summary, if rmsd_ensemble.py wrote it outside --results
         --skip-ss / --skip-rmsd / --skip-pae / --skip-plddt     run only some checks
Needs run_ss_plddt.py in the same folder (its functions are what is being checked).

Writes results\\verification\\:
    verification_report.xlsx   Summary tab + one tab per check (per model and per condition)
    verification_summary.txt   the key numbers, with a sentence for the Methods

Expected outcome, and how to read it:
    1  % identical residue calls per model. Differences come from small differences in hydrogen
       placement and from DSSP rules our 3-state code does not implement (e.g. beta-bulges);
       agreement in the high 90s % is typical for independent DSSP implementations.
    2  should match exactly (max difference 0.00).
    3  no-PTM conditions should match exactly. With PTMs, AF3 gives every atom of a modified
       residue its own token; our code averages those tokens into one residue value, the raw
       token mean does not, so the two differ slightly - both values are reported.
    4  close agreement, not identical: our pairwise values compare models in one common frame
       (each fitted to the ensemble average), here each pair is fitted directly to each other.
"""
import argparse
import csv
import json
import re
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
PTM_PARENT = {"SEP": "SER", "TPO": "THR", "PTR": "TYR"}
TO3 = {"H": "H", "S": "E", "L": "C"}          # our labels -> DSSP 3-state labels


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------------------- files
def find_models(folder):
    """AlphaFold Server layout: <job>_model_<n>.cif with <job>_full_data_<n>.json beside it."""
    out, seen = [], set()
    for p in sorted(folder.rglob("*.cif")):
        m = re.search(r"(.*)_model_(\d+)(?:_\d+)?\.cif$", p.name)
        if not m or (p.parent.name, m.group(2)) in seen:
            continue
        seen.add((p.parent.name, m.group(2)))
        full = sorted(p.parent.glob(f"{m.group(1)}_full_data_{m.group(2)}*.json"))
        out.append(dict(name=f"{p.parent.name}/model{m.group(2)}", job_dir=p.parent.name,
                        job=m.group(1), model=int(m.group(2)), seed=job_seed(p.parent), cif=p,
                        full=full[0] if full else None))
    return out


SEED_RE = re.compile(r"(?:^|[_\-\s])(?:sd|seed)[_\-]?(\d+)(?=$|[_\-\s])", re.I)


def job_seed(job_dir):
    for req in sorted(Path(job_dir).glob("*job_request*.json")):
        try:
            d = json.loads(req.read_text())
            d = d[0] if isinstance(d, list) else d
            if d.get("modelSeeds"):
                return str(d["modelSeeds"][0])
        except Exception:
            pass
    m = SEED_RE.search(Path(job_dir).name)
    return m.group(1) if m else None


def read_conditions(path):
    rows = [r for r in csv.DictReader(open(path, newline="", encoding="utf-8-sig"))
            if r.get("name") and not r["name"].startswith("#")]
    for r in rows:
        r["ctd_lo"], r["ctd_hi"] = (int(x) for x in r["ctd"].split("-"))
        fit = (r.get("fit") or "").strip()
        r["fit_lo"], r["fit_hi"] = (int(x) for x in fit.split("-")) if fit else (None, None)
    return rows


# ----------------------------------------------------------------------------- 1 secondary structure
def write_backbone_pdb(bb_chains, path):
    """Backbone-only PDB (N, CA, C, O) of the protein chains, phospho-residues written as their
    parent amino acid, so MDTraj recognises every residue as protein."""
    n = 1
    with open(path, "w") as f:
        for c, bb in bb_chains.items():
            for k, rn in enumerate(bb["resnum"]):
                res = PTM_PARENT.get(bb["resname"][k], bb["resname"][k])
                for a in ("N", "CA", "C", "O"):
                    x, y, z = bb[a][k]
                    if np.isnan(x):
                        continue
                    f.write(f"ATOM  {n % 100000:5d}  {a:<3s} {res:>3s} {c[0]}{int(rn) % 10000:4d}    "
                            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a[0]:>2s}\n")
                    n += 1
            f.write("TER\n")
        f.write("END\n")


def mdtraj_dssp(bb_chains, tmpdir):
    import mdtraj as md
    pdb = Path(tmpdir) / "backbone.pdb"
    write_backbone_pdb(bb_chains, pdb)
    t = md.load_pdb(str(pdb), standard_names=False)
    ss = md.compute_dssp(t, simplified=True)[0]
    calls = {}
    for res, s in zip(t.topology.residues, ss):
        s = s.decode() if isinstance(s, bytes) else str(s)
        calls[(res.chain.chain_id if hasattr(res.chain, "chain_id") else None, res.resSeq)] = s
    # older MDTraj versions have no chain_id: fall back to chain order
    if any(k[0] is None for k in calls):
        ids = list(bb_chains)
        calls = {(ids[res.chain.index], res.resSeq): (s.decode() if isinstance(s, bytes) else str(s))
                 for res, s in zip(t.topology.residues, ss)}
    return calls


def check_ss(cond, models, tmpdir):
    import run_ss_plddt as ours
    per_model, per_res = [], {}
    confusion = {(a, b): 0 for a in "HEC" for b in "HEC"}
    for m in models:
        bb = ours.read_cif_backbone(m["cif"])
        indep = mdtraj_dssp(bb, tmpdir)
        n_all = same_all = n_ctd = same_ctd = 0
        h_ours = h_ind = 0
        for c, b in bb.items():
            mine = ours.assign_ss(b)
            for rn, s in zip(b["resnum"], mine):
                o, i = TO3[s], indep.get((c, int(rn)), "NA")
                if i not in "HEC":
                    continue
                n_all += 1
                same_all += o == i
                if cond["ctd_lo"] <= rn <= cond["ctd_hi"]:
                    n_ctd += 1
                    same_ctd += o == i
                    h_ours += o == "H"
                    h_ind += i == "H"
                    confusion[(o, i)] += 1
                    d = per_res.setdefault((c, int(rn)), [0, 0])
                    d[0] += o == "H"
                    d[1] += i == "H"
        per_model.append(dict(condition=cond["name"], model=m["name"], residues_compared=n_all,
                              pct_identical_all=round(100 * same_all / max(n_all, 1), 2),
                              CTD_residues=n_ctd, pct_identical_CTD=round(100 * same_ctd / max(n_ctd, 1), 2),
                              CTD_pct_helix_ours=round(100 * h_ours / max(n_ctd, 1), 2),
                              CTD_pct_helix_MDTraj=round(100 * h_ind / max(n_ctd, 1), 2)))
    tot = sum(confusion.values())
    summ = dict(condition=cond["name"], models=len(per_model),
                pct_identical_all_mean=round(np.mean([r["pct_identical_all"] for r in per_model]), 2),
                pct_identical_CTD_mean=round(np.mean([r["pct_identical_CTD"] for r in per_model]), 2),
                pct_identical_CTD_min=round(min(r["pct_identical_CTD"] for r in per_model), 2),
                CTD_pct_helix_ours=round(np.mean([r["CTD_pct_helix_ours"] for r in per_model]), 2),
                CTD_pct_helix_MDTraj=round(np.mean([r["CTD_pct_helix_MDTraj"] for r in per_model]), 2))
    # consensus-level: residues whose helix majority (>= 8 of 15) agrees between the two
    nm = len(models)
    agree = sum((o >= nm / 2) == (i >= nm / 2) for o, i in per_res.values())
    summ["CTD_consensus_helix_agreement_pct"] = round(100 * agree / max(len(per_res), 1), 2)
    conf = [dict(condition=cond["name"], ours=a, MDTraj=b, count=v,
                 pct_of_CTD_calls=round(100 * v / max(tot, 1), 2)) for (a, b), v in confusion.items()]
    return per_model, summ, conf


# ----------------------------------------------------------------------------- 2 pLDDT
def biopython_atoms(cif):
    from Bio.PDB import MMCIFParser
    s = MMCIFParser(QUIET=True).get_structure("m", str(cif))
    model = next(iter(s))
    atoms = list(model.get_atoms())
    ca = {}
    for a in atoms:
        if a.get_id() == "CA":
            res = a.get_parent()
            ca[(res.get_parent().id, res.id[1])] = a.get_bfactor()
    return [a.get_bfactor() for a in atoms], ca


def check_plddt(cond, models):
    import run_ss_plddt as ours
    rows = []
    for m in models:
        bfac, ca_bio = biopython_atoms(m["cif"])
        bb = ours.read_cif_backbone(m["cif"])
        diffs = [abs(b["plddt"][k] - ca_bio[(c, int(rn))])
                 for c, b in bb.items() for k, rn in enumerate(b["resnum"]) if (c, int(rn)) in ca_bio]
        row = dict(condition=cond["name"], model=m["name"], residues_compared=len(diffs),
                   max_abs_diff_ours_vs_Biopython=round(max(diffs), 4) if diffs else "")
        if m["full"]:
            j = json.loads(Path(m["full"]).read_text())
            ap = j.get("atom_plddts") or []
            row["atoms_CIF"], row["atoms_JSON"] = len(bfac), len(ap)
            if len(ap) == len(bfac):
                row["max_abs_diff_CIF_vs_AF3_JSON"] = round(float(np.max(np.abs(np.array(bfac) - np.array(ap)))), 4)
            else:
                row["max_abs_diff_CIF_vs_AF3_JSON"] = "atom counts differ"
        rows.append(row)
    return rows


# ----------------------------------------------------------------------------- 3 PAE
def raw_ctd_pae(full, ctd_lo, ctd_hi, protein_chains):
    j = json.loads(Path(full).read_text())
    ch, rs = j["token_chain_ids"], j["token_res_ids"]
    idx = np.array([k for k, (c, r) in enumerate(zip(ch, rs))
                    if c in protein_chains and ctd_lo <= int(r) <= ctd_hi])
    pae = np.asarray(j["pae"], dtype=np.float32)
    return float(pae[np.ix_(idx, idx)].mean()), len(idx), len(set(zip([ch[i] for i in idx], [rs[i] for i in idx])))


def find_pae_table(results, cond):
    """results/<condition>/per_model_metrics.csv, else any per_model_metrics.csv under results
    whose parent folder is named after the condition."""
    f = Path(results) / cond["name"] / "per_model_metrics.csv"
    if f.exists():
        return f
    hits = [p for p in Path(results).rglob("per_model_metrics.csv")
            if p.parent.name.lower() == cond["name"].lower()]
    return hits[0] if hits else None


def norm_job(x):
    """Comparable job name: lower case, no leading 'fold_', no trailing '/seed..' part."""
    x = x.strip().lower().replace("\\", "/")
    x = re.sub(r"/seed[^/]*$", "", x)
    x = x.split("/")[-1] if "/" in x else x
    return x[5:] if x.startswith("fold_") else x


def ours_pae_table(results, cond):
    """Read our per-model PAE values. Model names can look like
         <job folder>/seed<N>_model<k>        (current run_all_pae.py)
         fold_<job>_model<k>                  (earlier versions)
         <anything>_model_<k> / model<k>
    Returns lookups by (job name, model) and by (seed, model)."""
    f = find_pae_table(results, cond)
    by_job, by_seed, names = {}, {}, []
    if not f:
        return by_job, by_seed, names, None
    for r in csv.DictReader(open(f, newline="", encoding="utf-8-sig")):
        name = (r.get("model") or "").strip()
        names.append(name)
        if not r.get("PAE_CTD_to_CTD") or name in ("MEAN", "SD"):
            continue
        m = re.match(r"^(.*?)[_/]?model_?(\d+)$", name, re.I)
        if not m:
            continue
        prefix, k = m.group(1), int(m.group(2))
        v = float(r["PAE_CTD_to_CTD"])
        seed = None
        sm = re.search(r"seed[_\-]?(\d+)$", prefix, re.I)          # '.../seed1'
        if sm:
            seed = sm.group(1)
        job = norm_job(prefix)
        if job and not re.fullmatch(r"seed[_\-]?\d*\??", job):
            by_job[(job, k)] = v
            if seed is None:
                sm = SEED_RE.search(job)
                seed = sm.group(1) if sm else None
        if seed:
            by_seed.setdefault((seed, k), v)
    return by_job, by_seed, names, f


def check_pae(cond, models, results, protein_chains):
    by_job, by_seed, names, f = ours_pae_table(results, cond)
    if not f:
        log(f"  PAE: no per_model_metrics.csv for {cond['name']} under {Path(results).resolve()} "
            "(run run_all_pae.py, or give --results the folder it wrote to)")
    rows = []
    for m in models:
        if not m["full"]:
            continue
        raw, ntok, nres = raw_ctd_pae(m["full"], cond["ctd_lo"], cond["ctd_hi"], protein_chains)
        mine = by_job.get((norm_job(m["job"]), m["model"]))
        if mine is None:
            mine = by_job.get((norm_job(m["job_dir"]), m["model"]))
        if mine is None and m["seed"]:
            mine = by_seed.get((m["seed"], m["model"]))
        rows.append(dict(condition=cond["name"], model=m["name"], CTD_tokens=ntok, CTD_residues=nres,
                         PAE_CTD_raw_token_mean=round(raw, 3),
                         PAE_CTD_ours=round(mine, 3) if mine is not None else "not found",
                         abs_diff=round(abs(raw - mine), 3) if mine is not None else ""))
    if f and rows and all(r["PAE_CTD_ours"] == "not found" for r in rows):
        log(f"  PAE: could not match models to {f}. Names there look like {names[:2]}; "
            f"model files here look like {[m['job'] + '_model' + str(m['model']) for m in models[:2]]}")
    return rows


# ----------------------------------------------------------------------------- 4 RMSD
def ca_coords(cif, chains=None):
    import run_ss_plddt as ours
    bb = ours.read_cif_backbone(cif)
    keys, xyz = [], []
    for c in (chains or list(bb)):
        b = bb[c]
        for k, rn in enumerate(b["resnum"]):
            keys.append((c, int(rn)))
            xyz.append(b["CA"][k])
    return keys, np.array(xyz)


def find_rmsd_summary(cond, results, rmsd_summary, data):
    names = {cond["name"].lower(), cond["folder"].lower()}
    files = []
    if rmsd_summary:
        files.append(Path(rmsd_summary))
    files.append(Path(results) / "all_conditions_RMSD_summary.csv")
    for root in (Path(results), Path(data), Path(results).parent):
        if root.exists():
            files += sorted(root.rglob("all_conditions_RMSD_summary.csv"),
                            key=lambda p: p.stat().st_mtime, reverse=True)
            files += [p for p in root.rglob("*_RMSD_summary.csv")
                      if p.name.lower()[:-len("_rmsd_summary.csv")] in names]
    seen = set()
    for f in files:
        if f in seen or not f.exists():
            continue
        seen.add(f)
        for r in csv.DictReader(open(f, newline="", encoding="utf-8-sig")):
            if (r.get("condition") or "").lower() in names and r.get("pairwise_CTD_RMSD_mean_A"):
                return r, f
    return {}, None


def check_rmsd(cond, models, results, rmsd_summary=None, data="data"):
    from Bio.SVDSuperimposer import SVDSuperimposer
    keys0, X = None, []
    for m in models:
        keys, xyz = ca_coords(m["cif"])
        if keys0 is None:
            keys0 = keys
        elif keys != keys0:
            raise ValueError(f"{m['name']}: CA atoms differ from first model")
        X.append(xyz)
    resi = np.array([r for _, r in keys0])
    ctd = (resi >= cond["ctd_lo"]) & (resi <= cond["ctd_hi"])
    core = ((resi >= cond["fit_lo"]) & (resi <= cond["fit_hi"]) & ~ctd) if cond["fit_lo"] else ~ctd
    ok = ~np.isnan(np.stack(X)).any(axis=(0, 2))
    ctd, core = ctd & ok, core & ok
    sup = SVDSuperimposer()
    core_r, ctd_r = [], []
    for i in range(len(X)):
        for j in range(i + 1, len(X)):
            sup.set(X[i][core], X[j][core])            # fixed = model i, moving = model j
            sup.run()
            rot, tran = sup.get_rotran()
            core_r.append(sup.get_rms())
            moved = np.dot(X[j][ctd], rot) + tran
            ctd_r.append(float(np.sqrt(((moved - X[i][ctd]) ** 2).sum(1).mean())))
    ours, src = find_rmsd_summary(cond, results, rmsd_summary, data)
    if not ours:
        log(f"  RMSD: no pairwise values for {cond['name']} found (looked for "
            f"all_conditions_RMSD_summary.csv under {Path(results).resolve()} and {Path(data).resolve()}, "
            f"and {cond['name']}_RMSD_summary.csv / {cond['folder']}_RMSD_summary.csv); "
            "use --rmsd-summary PATH")
    pw = lambda k: float(ours[k]) if ours.get(k) else None
    o_mean, o_sd = pw("pairwise_CTD_RMSD_mean_A"), pw("pairwise_CTD_RMSD_sd_A")
    b_mean, b_sd = float(np.mean(ctd_r)), float(np.std(ctd_r, ddof=1))
    return dict(condition=cond["name"], models=len(X), pairs=len(ctd_r),
                core_CA=int(core.sum()), CTD_CA=int(ctd.sum()),
                pairwise_core_RMSD_Biopython=f"{np.mean(core_r):.2f} ± {np.std(core_r, ddof=1):.2f}",
                pairwise_CTD_RMSD_Biopython=f"{b_mean:.2f} ± {b_sd:.2f}",
                pairwise_CTD_RMSD_ours=f"{o_mean:.2f} ± {o_sd:.2f}" if o_mean is not None else "not found",
                pct_difference_of_means=round(100 * (b_mean - o_mean) / o_mean, 1) if o_mean else "",
                our_values_from=str(src) if src else "not found")


# ----------------------------------------------------------------------------- report
def write_report(out_dir, tabs, summary_rows, text):
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    out_dir.mkdir(parents=True, exist_ok=True)
    wb = Workbook()
    ws = wb.active
    ws.title = "Summary"
    bold, head = Font(name="Arial", bold=True), PatternFill("solid", start_color="D9E1F2")
    ws.append(["Check", "Our code", "Independent tool", "Result"])
    for c in ws[1]:
        c.font, c.fill = bold, head
    for r in summary_rows:
        ws.append(r)
    for col, w in zip("ABCD", (26, 44, 52, 70)):
        ws.column_dimensions[col].width = w
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.alignment = Alignment(wrap_text=True, vertical="top")
    for name, rows in tabs.items():
        if not rows:
            continue
        t = wb.create_sheet(name[:31])
        cols = list(dict.fromkeys(k for r in rows for k in r))
        t.append(cols)
        for c in t[1]:
            c.font, c.fill = bold, head
        for r in rows:
            t.append([r.get(k, "") for k in cols])
        t.freeze_panes = "B2"
        for k in range(len(cols)):
            t.column_dimensions[t.cell(1, k + 1).column_letter].width = 18
    wb.save(out_dir / "verification_report.xlsx")
    (out_dir / "verification_summary.txt").write_text(text, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conditions", default=str(HERE / "conditions.csv"))
    ap.add_argument("--data", default=str(HERE / "data"))
    ap.add_argument("--results", default=str(HERE / "results"))
    ap.add_argument("--only", help="comma-separated condition names")
    ap.add_argument("--protein-chains", default="A,B")
    ap.add_argument("--rmsd-summary", help="path to all_conditions_RMSD_summary.csv if it is not under --results")
    for k in ("ss", "plddt", "pae", "rmsd"):
        ap.add_argument(f"--skip-{k}", action="store_true")
    a = ap.parse_args()
    missing = []
    if not a.skip_ss:
        try:
            import mdtraj  # noqa: F401
        except ImportError:
            missing.append("mdtraj")
    if not (a.skip_rmsd and a.skip_plddt):
        try:
            import Bio  # noqa: F401
        except ImportError:
            missing.append("biopython")
    if missing:
        raise SystemExit(f"Please install: python -m pip install {' '.join(missing)}")
    try:
        import run_ss_plddt  # noqa: F401
    except ImportError:
        raise SystemExit("run_ss_plddt.py must be in the same folder as this script.")

    conds = read_conditions(a.conditions)
    if a.only:
        want = {x.strip() for x in a.only.split(",")}
        conds = [c for c in conds if c["name"] in want]
    prot = {x.strip() for x in a.protein_chains.split(",")}
    tabs = {k: [] for k in ("SS by condition", "SS by model", "SS confusion (CTD)",
                            "pLDDT", "PAE", "RMSD pairwise")}
    with tempfile.TemporaryDirectory() as tmp:
        for c in conds:
            folder = Path(a.data) / c["folder"]
            models = find_models(folder) if folder.exists() else []
            if not models:
                log(f"skip {c['name']}: no model CIF files in {folder}")
                continue
            log(f"{c['name']}: {len(models)} models")
            if not a.skip_ss:
                pm, summ, conf = check_ss(c, models, tmp)
                tabs["SS by model"] += pm
                tabs["SS by condition"].append(summ)
                tabs["SS confusion (CTD)"] += conf
                log(f"  SS: {summ['pct_identical_CTD_mean']}% identical CTD calls vs MDTraj DSSP")
            if not a.skip_plddt:
                tabs["pLDDT"] += check_plddt(c, models)
            if not a.skip_pae:
                tabs["PAE"] += check_pae(c, models, a.results, prot)
            if not a.skip_rmsd:
                r = check_rmsd(c, models, a.results, a.rmsd_summary, a.data)
                tabs["RMSD pairwise"].append(r)
                log(f"  RMSD: pairwise CTD {r['pairwise_CTD_RMSD_Biopython']} (Biopython) "
                    f"vs {r['pairwise_CTD_RMSD_ours']} (ours)")

    # summary
    S, lines = [], []
    if tabs["SS by condition"]:
        v = [r["pct_identical_CTD_mean"] for r in tabs["SS by condition"]]
        va = [r["pct_identical_all_mean"] for r in tabs["SS by condition"]]
        n = len(tabs["SS by model"])
        res = (f"{n} models, {len(v)} conditions: identical 3-state calls for {np.mean(va):.1f}% of all "
               f"residues and {np.mean(v):.1f}% of CTD residues (condition means {min(v):.1f}–{max(v):.1f}%)")
        S.append(["Secondary structure", "run_ss_plddt.py (Kabsch–Sander H-bond criterion, 3-state)",
                  "MDTraj compute_dssp, simplified (DSSP 2.2.0 implementation)", res])
        lines.append(f"Secondary structure: {res}.")
    if tabs["pLDDT"]:
        d1 = [r["max_abs_diff_ours_vs_Biopython"] for r in tabs["pLDDT"] if isinstance(r.get("max_abs_diff_ours_vs_Biopython"), float)]
        d2 = [r.get("max_abs_diff_CIF_vs_AF3_JSON") for r in tabs["pLDDT"]]
        d2n = [x for x in d2 if isinstance(x, float)]
        res = (f"{len(tabs['pLDDT'])} models: max |difference| ours vs Biopython = {max(d1) if d1 else 'n/a'}; "
               f"CIF vs AF3 JSON atom_plddts = {max(d2n) if d2n else 'n/a'}"
               + ("" if len(d2n) == len(d2) else f" ({len(d2) - len(d2n)} models not compared: no JSON or atom counts differ)"))
        S.append(["pLDDT", "CA B-factor read by run_ss_plddt.py", "Biopython MMCIFParser; AF3 atom_plddts", res])
        lines.append(f"pLDDT: {res}.")
    if tabs["PAE"]:
        found = [r for r in tabs["PAE"] if isinstance(r["abs_diff"], float)]
        same_tok = [r for r in found if r["CTD_tokens"] == r["CTD_residues"]]
        diff_tok = [r for r in found if r["CTD_tokens"] != r["CTD_residues"]]
        res = (f"no modified residues ({len(same_tok)} models): max |difference| = "
               f"{max((r['abs_diff'] for r in same_tok), default='n/a')} Å; with PTMs ({len(diff_tok)} models): "
               f"max |difference| = {max((r['abs_diff'] for r in diff_tok), default='n/a')} Å "
               "(expected: tokens of modified residues averaged per residue in our code)")
        S.append(["PAE, CTD block", "run_all_pae.py (residue-collapsed)", "raw AF3 token-level pae matrix (numpy)", res])
        lines.append(f"PAE (CTD–CTD mean): {res}.")
    if tabs["RMSD pairwise"]:
        p = [r["pct_difference_of_means"] for r in tabs["RMSD pairwise"] if r["pct_difference_of_means"] != ""]
        res = (f"{len(tabs['RMSD pairwise'])} conditions: Biopython pairwise CTD RMSD differs from ours by "
               f"{min(p):+.1f}% to {max(p):+.1f}% (mean of 105 pairs)" if p else "our summary file not found")
        S.append(["Pairwise CTD RMSD", "rmsd_ensemble.py (common frame, iterative average)",
                  "Biopython SVDSuperimposer (each pair fitted directly on core CA)", res])
        lines.append(f"Pairwise CTD RMSD: {res}.")
    text = "\n".join(lines)
    text += ("\n\nSuggested Methods sentence (edit numbers to match):\n"
             "Analysis scripts were checked against established tools on the same model files: "
             "secondary-structure assignments were compared with MDTraj (compute_dssp, simplified "
             "three-state DSSP), pLDDT values with Biopython and the AF3 JSON output, CTD PAE values "
             "with the raw AF3 PAE matrices, and pairwise CTD RMSD with Biopython (SVDSuperimposer); "
             "results are given in Supplementary Table S_.\n")
    out = Path(a.results) / "verification"
    write_report(out, tabs, S, text)
    print("\n" + text)
    log(f"wrote {out / 'verification_report.xlsx'}")


if __name__ == "__main__":
    main()

"""
rmsd_ensemble.py — CTD and core RMSD for an ensemble of AF3 models (e.g. 3 seeds x 5 models = 15),
measured against the ensemble's AVERAGE structure.

What it calculates (CA atoms of the protein chains)
---------------------------------------------------
1. Superimposes every model on the ordered CORE (all residues outside the CTD), iteratively
   refining an average structure (Kabsch fit, no outlier rejection, so every atom counts).
2. For each model, RMSD to the average structure:
     core_RMSD           core residues (reference: how well the ordered part agrees)
     CTD_RMSD            CTD residues, models still superimposed on the core
                         (variation in CTD shape AND position relative to the enzyme)
     CTD_RMSD_ctd_fit    CTD residues after refitting on the CTD itself (shape only)
   Reported as mean +/- SD over the models.
3. Per-residue RMSF (root-mean-square deviation of each CA from its average position).
4. All-against-all pairwise CTD RMSD (105 pairs for 15 models), mean +/- SD, as a check
   that does not depend on the average structure.

Note: PyMOL's own align/super commands discard "outlier" atoms by default, which makes a
disordered CTD look far more similar than it is. This script does its own fitting instead.

How to run
----------
In PyMOL (File > Run Script..., or type):
    run C:/AF3/rmsd_ensemble.py
    ensemble_rmsd C:/AF3/data/TOP2B_noPTM_noDNA, ctd=1200-1626
    ensemble_rmsd C:/AF3/data/TOP2A_PTM_50_75bp, ctd=1172-1531

All conditions in conditions.csv (same file as the PAE / SS scripts):
    ensemble_rmsd_all C:/AF3/conditions.csv, C:/AF3/data, C:/AF3/results

Without opening the PyMOL window (PowerShell, from C:\\AF3):
    pymol -cq rmsd_ensemble.py -- all conditions.csv data results
    pymol -cq rmsd_ensemble.py -- one data\\TOP2B_noPTM_noDNA 1200-1626 results

Options for ensemble_rmsd:
    ctd=START-END      CTD residue range (required)
    chains=A+B         protein chains to use (default: every chain with CA atoms)
    out=FOLDER         where to write results (default: <condition folder>_RMSD)
    keep=1             keep the superimposed models loaded (for images); default 0 deletes them
    fit=START-END      superimpose on this residue range only (e.g. fit=446-1199 for the TOP2B
                       DNA-binding/cleavage core); residues outside it and outside the CTD (e.g. the
                       ATPase domains) are then reported separately as "other_RMSD".
                       In ensemble_rmsd_all, add fit= at the end, or a "fit" column in conditions.csv

Output (in out folder; name = condition folder name)
    <name>_RMSD_per_model.csv      one row per model: core, CTD, CTD (CTD-fit) RMSD to average
    <name>_RMSD_summary.csv        mean +/- SD of each, plus pairwise CTD RMSD mean +/- SD
    <name>_RMSF_per_residue.csv    chain, residue, RMSF (core-aligned), RMSF (CTD-fit)
    <name>_pairwise_CTD_RMSD.csv   full model x model matrix
"""
import csv
import glob
import os
import re
import sys

import numpy as np

try:
    from pymol import cmd
except ImportError:          # allows the maths to be imported/tested outside PyMOL
    cmd = None


# ----------------------------------------------------------------------------- maths
def kabsch(P, Q):
    """Rotation R and translation t that best map P onto Q (both n x 3): P @ R.T + t."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T


def rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def superimpose_ensemble(X, core, max_iter=20, tol=1e-4):
    """X: models x atoms x 3. Fit every model on `core` atoms to an iteratively refined average.
    Returns (aligned X, average structure, list of (R, t) per model)."""
    Y = X.copy()
    ref = X[0]
    transforms = [None] * len(X)
    for _ in range(max_iter):
        for m in range(len(X)):
            R, t = kabsch(X[m][core], ref[core])
            Y[m] = X[m] @ R.T + t
            transforms[m] = (R, t)
        new_ref = Y.mean(0)
        R, t = kabsch(new_ref[core], ref[core])          # keep the frame from drifting
        new_ref = new_ref @ R.T + t
        shift = rmsd(new_ref[core], ref[core])
        ref = new_ref
        if shift < tol:
            break
    # final pass so every model is fitted to the final average
    for m in range(len(X)):
        R, t = kabsch(X[m][core], ref[core])
        Y[m] = X[m] @ R.T + t
        transforms[m] = (R, t)
    return Y, Y.mean(0), transforms


def analyse(X, keys, ctd_lo, ctd_hi, fit=None):
    """X: models x atoms x 3 (CA), keys: list of (chain, resi). Returns dict of results."""
    resi = np.array([r for _, r in keys])
    ctd = (resi >= ctd_lo) & (resi <= ctd_hi)
    if fit:                                    # superimpose on a chosen region (e.g. DNA-binding/cleavage core)
        core = (resi >= fit[0]) & (resi <= fit[1]) & ~ctd
    else:                                      # default: everything outside the CTD
        core = ~ctd
    other = ~ctd & ~core                       # e.g. ATPase domains when fit= is used
    if ctd.sum() < 3 or core.sum() < 3:
        raise ValueError("CTD range leaves too few CTD or core atoms; check ctd=START-END")
    Y, avg, transforms = superimpose_ensemble(X, core)

    per_model = []
    for m in range(len(Y)):
        Rc, tc = kabsch(Y[m][ctd], avg[ctd])
        per_model.append({
            "core_RMSD": rmsd(Y[m][core], avg[core]),
            "other_RMSD": rmsd(Y[m][other], avg[other]) if other.any() else float("nan"),
            "CTD_RMSD": rmsd(Y[m][ctd], avg[ctd]),
            "CTD_RMSD_ctd_fit": rmsd(Y[m][ctd] @ Rc.T + tc, avg[ctd]),
        })

    rmsf_core_aligned = np.sqrt(((Y - avg) ** 2).sum(2).mean(0))
    # CTD-only fit for per-residue shape variation
    Yc, avg_c, _ = superimpose_ensemble(X[:, ctd, :], np.ones(ctd.sum(), dtype=bool))
    rmsf_ctd_fit = np.full(len(keys), np.nan)
    rmsf_ctd_fit[ctd] = np.sqrt(((Yc - avg_c) ** 2).sum(2).mean(0))

    n = len(Y)
    pair = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            pair[i, j] = pair[j, i] = rmsd(Y[i][ctd], Y[j][ctd])
    iu = np.triu_indices(n, 1)
    return dict(Y=Y, avg=avg, transforms=transforms, per_model=per_model,
                rmsf_core_aligned=rmsf_core_aligned, rmsf_ctd_fit=rmsf_ctd_fit,
                pairwise=pair, pairwise_values=pair[iu], ctd=ctd)


def mean_sd(v):
    v = np.asarray(v, dtype=float)
    return float(v.mean()), float(v.std(ddof=1)) if len(v) > 1 else 0.0


# ----------------------------------------------------------------------------- files
def find_model_cifs(folder):
    found = []
    for p in sorted(glob.glob(os.path.join(folder, "**", "*.cif"), recursive=True)):
        name = os.path.basename(p)
        m = re.search(r"(.*)_model_(\d+)(?:_\d+)?\.cif$", name)            # AlphaFold Server
        if m:
            found.append((f"{m.group(1)}_model{m.group(2)}", p))
            continue
        parent = os.path.basename(os.path.dirname(p))
        if re.search(r"seed-\d+_sample-\d+", parent):                       # local AF3
            found.append((f"{os.path.basename(os.path.dirname(os.path.dirname(p)))}_{parent}", p))
    seen, out = set(), []
    for k, p in found:
        if k not in seen:
            seen.add(k)
            out.append((k, p))
    return out


def safe(name):
    return re.sub(r"[^A-Za-z0-9_]", "_", name)[:60]


# ----------------------------------------------------------------------------- PyMOL part
def _ca_coords(obj, chains):
    sel = f"{obj} and name CA"
    if chains:
        sel += f" and chain {chains}"
    model = cmd.get_model(sel)
    keys = [(a.chain, int(re.sub(r"[^0-9-]", "", a.resi))) for a in model.atom]
    xyz = np.array([a.coord for a in model.atom], dtype=float)
    return keys, xyz


def ensemble_rmsd(folder, ctd, chains="", out="", keep=0, fit="", quiet=0):
    """PyMOL command: ensemble_rmsd folder, ctd=START-END [, chains=A+B, out=FOLDER, keep=1]"""
    folder = folder.strip().strip('"')
    lo, hi = (int(x) for x in str(ctd).split("-"))
    keep = int(keep)
    cifs = find_model_cifs(folder)
    if not cifs:
        print(f"ensemble_rmsd: no model CIF files found in {folder}")
        return None
    name = safe(os.path.basename(os.path.normpath(folder)))
    out = out or os.path.normpath(folder) + "_RMSD"
    os.makedirs(out, exist_ok=True)

    objs, keys0, coords = [], None, []
    short = safe(name)[:30]
    for i, (k, path) in enumerate(cifs, start=1):
        # short, unique object names (long AF3 file names were truncated into duplicates)
        obj = f"{short}_m{i:02d}"
        if obj in cmd.get_names("objects"):
            cmd.delete(obj)
        cmd.load(path, obj, zoom=0)
        objs.append(obj)
        keys, xyz = _ca_coords(obj, chains)
        if keys0 is None:
            keys0 = keys
        elif keys != keys0:
            print(f"ensemble_rmsd: {k} has a different set of CA atoms from the first model; "
                  "do not mix conditions in one folder. Stopping.")
            for o in objs:
                cmd.delete(o)
            return None
        coords.append(xyz)
    X = np.stack(coords)
    print(f"ensemble_rmsd: {name}: {len(objs)} models, {X.shape[1]} CA atoms, CTD {lo}-{hi}")
    if len(objs) != 15:
        print(f"  note: expected 15 models (3 seeds x 5), found {len(objs)}")

    fit_rng = tuple(int(x) for x in str(fit).split("-")) if fit else None
    res = analyse(X, keys0, lo, hi, fit_rng)

    # per-model table
    with open(os.path.join(out, f"{name}_RMSD_per_model.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "pymol_object", "core_RMSD_A", "other_noncore_RMSD_A", "CTD_RMSD_A", "CTD_RMSD_ctd_fit_A"])
        for (k, _), obj, pm in zip(cifs, objs, res["per_model"]):
            w.writerow([k, obj, f"{pm['core_RMSD']:.2f}", f"{pm['other_RMSD']:.2f}", f"{pm['CTD_RMSD']:.2f}", f"{pm['CTD_RMSD_ctd_fit']:.2f}"])
    # summary
    summary = {"condition": name, "n_models": len(objs), "ctd": f"{lo}-{hi}"}
    summary["fit_region"] = f"{fit_rng[0]}-{fit_rng[1]}" if fit_rng else "all non-CTD"
    for key in ("core_RMSD", "other_RMSD", "CTD_RMSD", "CTD_RMSD_ctd_fit"):
        m, s = mean_sd([pm[key] for pm in res["per_model"]])
        summary[f"{key}_mean_A"], summary[f"{key}_sd_A"] = round(m, 2), round(s, 2)
    m, s = mean_sd(res["pairwise_values"])
    summary["pairwise_CTD_RMSD_mean_A"], summary["pairwise_CTD_RMSD_sd_A"] = round(m, 2), round(s, 2)
    summary["n_pairs"] = len(res["pairwise_values"])
    with open(os.path.join(out, f"{name}_RMSD_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary))
        w.writeheader()
        w.writerow(summary)
    # per-residue RMSF
    with open(os.path.join(out, f"{name}_RMSF_per_residue.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["chain", "residue", "in_CTD", "RMSF_core_aligned_A", "RMSF_ctd_fit_A"])
        for (c, r), a, b, inctd in zip(keys0, res["rmsf_core_aligned"], res["rmsf_ctd_fit"], res["ctd"]):
            w.writerow([c, r, "yes" if inctd else "", f"{a:.2f}", "" if np.isnan(b) else f"{b:.2f}"])
    # pairwise matrix
    with open(os.path.join(out, f"{name}_pairwise_CTD_RMSD.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + [k for k, _ in cifs])
        for (k, _), row in zip(cifs, res["pairwise"]):
            w.writerow([k] + [f"{v:.2f}" for v in row])

    if keep:
        # move whole models (all atoms, DNA and ligands too) into the core-superimposed frame
        for obj, (R, t) in zip(objs, res["transforms"]):
            xyz = cmd.get_coords(obj)
            cmd.load_coords(xyz @ R.T + t, obj)
        cmd.group(f"{name}_aligned", " ".join(objs))
        cmd.color("grey80", f"{name}_aligned")
        cmd.color("orange", f"{name}_aligned and resi {lo}-{hi}")
    else:
        for o in objs:
            cmd.delete(o)

    print(f"  core RMSD to average      {summary['core_RMSD_mean_A']:.2f} +/- {summary['core_RMSD_sd_A']:.2f} A (fit: {summary['fit_region']})")
    if fit_rng:
        print(f"  other non-CTD residues    {summary['other_RMSD_mean_A']:.2f} +/- {summary['other_RMSD_sd_A']:.2f} A")
    print(f"  CTD RMSD to average       {summary['CTD_RMSD_mean_A']:.2f} +/- {summary['CTD_RMSD_sd_A']:.2f} A (core-superimposed)")
    print(f"  CTD RMSD, CTD refit       {summary['CTD_RMSD_ctd_fit_mean_A']:.2f} +/- {summary['CTD_RMSD_ctd_fit_sd_A']:.2f} A")
    print(f"  pairwise CTD RMSD         {summary['pairwise_CTD_RMSD_mean_A']:.2f} +/- {summary['pairwise_CTD_RMSD_sd_A']:.2f} A "
          f"({summary['n_pairs']} pairs)")
    print(f"  results written to {out}")
    return summary


def ensemble_rmsd_all(conditions, data, results, chains="", fit=""):
    """PyMOL command: ensemble_rmsd_all conditions.csv, data_folder, results_folder"""
    rows = [r for r in csv.DictReader(open(conditions.strip(), newline="", encoding="utf-8-sig"))
            if r.get("name") and not r["name"].startswith("#")]
    summaries = []
    for r in rows:
        folder = os.path.join(data.strip(), r["folder"])
        if not os.path.isdir(folder):
            print(f"skip {r['name']}: no folder {folder}")
            continue
        out = os.path.join(results.strip(), r["name"], "rmsd")
        try:
            s = ensemble_rmsd(folder, r["ctd"], chains=chains, out=out,
                              fit=r.get("fit") or fit)
        except Exception as e:                                     # keep going
            print(f"ERROR {r['name']}: {e}")
            s = None
        if s:
            s["condition"] = r["name"]
            summaries.append(s)
    if summaries:
        path = os.path.join(results.strip(), "all_conditions_RMSD_summary.csv")
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summaries[0]))
            w.writeheader()
            w.writerows(summaries)
        print(f"summary of {len(summaries)} conditions: {path}")


if cmd is not None:
    cmd.extend("ensemble_rmsd", ensemble_rmsd)
    cmd.extend("ensemble_rmsd_all", ensemble_rmsd_all)
    # command-line use:  pymol -cq rmsd_ensemble.py -- all conditions.csv data results
    #                    pymol -cq rmsd_ensemble.py -- one <folder> <START-END> [out]
    if __name__ in ("__main__", "pymol") and len(sys.argv) > 1:
        args = sys.argv[1:]
        if args[0] == "all" and len(args) >= 4:
            ensemble_rmsd_all(args[1], args[2], args[3])
        elif args[0] == "one" and len(args) >= 3:
            ensemble_rmsd(args[1], args[2], out=args[3] if len(args) > 3 else "")

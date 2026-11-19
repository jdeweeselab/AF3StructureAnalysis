#!/usr/bin/env python3
"""
combine_ss_tables.py — one spreadsheet of CTD secondary structure by residue, across all conditions.

Reads the per-residue tables that run_ss_plddt.py already wrote
(results\\<condition>\\secondary_structure\\<condition>_chain<X>_SS.csv) and writes
results\\CTD_SS_by_position.xlsx with:

  <protein> chain <X>   one row per CTD residue; for every condition three columns:
                        SS (consensus H/S/L), Agree (models agreeing, e.g. 14/15), pLDDT (mean).
                        SS cells are coloured red/blue/grey, paler when agreement is low;
                        pLDDT cells use the AlphaFold colour bands.
  Segments              every run of consecutive residues with the same consensus H or S:
                        protein, chain, condition, type, start, end, length, sequence,
                        mean agreement, min agreement, mean pLDDT
  Strands               every residue that ANY model calls strand, per condition:
                        how many models (of 15) call it strand, and the consensus

Usage (from C:\\AF3):   python combine_ss_tables.py
Options:  --conditions conditions.csv  --results results  --out results\\CTD_SS_by_position.xlsx
          --all-residues   include the whole chain, not only the CTD
Requires: python -m pip install openpyxl
"""
import argparse
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
SS_RGB = {"H": (214, 39, 40), "S": (31, 119, 180), "L": (199, 199, 199)}


def faded_hex(ss, frac):
    r, g, b = SS_RGB[ss]
    f = (1 - frac) * 0.85
    return "".join(f"{round(c * (1 - f) + 255 * f):02X}" for c in (r, g, b))


def plddt_hex(v):
    if v >= 90:
        return "0053D6"
    if v >= 70:
        return "65CBF3"
    if v >= 50:
        return "FFDB13"
    return "FF7D45"


def load(results, cond):
    d = Path(results) / cond["name"] / "secondary_structure"
    tables = {}
    for f in sorted(d.glob(f"{cond['name']}_chain*_SS.csv")):
        chain = f.stem[len(cond["name"]) + len("_chain"):-len("_SS")]
        rows = list(csv.DictReader(open(f, newline="")))
        nm = sum(1 for k in rows[0] if k.startswith("m") and k[1:].isdigit())
        tables[chain] = (rows, nm)
    return tables


def segments(rows, nm):
    out, i = [], 0
    while i < len(rows):
        ss = rows[i]["consensus"]
        if ss in ("H", "S"):
            j = i
            while j + 1 < len(rows) and rows[j + 1]["consensus"] == ss and \
                    int(rows[j + 1]["residue"]) == int(rows[j]["residue"]) + 1:
                j += 1
            seg = rows[i:j + 1]
            ag = [int(r["agreement"]) for r in seg]
            pl = [float(r["pLDDT_mean"]) for r in seg]
            out.append(dict(type="helix" if ss == "H" else "strand",
                            start=int(seg[0]["residue"]), end=int(seg[-1]["residue"]),
                            length=len(seg), sequence="".join(r["aa"] for r in seg),
                            mean_agree=f"{sum(ag) / len(ag):.1f}/{nm}", min_agree=f"{min(ag)}/{nm}",
                            mean_pLDDT=round(sum(pl) / len(pl), 1)))
            i = j + 1
        else:
            i += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conditions", default=str(HERE / "conditions.csv"))
    ap.add_argument("--results", default=str(HERE / "results"))
    ap.add_argument("--out", default="")
    ap.add_argument("--all-residues", action="store_true")
    a = ap.parse_args()
    from openpyxl import Workbook
    from openpyxl.styles import PatternFill, Font, Alignment
    from openpyxl.utils import get_column_letter

    fill = lambda h: PatternFill(start_color=h, end_color=h, fill_type="solid")
    conds = [r for r in csv.DictReader(open(a.conditions, newline="", encoding="utf-8-sig"))
             if r.get("name") and not r["name"].startswith("#")]
    data = {}                                  # protein -> list of (cond, tables)
    for c in conds:
        t = load(a.results, c)
        if t:
            data.setdefault(c["protein"], []).append((c, t))
        else:
            print(f"skip {c['name']}: no per-residue tables in {a.results}")

    wb = Workbook()
    wb.remove(wb.active)
    seg_rows, strand_rows = [], []
    for prot, items in data.items():
        chains = sorted({ch for _, t in items for ch in t})
        for ch in chains:
            ws = wb.create_sheet(f"{prot} chain {ch}"[:31])
            present = [(c, t[ch]) for c, t in items if ch in t]
            lo, hi = (int(x) for x in present[0][0]["ctd"].split("-"))
            ref_rows = present[0][1][0]
            keep = [r for r in ref_rows if a.all_residues or lo <= int(r["residue"]) <= hi]
            # header rows
            ws.cell(1, 1, "Residue").font = Font(bold=True)
            ws.cell(1, 2, "AA").font = Font(bold=True)
            col = 3
            for c, _ in present:
                ws.merge_cells(start_row=1, start_column=col, end_row=1, end_column=col + 2)
                h = ws.cell(1, col, c["name"])
                h.font, h.alignment = Font(bold=True), Alignment(horizontal="center")
                for k, lab in enumerate(("SS", "Agree", "pLDDT")):
                    ws.cell(2, col + k, lab).font = Font(bold=True)
                col += 3
            by_cond = [{int(r["residue"]): r for r in rows} for _, (rows, _) in present]
            for i, r in enumerate(keep, start=3):
                res = int(r["residue"])
                ws.cell(i, 1, res)
                ws.cell(i, 2, r["aa"])
                col = 3
                for (c, (rows, nm)), idx in zip(present, by_cond):
                    rr = idx.get(res)
                    if rr:
                        ss, ag, pl = rr["consensus"], int(rr["agreement"]), float(rr["pLDDT_mean"])
                        x = ws.cell(i, col, ss)
                        x.fill = fill(faded_hex(ss, ag / nm))
                        x.alignment = Alignment(horizontal="center")
                        if ss in ("H", "S") and ag / nm >= 0.8:
                            x.font = Font(bold=True, color="FFFFFF")
                        ws.cell(i, col + 1, f"{ag}/{nm}")
                        p = ws.cell(i, col + 2, pl)
                        p.fill = fill(plddt_hex(pl))
                        if pl >= 90:
                            p.font = Font(color="FFFFFF")
                    col += 3
            ws.freeze_panes = "C3"
            for k in range(3, col):
                ws.column_dimensions[get_column_letter(k)].width = 7
            # segments + strands for this chain
            for c, (rows, nm) in present:
                ctd_rows = [r for r in rows if a.all_residues or lo <= int(r["residue"]) <= hi]
                for s in segments(ctd_rows, nm):
                    seg_rows.append(dict(protein=prot, chain=ch, condition=c["name"], **s))
                for r in ctd_rows:
                    if int(r["n_S"]) > 0:
                        strand_rows.append(dict(protein=prot, chain=ch, condition=c["name"],
                                                residue=int(r["residue"]), aa=r["aa"],
                                                models_calling_strand=f"{r['n_S']}/{nm}",
                                                consensus=r["consensus"], agreement=f"{r['agreement']}/{nm}",
                                                pLDDT=float(r["pLDDT_mean"])))

    for title, rows in (("Segments", seg_rows), ("Strands", strand_rows)):
        ws = wb.create_sheet(title)
        if not rows:
            ws.append(["none found"])
            continue
        heads = list(rows[0])
        ws.append(heads)
        for c in ws[1]:
            c.font = Font(bold=True)
        for r in rows:
            ws.append([r[h] for h in heads])
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
    out = a.out or str(Path(a.results) / "CTD_SS_by_position.xlsx")
    wb.save(out)
    print(f"wrote {out}: {len(seg_rows)} segments, {len(strand_rows)} strand rows")


if __name__ == "__main__":
    main()

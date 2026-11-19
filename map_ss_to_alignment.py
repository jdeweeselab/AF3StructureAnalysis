#!/usr/bin/env python3
"""
map_ss_to_alignment.py — colour AF3 consensus secondary structure onto a TOP2A/TOP2B CTD alignment.

Inputs
  --aln     SnapGene protein alignment (.praln) holding a global alignment of two sequences,
            or an aligned FASTA (two entries, gaps as '-').
  --panel   LABEL  TOP2A_chainA_SS.csv  TOP2B_chainA_SS.csv     (repeat for each panel)
  --trackA  LABEL  TOP2A_chainA_SS.csv   /  --trackB  LABEL  TOP2B_chainA_SS.csv
            tracks for one protein only, so the two sides can differ (e.g. 1 TOP2A PTM set
            above, 3 TOP2B PTM sets below); repeat as needed
            the per-residue tables written by run_ss_plddt.py (any chain file works)
  Optional: --top2a-start / --top2b-start  first residue number of each aligned sequence
            (read from the .praln names, e.g. "TOP2A (p1175 .. p1531)", if not given)
            --names "TOP2A" "TOP2B"   labels;  --width 60  alignment columns per line
            --out PREFIX               output file prefix (default: CTD_alignment_SS)

Outputs  (SVG / PDF / PNG for each figure)
  PREFIX_tracks       main-text figure: per line, the two sequences face each other (identities
                      marked |), with one thin secondary-structure track per panel above TOP2A and
                      below TOP2B. Colour = consensus (red helix, blue strand, grey loop), paler =
                      fewer models agree; black dot = mean pLDDT >= 70; triangle = phosphorylated
                      residue in that panel's model.
  PREFIX_<panel>      one supplementary map per panel: each letter on its secondary-structure
                      colour (same fading, dots and triangles)
  PREFIX.xlsx         one row per alignment column: residue numbers, amino acids, identity, and
                      for every panel the consensus, agreement (n/15) and pLDDT for both proteins

Example (PowerShell, one line):
  python map_ss_to_alignment.py --aln TOP2A_TOP2B_CTD.praln
     --panel "No PTM, no DNA" A_noPTM_noDNA.csv B_noPTM_noDNA.csv
     --panel "No PTM, 50/75 bp" A_noPTM_50_75.csv B_noPTM_50_75.csv
     --panel "PTM, 50/75 bp" A_PTM_50_75.csv B_mitotic_50_75.csv
Requires numpy, matplotlib, openpyxl.
"""
import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
SANS = "Arial" if any("Arial" in f.name for f in matplotlib.font_manager.fontManager.ttflist) else "DejaVu Sans"
MONO = "Courier New" if any("Courier New" in f.name for f in matplotlib.font_manager.fontManager.ttflist) else "DejaVu Sans Mono"
SS_COL = {"H": "#d62728", "S": "#1f77b4", "L": "#c7c7c7"}
FADE_MAX = 0.85


# ----------------------------------------------------------------------------- reading
def read_alignment(path):
    """Return [(name, start, gapped_seq), (name, start, gapped_seq)] for the global alignment."""
    p = Path(path)
    raw = p.read_bytes()
    if p.suffix.lower() in (".fa", ".fasta", ".aln", ".txt") or raw.lstrip().startswith(b">"):
        entries, name, seq = [], None, []
        for line in raw.decode().splitlines():
            if line.startswith(">"):
                if name:
                    entries.append((name, "".join(seq)))
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
        entries.append((name, "".join(seq)))
        out = []
        for n, s in entries[:2]:
            m = re.search(r"p?(\d+)\s*\.\.", n)
            out.append((n.split()[0], int(m.group(1)) if m else None, s))
        return out
    # SnapGene .praln: names/starts in the JSON header, sequences as long letter runs
    names = []
    m = re.search(rb'"sequences":(\[.*?\])', raw)
    if m:
        for s in json.loads(m.group(1)):
            mm = re.search(r"p(\d+)", s["name"])
            names.append((s["name"].split()[0], int(mm.group(1)) if mm else None))
    runs = [r.decode() for r in re.findall(rb"[A-Z\-]{50,}", raw)]
    plain = [r for r in runs if "-" not in r]
    gapped = [r for r in runs if "-" in r]
    # the global alignment = two equal-length gapped strings whose ungapped forms are the full sequences
    full = set(plain)
    for i in range(len(gapped) - 1, 0, -1):
        a, b = gapped[i - 1], gapped[i]
        if len(a) == len(b) and a.replace("-", "") in full and b.replace("-", "") in full:
            pair = [a, b]
            break
    else:
        raise SystemExit("Could not find a global alignment of the two full sequences in the .praln file. "
                         "Export the alignment as aligned FASTA and pass that instead.")
    out = []
    for s in pair:
        ung = s.replace("-", "")
        idx = plain.index(ung) if ung in plain else None
        nm, st = names[idx] if idx is not None and idx < len(names) else (f"seq{len(out) + 1}", None)
        out.append((nm, st, s))
    return out


def read_ss(path):
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8-sig")))
    nm = sum(1 for k in rows[0] if re.fullmatch(r"m\d+", k))
    return {int(r["residue"]): r for r in rows}, nm


def map_columns(gapped, start, ss):
    """Per alignment column: dict(aa, res, ss, frac, plddt, phospho) or None for a gap."""
    table, nm = ss
    cols, r = [], start
    for ch in gapped:
        if ch == "-":
            cols.append(None)
            continue
        row = table.get(r)
        if row is None:
            cols.append(dict(aa=ch, res=r, ss=None, frac=0, plddt=np.nan, phospho=False))
        else:
            if row["aa"].upper() != ch:
                raise SystemExit(f"Residue {r}: alignment has {ch} but the SS table has {row['aa']} - "
                                 "check the start numbers (--top2a-start / --top2b-start).")
            cols.append(dict(aa=ch, res=r, ss=row["consensus"], frac=int(row["agreement"]) / nm,
                             agree=f"{row['agreement']}/{nm}", plddt=float(row["pLDDT_mean"]),
                             phospho=row["aa"].islower()))
        r += 1
    return cols


def faded(ss, frac):
    c = np.array(to_rgb(SS_COL[ss]))
    f = (1 - frac) * FADE_MAX
    return tuple(c * (1 - f) + f)


# ----------------------------------------------------------------------------- drawing
def draw_cell(ax, x, y, h, col, show_dot, show_tri):
    if col is None:
        return
    if col["ss"]:
        ax.add_patch(Rectangle((x, y), 1, h, facecolor=faded(col["ss"], col["frac"]), lw=0))
    if show_dot and col["plddt"] >= 70:
        ax.add_patch(plt.Circle((x + 0.5, y + h / 2), min(0.16, h * 0.28), color="black", lw=0))
    if show_tri and col["phospho"]:
        ax.add_patch(Polygon([(x + 0.2, y + h), (x + 0.8, y + h), (x + 0.5, y + h - 0.35 * h)],
                             closed=True, color="#7b2cbf", lw=0))


def first_res(cols, a, b):
    for c in cols[a:b]:
        if c:
            return c["res"]
    return ""


def legend(ax, x0, y0, extra=True):
    items = [("Helix", SS_COL["H"]), ("Strand", SS_COL["S"]), ("Loop", SS_COL["L"])]
    x = x0
    for lab, c in items:
        ax.add_patch(Rectangle((x, y0), 1.2, 0.7, facecolor=c, lw=0))
        ax.text(x + 1.5, y0 + 0.35, lab, va="center", fontsize=5.5, family=SANS)
        x += 6.5
    ax.add_patch(Rectangle((x, y0), 1.2, 0.7, facecolor=faded("H", 8 / 15), lw=0))
    ax.text(x + 1.5, y0 + 0.35, "paler = fewer of 15 models agree", va="center", fontsize=5.5, family=SANS)
    x += 19
    ax.add_patch(plt.Circle((x + 0.6, y0 + 0.35), 0.2, color="black"))
    ax.text(x + 1.3, y0 + 0.35, "pLDDT ≥ 70", va="center", fontsize=5.5, family=SANS)
    x += 8
    if extra:
        ax.add_patch(Polygon([(x + 0.2, y0 + 0.7), (x + 1.0, y0 + 0.7), (x + 0.6, y0 + 0.2)], color="#7b2cbf"))
        ax.text(x + 1.3, y0 + 0.35, "phosphorylated in model", va="center", fontsize=5.5, family=SANS)
    x += 13
    ax.text(x, y0 + 0.35, "| identical", va="center", fontsize=5.5, family=SANS)


def tracks_figure(seqA, seqB, tracksA, tracksB, width, out, formats):
    """tracksA / tracksB: lists of (label, cols) drawn above TOP2A / below TOP2B (any number each)."""
    nA, sA, gA = seqA
    nB, sB, gB = seqB
    n = len(gA)
    blocks = (n + width - 1) // width
    th, sh, gap_id, bgap = 0.55, 1.0, 0.45, 1.3
    kA, kB = len(tracksA), len(tracksB)
    block_h = (kA + kB) * th + 2 * sh + gap_id + bgap
    lab_w = 15
    fig_h = 0.105 * (blocks * block_h + 2.5)
    fig, ax = plt.subplots(figsize=(7.0, fig_h))
    y = 0
    for bi in range(blocks):
        a, b = bi * width, min(n, (bi + 1) * width)
        for pi, (lab, cols) in enumerate(tracksA):          # first track at the top
            ty = y - (pi + 1) * th
            ax.text(-1.2, ty + th / 2, lab, ha="right", va="center", fontsize=4.6, family=SANS, color="#444")
            for j, c in enumerate(cols[a:b]):
                draw_cell(ax, j, ty + 0.06, th - 0.12, c, True, True)
        ya = y - kA * th - sh
        ax.text(-lab_w, ya + sh / 2, nA, ha="left", va="center", fontsize=6, family=SANS, weight="bold")
        ax.text(-1.2, ya + sh / 2, str(first_res(tracksA[0][1], a, b)), ha="right", va="center", fontsize=5.5, family=SANS)
        for j in range(b - a):
            ax.text(j + 0.5, ya + sh / 2, gA[a + j], ha="center", va="center", fontsize=5.6, family=MONO)
        yi = ya - gap_id
        for j in range(b - a):
            if gA[a + j] == gB[a + j] and gA[a + j] != "-":
                ax.add_line(Line2D([j + 0.5, j + 0.5], [yi + 0.08, yi + gap_id - 0.08], color="#666", lw=0.5))
        yb = yi - sh
        ax.text(-lab_w, yb + sh / 2, nB, ha="left", va="center", fontsize=6, family=SANS, weight="bold")
        ax.text(-1.2, yb + sh / 2, str(first_res(tracksB[0][1], a, b)), ha="right", va="center", fontsize=5.5, family=SANS)
        for j in range(b - a):
            ax.text(j + 0.5, yb + sh / 2, gB[a + j], ha="center", va="center", fontsize=5.6, family=MONO)
        for pi, (lab, cols) in enumerate(tracksB):          # first track nearest the sequence
            ty = yb - (pi + 1) * th
            ax.text(-1.2, ty + th / 2, lab, ha="right", va="center", fontsize=4.6, family=SANS, color="#444")
            for j, c in enumerate(cols[a:b]):
                draw_cell(ax, j, ty + 0.06, th - 0.12, c, True, True)
        y = yb - kB * th - bgap
    legend(ax, 0, y - 0.6)
    ax.set_xlim(-lab_w - 0.5, width + 0.5)
    ax.set_ylim(y - 1.2, 0.3)
    ax.axis("off")
    fig.tight_layout(pad=0.2)
    for f in formats:
        fig.savefig(f"{out}.{f}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def letters_figure(seqA, seqB, panel, width, out, formats):
    lab, colsA, colsB = panel
    nA, _, gA = seqA
    nB, _, gB = seqB
    n = len(gA)
    blocks = (n + width - 1) // width
    sh, gap_id, bgap = 1.0, 0.4, 1.2
    block_h = 2 * sh + gap_id + bgap
    fig, ax = plt.subplots(figsize=(7.0, 0.11 * (blocks * block_h + 3)))
    y = 0
    for bi in range(blocks):
        a, b = bi * width, min(n, (bi + 1) * width)
        for nm, g, cols, yy in ((nA, gA, colsA, y - sh), (nB, gB, colsB, y - 2 * sh - gap_id)):
            ax.text(-9, yy + sh / 2, nm, ha="left", va="center", fontsize=6, family=SANS, weight="bold")
            ax.text(-1.2, yy + sh / 2, str(first_res(cols, a, b)), ha="right", va="center", fontsize=5.5, family=SANS)
            for j in range(b - a):
                c = cols[a + j]
                if c and c["ss"]:
                    ax.add_patch(Rectangle((j, yy + 0.05), 1, sh - 0.1, facecolor=faded(c["ss"], c["frac"]), lw=0))
                dark = c and c["ss"] in ("H", "S") and c["frac"] >= 0.6
                ax.text(j + 0.5, yy + sh / 2 + 0.05, g[a + j], ha="center", va="center", fontsize=5.6,
                        family=MONO, color="white" if dark else "black")
                if c and c["plddt"] >= 70:
                    ax.add_patch(plt.Circle((j + 0.5, yy + 0.16), 0.09, color="black" if not dark else "white"))
                if c and c["phospho"]:
                    ax.add_patch(Polygon([(j + 0.22, yy + sh + 0.28), (j + 0.78, yy + sh + 0.28),
                                          (j + 0.5, yy + sh + 0.03)], color="#7b2cbf", lw=0))
        yi = y - sh - gap_id
        for j in range(b - a):
            if gA[a + j] == gB[a + j] and gA[a + j] != "-":
                ax.add_line(Line2D([j + 0.5, j + 0.5], [yi + 0.06, yi + gap_id - 0.06], color="#666", lw=0.5))
        y -= block_h
    ax.text(0, 1.0, lab, fontsize=7, family=SANS, weight="bold")
    legend(ax, 0, y - 0.3)
    ax.set_xlim(-9.5, width + 0.5)
    ax.set_ylim(y - 1.0, 1.8)
    ax.axis("off")
    fig.tight_layout(pad=0.2)
    for f in formats:
        fig.savefig(f"{out}.{f}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_xlsx(seqA, seqB, tracksA, tracksB, out):
    try:
        from openpyxl import Workbook
        from openpyxl.styles import PatternFill, Font, Alignment
    except ImportError:
        print("openpyxl not installed: skipping Excel")
        return
    hexc = lambda rgb: "".join(f"{round(v * 255):02X}" for v in rgb)
    fill = lambda h: PatternFill(start_color=h, end_color=h, fill_type="solid")
    wb = Workbook()
    ws = wb.active
    ws.title = "Aligned SS"
    nA, _, gA = seqA
    nB, _, gB = seqB
    head = ["Aln col", f"{nA} res", f"{nA} aa", "Identical", f"{nB} aa", f"{nB} res"]
    tracks = [(nA, lab, cols) for lab, cols in tracksA] + [(nB, lab, cols) for lab, cols in tracksB]
    for prot, lab, _ in tracks:
        head += [f"{prot} {lab} SS", "agree", "pLDDT"]
    ws.append(head)
    for c in ws[1]:
        c.font = Font(bold=True)
        c.alignment = Alignment(wrap_text=True, vertical="top")
    for i in range(len(gA)):
        cA, cB = tracksA[0][1][i], tracksB[0][1][i]
        row = [i + 1, cA["res"] if cA else "", gA[i], "|" if gA[i] == gB[i] and gA[i] != "-" else "",
               gB[i], cB["res"] if cB else ""]
        cells = []
        for prot, lab, cols in tracks:
            c = cols[i]
            if c and c["ss"]:
                row += [c["ss"], c["agree"], c["plddt"]]
                cells.append((prot, c))
            else:
                row += ["", "", ""]
                cells.append((prot, None))
        ws.append(row)
        r = ws.max_row
        for k, (prot, c) in enumerate(cells):
            if c:
                cell = ws.cell(r, 7 + 3 * k)
                cell.fill = fill(hexc(faded(c["ss"], c["frac"])))
                cell.alignment = Alignment(horizontal="center")
                if c["phospho"]:
                    cell.font = Font(color="7B2CBF", bold=True, underline="single")
    ws.freeze_panes = "G2"
    ws.auto_filter.ref = ws.dimensions
    ws2 = wb.create_sheet("Key")
    for line in ("SS: H helix, S strand, L loop (consensus of the models); fill paler = fewer models agree.",
                 "agree: models giving the consensus call; pLDDT: mean over models.",
                 "Purple underlined SS = residue phosphorylated in that condition's model.",
                 "Identical: | where the two aligned residues are the same."):
        ws2.append([line])
    wb.save(f"{out}.xlsx")


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--panel", nargs=3, action="append", metavar=("LABEL", "TOP2A_CSV", "TOP2B_CSV"),
                    help="same condition for both proteins: a track above TOP2A and below TOP2B, plus a single map")
    ap.add_argument("--trackA", nargs=2, action="append", metavar=("LABEL", "TOP2A_CSV"),
                    help="extra track above TOP2A only (repeatable)")
    ap.add_argument("--trackB", nargs=2, action="append", metavar=("LABEL", "TOP2B_CSV"),
                    help="extra track below TOP2B only (repeatable)")
    ap.add_argument("--top2a-start", type=int)
    ap.add_argument("--top2b-start", type=int)
    ap.add_argument("--names", nargs=2)
    ap.add_argument("--width", type=int, default=60)
    ap.add_argument("--out", default="CTD_alignment_SS")
    ap.add_argument("--format", default="svg,pdf,png")
    a = ap.parse_args()
    formats = [f.strip() for f in a.format.split(",")]

    (n1, s1, g1), (n2, s2, g2) = read_alignment(a.aln)
    # first sequence = TOP2A, second = TOP2B unless the names say otherwise
    if "2B" in n1.upper() and "2A" in n2.upper():
        (n1, s1, g1), (n2, s2, g2) = (n2, s2, g2), (n1, s1, g1)
    s1 = a.top2a_start or s1
    s2 = a.top2b_start or s2
    if not s1 or not s2:
        raise SystemExit("Start residue numbers not found in the alignment file; give --top2a-start and --top2b-start.")
    if a.names:
        n1, n2 = a.names
    seqA, seqB = (n1, s1, g1), (n2, s2, g2)
    print(f"alignment: {n1} from {s1} ({len(g1.replace('-', ''))} aa), {n2} from {s2} "
          f"({len(g2.replace('-', ''))} aa), {len(g1)} columns")

    tracksA, tracksB, panels = [], [], []
    for lab, fa, fb in (a.panel or []):
        p = (lab, map_columns(g1, s1, read_ss(fa)), map_columns(g2, s2, read_ss(fb)))
        panels.append(p)
        tracksA.append((lab, p[1]))
        tracksB.append((lab, p[2]))
    for lab, fa in (a.trackA or []):
        tracksA.append((lab, map_columns(g1, s1, read_ss(fa))))
    for lab, fb in (a.trackB or []):
        tracksB.append((lab, map_columns(g2, s2, read_ss(fb))))
    if not tracksA or not tracksB:
        raise SystemExit("Give at least one track for each protein (--panel, or --trackA and --trackB).")
    tracks_figure(seqA, seqB, tracksA, tracksB, a.width, f"{a.out}_tracks", formats)
    for p in panels:
        tag = re.sub(r"[^A-Za-z0-9]+", "_", p[0]).strip("_")
        letters_figure(seqA, seqB, p, a.width, f"{a.out}_{tag}", formats)
    write_xlsx(seqA, seqB, tracksA, tracksB, a.out)
    print(f"wrote {a.out}_tracks.*" + (", one map per panel" if panels else "") + f", and {a.out}.xlsx")


if __name__ == "__main__":
    main()

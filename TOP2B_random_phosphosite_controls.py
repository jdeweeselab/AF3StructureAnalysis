#!/usr/bin/env python3
"""
Random control phosphosite sets for the human TOP2B C-terminal domain.

Purpose
-------
Generates random sets of Ser/Thr/Tyr residues in the C-terminal region of
human DNA topoisomerase II beta (UniProt Q02880, canonical isoform beta-2,
1626 aa). These serve as controls for AlphaFold3 models phosphorylated at
experimentally observed sites.

Method
------
1. Candidate pool: every Ser, Thr and Tyr in residues 1200-1626 of Q02880
   (76 residues).
2. Excluded: every residue in the experimentally derived sets A, B and C
   defined below (23 distinct residues). Pool = 53 residues.
   Residues reported phosphorylated in other sources but not in sets A-C are
   NOT excluded, so the pool is "random CTD S/T/Y outside the modelled sets",
   not "residues never observed to be phosphorylated".
3. One pseudorandom generator, random.Random(seed=20261014), draws sets of
   10, 18 and 19 residues in that order with random.sample() (sampling without
   replacement within each set; sets are drawn independently and may share
   residues with each other).

Because the draws are sequential from one generator, the first set (10
sites) is the same whether or not the later sets are drawn.

Sets used in the manuscript
---------------------------
Random control (Table 1): first 10-site draw with the default seed 20261014
    python3 TOP2B_random_phosphosite_controls.py --sizes 10
Random A, B, C (Supplementary Information): one 10-site draw per seed
    python3 TOP2B_random_phosphosite_controls.py --sizes 10 --seed 20261116
    python3 TOP2B_random_phosphosite_controls.py --sizes 10 --seed 20261117
    python3 TOP2B_random_phosphosite_controls.py --sizes 10 --seed 20261118
The 18- and 19-site draws (default --sizes 10 18 19) come from earlier
exploratory runs and were not used in the manuscript; they are kept so that the
default run reproduces the original output exactly.

Note on the excluded sites: Sets A-C are the sites excluded from the pool.
The normal-liver site T1575 (Bian et al. 2014) is not in Sets A-C and so was
not excluded; none of the random sets drawn contains it.

Experimental sets excluded from the pool (Q02880 numbering)
-----------------------------------------------------------
Set A (18): Detected in nocodazole-arrested mitotic HeLa cells and localized
            (Ascore >= 13); Dephoure et al. 2008, PNAS 105:10762-10767,
            Supplementary Table S2.
Set B (8):  Set A sites also localized (class 1) and quantified in mitosis by
            Olsen et al. 2010, Sci Signal 3:ra3, Supplementary Table S2.
Set C (19): Sites localized in both Dephoure 2008 and Olsen 2010 (any
            cell-cycle stage).

Reproducibility
---------------
Uses only the Python standard library. random.Random with an integer seed and
random.sample() on a list give identical results across CPython 3.x
releases; tested with CPython 3.12 and 3.13. The script compares its draws with the
published sets and reports whether they match.

Usage
-----
    python3 TOP2B_random_phosphosite_controls.py
    python3 TOP2B_random_phosphosite_controls.py --csv random_sets.csv
    python3 TOP2B_random_phosphosite_controls.py --sizes 10
"""

import argparse
import csv
import random
import sys
from collections import Counter

# ---------------------------------------------------------------------------
# Human TOP2B, UniProt Q02880 (canonical isoform beta-2), 1626 residues
# ---------------------------------------------------------------------------
Q02880 = (
    "MAKSGGCGAGAGVGGGNGALTWVTLFDQNNAAKKEESETANKNDSSKKLSVERVYQKKTQ"
    "LEHILLRPDTYIGSVEPLTQFMWVYDEDVGMNCREVTFVPGLYKIFDEILVNAADNKQRD"
    "KNMTCIKVSIDPESNIISIWNNGKGIPVVEHKVEKVYVPALIFGQLLTSSNYDDDEKKVT"
    "GGRNGYGAKLCNIFSTKFTVETACKEYKHSFKQTWMNNMMKTSEAKIKHFDGEDYTCITF"
    "QPDLSKFKMEKLDKDIVALMTRRAYDLAGSCRGVKVMFNGKKLPVNGFRSYVDLYVKDKL"
    "DETGVALKVIHELANERWDVCLTLSEKGFQQISFVNSIATTKGGRHVDYVVDQVVGKLIE"
    "VVKKKNKAGVSVKPFQVKNHIWVFINCLIENPTFDSQTKENMTLQPKSFGSKCQLSEKFF"
    "KAASNCGIVESILNWVKFKAQTQLNKKCSSVKYSKIKGIPKLDDANDAGGKHSLECTLIL"
    "TEGDSAKSLAVSGLGVIGRDRYGVFPLRGKILNVREASHKQIMENAEINNIIKIVGLQYK"
    "KSYDDAESLKTLRYGKIMIMTDQDQDGSHIKGLLINFIHHNWPSLLKHGFLEEFITPIVK"
    "ASKNKQELSFYSIPEFDEWKKHIENQKAWKIKYYKGLGTSTAKEAKEYFADMERHRILFR"
    "YAGPEDDAAITLAFSKKKIDDRKEWLTNFMEDRRQRRLHGLPEQFLYGTATKHLTYNDFI"
    "NKELILFSNSDNERSIPSLVDGFKPGQRKVLFTCFKRNDKREVKVAQLAGSVAEMSAYHH"
    "GEQALMMTIVNLAQNFVGSNNINLLQPIGQFGTRLHGGKDAASPRYIFTMLSTLARLLFP"
    "AVDDNLLKFLYDDNQRVEPEWYIPIIPMVLINGAEGIGTGWACKLPNYDAREIVNNVRRM"
    "LDGLDPHPMLPNYKNFKGTIQELGQNQYAVSGEIFVVDRNTVEITELPVRTWTQVYKEQV"
    "LEPMLNGTDKTPALISDYKEYHTDTTVKFVVKMTEEKLAQAEAAGLHKVFKLQTTLTCNS"
    "MVLFDHMGCLKKYETVQDILKEFFDLRLSYYGLRKEWLVGMLGAESTKLNNQARFILEKI"
    "QGKITIENRSKKDLIQMLVQRGYESDPVKAWKEAQEKAAEEDETQNQHDDSSSDSGTPSG"
    "PDFNYILNMSLWSLTKEKVEELIKQRDAKGREVNDLKRKSPSDLWKEDLAAFVEELDKVE"
    "SQEREDVLAGMSGKAIKGKVGKPKVKKLQLEETMPSPYGRRIIPEITAMKADASKKLLKK"
    "KKGDLDTAAVKVEFDEEFSGAPVEGAGEEALTPSVPINKGPKPKREKKEPGTRVRKTPTS"
    "SGKPSAKKVKKRNPWSDDESKSESDLEETEPVVIPRDSLLRRAAAERPKYTFDFSEEEDD"
    "DADDDDDDNNDLEELKVKASPITNDGEDEFVPSDGLDKDEYTFSPGKSKATPEKSLHDKK"
    "SQDFGNLFSFPSYSQKSEDDSAKFDSNEEDSASVFSPSFGLKQTDKVPSKTVAAKKGKPS"
    "SDTVPKPKRAPKQKKVVEAVNSDSDSEFGIPKKTTTPKGKGRGAKKRKASGSENEGDYNP"
    "GRKTSKTTSKKPKKTSFDQDSDVDIFPSDFPTEPPSLPRTGRARKEVKYFAESDEEEDDV"
    "DFAMFN"
)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
REGION_START = 1200          # first residue of the region sampled (inclusive)
REGION_END = 1626            # last residue (C terminus)
TARGET_RESIDUES = "STY"
SEED = 20261014
DEFAULT_SIZES = (10, 18, 19)

# Experimentally derived sets excluded from the pool (Q02880 numbering)
SET_A = [1236, 1375, 1400, 1413, 1421, 1422, 1424, 1466, 1522, 1524, 1526,
         1550, 1552, 1576, 1581, 1592, 1596, 1613]
SET_B = [1236, 1522, 1524, 1526, 1550, 1552, 1581, 1596]
SET_C = [1236, 1336, 1340, 1342, 1344, 1400, 1413, 1422, 1424, 1466, 1471,
         1522, 1524, 1526, 1550, 1552, 1581, 1596, 1613]

# Draws from the original run with the default seed and sizes (used only to
# confirm reproduction; only the 10-site draw is used in the manuscript)
PUBLISHED = {
    10: [1294, 1317, 1428, 1461, 1476, 1484, 1500, 1534, 1558, 1567],
    18: [1238, 1294, 1325, 1358, 1371, 1403, 1431, 1449, 1453, 1454, 1476,
         1478, 1484, 1501, 1564, 1567, 1600, 1609],
    19: [1201, 1238, 1247, 1267, 1279, 1319, 1320, 1349, 1358, 1428, 1435,
         1441, 1449, 1453, 1491, 1536, 1564, 1569, 1588],
}


def label(pos):
    """One-letter residue plus position, e.g. S1294 (1-based numbering)."""
    return f"{Q02880[pos - 1]}{pos}"


def window(pos, flank=7):
    """Sequence around a residue: flanks in lower case, the residue in upper case."""
    left = Q02880[max(0, pos - 1 - flank):pos - 1].lower()
    right = Q02880[pos:pos + flank].lower()
    return left + Q02880[pos - 1] + right


def check_inputs():
    assert len(Q02880) == 1626, "Q02880 sequence must be 1626 residues"
    for name, sites in (("Set A", SET_A), ("Set B", SET_B), ("Set C", SET_C)):
        bad = [p for p in sites if Q02880[p - 1] not in TARGET_RESIDUES]
        assert not bad, f"{name} contains non-S/T/Y positions: {bad}"
        assert all(REGION_START <= p <= REGION_END for p in sites), f"{name} outside region"
    assert set(SET_B) <= set(SET_A), "Set B must be a subset of Set A"


def build_pool():
    region = [p for p in range(REGION_START, REGION_END + 1)
              if Q02880[p - 1] in TARGET_RESIDUES]
    excluded = set(SET_A) | set(SET_B) | set(SET_C)
    pool = [p for p in region if p not in excluded]   # ascending order
    return region, excluded, pool


def draw_sets(pool, sizes, seed=SEED):
    rng = random.Random(seed)
    return [(n, sorted(rng.sample(pool, n))) for n in sizes]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sizes", type=int, nargs="+", default=list(DEFAULT_SIZES),
                        help="set sizes to draw, in order (default: 10 18 19)")
    parser.add_argument("--seed", type=int, default=SEED,
                        help=f"random seed (default: {SEED})")
    parser.add_argument("--csv", metavar="PATH",
                        help="also write one row per drawn residue to this CSV file")
    args = parser.parse_args()

    check_inputs()
    region, excluded, pool = build_pool()
    draws = draw_sets(pool, args.sizes, args.seed)

    def comp(sites):
        c = Counter(Q02880[p - 1] for p in sites)
        return f"S={c.get('S', 0)} T={c.get('T', 0)} Y={c.get('Y', 0)}"

    print(f"Python {sys.version.split()[0]}")
    print(f"Protein: human TOP2B, UniProt Q02880 ({len(Q02880)} aa)")
    print(f"Region sampled: residues {REGION_START}-{REGION_END}; residue types: {TARGET_RESIDUES}")
    print(f"S/T/Y in region: {len(region)} ({comp(region)})")
    print(f"Excluded (Sets A, B, C): {len(excluded)}")
    print(f"Pool: {len(pool)} ({comp(pool)})")
    print(f"Seed: {args.seed}\n")

    all_reproduced = True
    rows = []
    for i, (n, sites) in enumerate(draws):
        overlap = {name: len(set(sites) & set(s))
                   for name, s in (("A", SET_A), ("B", SET_B), ("C", SET_C))}
        assert not any(overlap.values()), "draw overlaps an excluded set"
        print(f"Random-{n} ({comp(sites)}): " + ", ".join(label(p) for p in sites))
        print("  overlap with Set A/B/C: "
              + ", ".join(f"{k}={v}" for k, v in overlap.items()))
        # Draws are sequential, so a published set is only comparable when the
        # seed and every earlier draw size match the published run.
        comparable = (args.seed == SEED
                      and list(args.sizes[:i + 1]) == list(DEFAULT_SIZES[:i + 1]))
        if comparable:
            same = sites == PUBLISHED[n]
            all_reproduced &= same
            print(f"  matches published Random-{n}: {'yes' if same else 'NO'}")
        for p in sites:
            rows.append({"set": f"Random-{n}", "site": label(p), "position": p,
                         "residue": Q02880[p - 1], "sequence_window_pm7": window(p)})
        print()

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {args.csv}")

    if not all_reproduced:
        print("WARNING: draws differ from the published sets.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

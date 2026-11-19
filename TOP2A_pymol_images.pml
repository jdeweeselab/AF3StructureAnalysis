# TOP2A MEGA CODE
# Folds 0-14 | Alignment, pLDDT coloring, and CSV export
# ============================================================


# ============================================================
# SECTION 1: ALIGNMENT
# Align all models to TOP2AFold0
# ============================================================

align TOP2AFold1, TOP2AFold0
align TOP2AFold2, TOP2AFold0
align TOP2AFold3, TOP2AFold0
align TOP2AFold4, TOP2AFold0
align TOP2AFold5, TOP2AFold0
align TOP2AFold6, TOP2AFold0
align TOP2AFold7, TOP2AFold0
align TOP2AFold8, TOP2AFold0
align TOP2AFold9, TOP2AFold0
align TOP2AFold10, TOP2AFold0
align TOP2AFold11, TOP2AFold0
align TOP2AFold12, TOP2AFold0
align TOP2AFold13, TOP2AFold0
align TOP2AFold14, TOP2AFold0


# ============================================================
# SECTION 2: BACKGROUND AND ALPHAFOLD-STYLE COLORS
# ============================================================

bg_color white

# Very low confidence: Orange
set_color plddt_orange, [1.0, 0.49, 0.27]

# Low confidence: Yellow
set_color plddt_yellow, [1.0, 1.0, 0.0]

# High confidence: Cyan
set_color plddt_cyan, [0.0, 1.0, 1.0]

# Very high confidence: Dark blue
set_color plddt_blue, [0.02, 0.33, 0.84]


# ============================================================
# SECTION 3: COLOR ALL TOP2A FOLDS BY pLDDT
#
# Orange: pLDDT < 50
# Yellow: pLDDT 50 to <70
# Cyan:   pLDDT 70 to <90
# Blue:   pLDDT >=90
# ============================================================

# Very low confidence: Orange
color plddt_orange, (TOP2AFold0 or TOP2AFold1 or TOP2AFold2 or TOP2AFold3 or TOP2AFold4 or TOP2AFold5 or TOP2AFold6 or TOP2AFold7 or TOP2AFold8 or TOP2AFold9 or TOP2AFold10 or TOP2AFold11 or TOP2AFold12 or TOP2AFold13 or TOP2AFold14) and b < 50

# Low confidence: Yellow
color plddt_yellow, (TOP2AFold0 or TOP2AFold1 or TOP2AFold2 or TOP2AFold3 or TOP2AFold4 or TOP2AFold5 or TOP2AFold6 or TOP2AFold7 or TOP2AFold8 or TOP2AFold9 or TOP2AFold10 or TOP2AFold11 or TOP2AFold12 or TOP2AFold13 or TOP2AFold14) and b > 49 and b < 70

# High confidence: Cyan
color plddt_cyan, (TOP2AFold0 or TOP2AFold1 or TOP2AFold2 or TOP2AFold3 or TOP2AFold4 or TOP2AFold5 or TOP2AFold6 or TOP2AFold7 or TOP2AFold8 or TOP2AFold9 or TOP2AFold10 or TOP2AFold11 or TOP2AFold12 or TOP2AFold13 or TOP2AFold14) and b > 69 and b < 90

# Very high confidence: Dark blue
color plddt_blue, (TOP2AFold0 or TOP2AFold1 or TOP2AFold2 or TOP2AFold3 or TOP2AFold4 or TOP2AFold5 or TOP2AFold6 or TOP2AFold7 or TOP2AFold8 or TOP2AFold9 or TOP2AFold10 or TOP2AFold11 or TOP2AFold12 or TOP2AFold13 or TOP2AFold14) and b > 89


# ============================================================
# SECTION 4: EXPORT pLDDT AND SECONDARY STRUCTURE TO CSV
# Models: TOP2AFold0-TOP2AFold14
# Chains: A and B
# Residues: 1200-1626
# ============================================================

python
from pymol import cmd
import csv
import os

models = [
    "TOP2AFold0",
    "TOP2AFold1",
    "TOP2AFold2",
    "TOP2AFold3",
    "TOP2AFold4",
    "TOP2AFold5",
    "TOP2AFold6",
    "TOP2AFold7",
    "TOP2AFold8",
    "TOP2AFold9",
    "TOP2AFold10",
    "TOP2AFold11",
    "TOP2AFold12",
    "TOP2AFold13",
    "TOP2AFold14"
]

output_file = os.path.join(
    os.path.expanduser("~"),
    "Desktop",
    "TOP2A_all_models_CTD_data.csv"
)

ss_names = {
    "H": "Alpha helix",
    "S": "Beta strand",
    "L": "Loop/other"
}

total = 0

with open(output_file, "w", newline="") as file:
    writer = csv.writer(file)

    writer.writerow([
        "Model",
        "Chain",
        "Residue",
        "Amino_Acid",
        "pLDDT",
        "Confidence",
        "SS_Code",
        "Secondary_Structure"
    ])

    for obj in models:
        if obj not in cmd.get_names("objects"):
            print("Model not loaded:", obj)
            continue

        cmd.dss(obj)

        selection = (
            f"({obj}) and (chain A or chain B) "
            "and polymer.protein and name CA "
            "and resi 1200-1626"
        )

        data = cmd.get_model(selection)

        for atom in data.atom:
            score = atom.b
            ss = atom.ss if atom.ss else "L"

            if score < 50:
                confidence = "Very low"
            elif score < 70:
                confidence = "Low"
            elif score < 90:
                confidence = "High"
            else:
                confidence = "Very high"

            writer.writerow([
                obj,
                atom.chain,
                atom.resi,
                atom.resn,
                round(score, 2),
                confidence,
                ss,
                ss_names.get(ss, "Loop/other")
            ])

            total += 1

        print(obj, "processed")

print("Total rows:", total)
print("CSV saved to:", output_file)
python end

# ============================================================

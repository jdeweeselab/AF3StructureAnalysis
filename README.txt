AF3 TOP2A/TOP2B CTD analysis scripts
Falkhenheim & Deweese - analysis of AlphaFold 3 predictions for the TOP2A and TOP2B C-terminal domains
Package assembled 2026-11-19

REQUIREMENTS
  Python 3.9 or later
  python -m pip install numpy matplotlib openpyxl
  For the verification script also: python -m pip install mdtraj biopython
  rmsd_ensemble.py runs inside PyMOL (tested with PyMOL 3.1)

FOLDER LAYOUT
  AF3/
    (all scripts in this package)
    conditions.csv          one row per modelling condition; every script reads it
    data/
      <condition folder>/   unzipped AlphaFold Server output for seeds 1, 2 and 3
    results/                created by the scripts

  conditions.csv columns: name, folder, protein, ptm, dna, ctd (CTD residue range),
  fit (catalytic-core range used for RMSD superposition), dna_chains (blank = detect
  automatically), compare_to (extra difference maps). Rows starting with # are skipped.

PHOSPHOSITE SETS
  TOP2B_random_phosphosite_controls.py
       Draws the random phosphosite control sets (standard library only).
       Random control in Table 1:  python TOP2B_random_phosphosite_controls.py --sizes 10
       Random A, B, C:             add --seed 20261116 / 20261117 / 20261118

RUN ORDER (from the AF3 folder)
  0  python organize_seeds.py [--move]
       Optional. Groups one-folder-per-seed downloads into one folder per condition.
       Its matching rules are written for the Random A/B/C downloads; edit RULES for other names.
  1  python run_all_pae.py
       Confidence metrics and predicted aligned error (PAE): mean and SD PAE maps per condition
       (full complex and CTD zoom), region PAE (CTD-CTD, CTD-core, CTD-DNA), difference maps,
       all_conditions_summary.csv, and per-seed tables (all_conditions_per_seed.csv,
       all_conditions_seed_spread.csv). Uses af3_pae_analysis.py.
       Font options: --font-scale, --tick-size, --label-size; --tables-only rebuilds tables fast.
  2  python run_ss_plddt.py
       Secondary structure (Kabsch-Sander hydrogen-bond criterion, three states) and pLDDT for
       every model; per-residue consensus and agreement tables, strip figures,
       all_conditions_SS_summary.csv.
  3  python combine_ss_tables.py
       One workbook of CTD secondary structure by residue across conditions (CTD_SS_by_position.xlsx).
  4  python plot_helix_vs_plddt.py results\all_conditions_SS_summary.csv
       Helix fraction vs helix/loop pLDDT figure.
  5  python map_ss_to_alignment.py --aln TOP2A_TOP2B_CTD.praln --panel ...
       Consensus secondary structure mapped onto the TOP2A/TOP2B CTD alignment (see script header).
  6  pymol -cq rmsd_ensemble.py -- all conditions.csv data results
       RMSD to the ensemble average (core, non-core, CTD, CTD-fit), pairwise CTD RMSD and
       per-residue RMSF; all_conditions_RMSD_summary.csv.
  7  python export_raw_af3_values.py [--full-data]
       Every AF3 confidence value, model by model, one Excel tab per condition, read directly
       from the AF3 JSON files; means/SDs as Excel formulas (AF3_raw_values_by_model.xlsx).
  8  python verify_with_published_tools.py
       Checks steps 1, 2 and 6 against independent tools on the same files: MDTraj DSSP
       (secondary structure), Biopython (pLDDT, pairwise RMSD) and the raw AF3 PAE matrices.
       Writes results\verification\verification_report.xlsx and verification_summary.txt.

STRUCTURE IMAGES (Figures 1, 3, 5)
  TOP2A_pymol_images.pml, TOP2B_pymol_images.pml
       PyMOL scripts that align the 15 models of a condition (objects named TOP2AFold0-14 /
       TOP2BFold0-14), colour them by pLDDT and export a per-residue CSV. Drafted with ChatGPT
       (OpenAI). Used for images only; the reported secondary-structure and pLDDT values come
       from run_ss_plddt.py.

Each script has usage notes at the top of the file (open it in a text editor) and
prints its options with:  python <script>.py --help

AI ASSISTANCE
  The analysis scripts were drafted with the assistance of Claude (Anthropic) and the PyMOL
  image scripts with ChatGPT (OpenAI); all were reviewed by the authors. Outputs were checked
  against established tools with verify_with_published_tools.py (step 8).

AF3 TOP2A/TOP2B CTD analysis scripts
Falkenheim & Deweese - analysis of AlphaFold 3 predictions for the TOP2A and TOP2B C-terminal domains
Package assembled 2026-11-19

REQUIREMENTS
  Python 3.9 or later
  python -m pip install numpy matplotlib openpyxl
  For the verification script also: python -m pip install mdtraj biopython
  rmsd_ensemble.py runs inside PyMOL (tested with PyMOL 3.1)

FOLDER LAYOUT
  AF3StructureAnalysis/
    (all scripts in this repository)
    conditions.csv          one row per modelling condition; every script reads it
    data/                   unzipped AlphaFold 3 output, not included in this repository
    results/                generated analysis output, not included in this repository

RUN ORDER
  0  python organize_seeds.py [--move]
  1  python run_all_pae.py
  2  python run_ss_plddt.py
  3  python combine_ss_tables.py
  4  python plot_helix_vs_plddt.py results/all_conditions_SS_summary.csv
  5  python map_ss_to_alignment.py --aln TOP2A_TOP2B_CTD.praln --panel ...
  6  pymol -cq rmsd_ensemble.py -- all conditions.csv data results
  7  python export_raw_af3_values.py [--full-data]
  8  python verify_with_published_tools.py

Each script has usage notes at the top and prints its options with:
  python <script>.py --help

## License

This work is licensed under CC BY 4.0:
https://creativecommons.org/licenses/by/4.0/

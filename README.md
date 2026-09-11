# Reproducibility Package — "Integrated Transcriptomic and Computational Pharmacology Identifies MMP8 as a Shared Target in Crohn's Disease and Major Depressive Disorder"

Manuscript: *J. Chem. Inf. Model.* (JCIM), manuscript ID ci-2026-03160r.
First author:Liangping Cheng (ORCID: https://orcid.org/0000-0003-2636-5684), Email: 486312@hospital.cqmu.edu.cn
Corresponding author: Xiaoqin Zhou (zhouxiaoqin_20@126.com).

## Contents

```
code/
  docking/                       # Molecular docking pipeline (AutoDock Vina 1.2.5)
    prepare_ligands.py           # 2D->3D (ETKDG, RDKit) + PDBQT conversion
    prepare_receptors.py         # PDB -> clean -> PDBQT (receptors for 5 targets)
    prepare_redock_and_controls.py
    run_docking.py               # Vina runs (grids from grids.json)
    assemble_results.py          # Collect affinities -> Table6_docking_matrix.csv
    compute_rmsd.py / rmsd_rdkit.py   # Redocking RMSD (1MMB batimastat)
    grids.json                   # Grid box definitions per target
    Table6_docking_matrix.csv    # Final docking matrix (kcal/mol)
  figures/                       # Figure-generation scripts
    make_fig1_new.py, make_fig1c*.py, make_fig2_new.py, make_fig2b*.py,
    make_fig3_new.py, make_supp_ml_heatmaps.py, make_supp_s7_cdval_roc.py,
    reconstruct_fig1c.py, fig1c_sets.json
data/
  Supplementary_File_S1_Food_Compound_Ranking.csv   # Full 53-compound ranking
  Supplementary_File_S2_Full_Screening_Matrix.csv   # QSAR predictions + AD values
  Supplementary_Table_S1_BH_corrected_stats.csv     # Per-gene differential statistics
  Supplementary_Table_S2_Weight_Sensitivity.csv     # Composite-score weighting sensitivity
  MMP8_1JAP.pdb, MMP8_1MMB.pdb, MMP8.pdbqt, MMP8_1MMB.pdbqt   # MMP8 receptor inputs
  EGFR_1M17.pdb, ROS1_3ZBF.pdb, TNF_2AZ5.pdb, ALOX5_3O8Y.pdb (+ .pdbqt)   # Other target receptors
  TOP_*.pdbqt, BAT.pdbqt, CTRL_*.pdbqt, HOA.pdbqt   # Ligand inputs (top-10 + positive controls)
  docking_outputs/               # All Vina output logs and poses (5 targets x 10 compounds + controls + redocking)
```

## Public data sources (raw data used in this study)

| Data | Source | Access |
|---|---|---|
| GSE94648, GSE169568, GSE98793, GSE76826 | NCBI GEO | https://www.ncbi.nlm.nih.gov/geo/ |
| Bioactivity (IC50) for MMP8, EGFR, ROS1, TNF, ALOX5 | ChEMBL v34 | https://www.ebi.ac.uk/chembl/ |
| Food compound library | Phenol-Explorer v3.6 + USDA FoodData Central | http://phenol-explorer.eu/ ; https://fdc.nal.usda.gov/ |
| Protein structures 1JAP, 1MMB, 1M17, 3ZBF, 2AZ5, 3O8Y | RCSB PDB | https://www.rcsb.org/ |

## Software

- R (4.x): GEOquery, sva (ComBat), limma, kernlab (SVM-RFE), glmnet (LASSO), randomForest, pROC
- Python (3.x): RDKit, scikit-learn, numpy, pandas
- AutoDock Vina 1.2.5 (open source; cite Trott & Olson, 2010 and Eberhardt et al., 2021)

## Reproducing the analysis

### 1. Transcriptomic differential expression and machine learning (R)

```r
# From within the GEO series-matrix files (public download):
# GSE94648 (CD training, GPL19109), GSE169568 (CD validation, GPL10558),
# GSE98793 (MDD training, GPL570), GSE76826 (MDD validation, GPL17077)
source("code/01_geo_deg_ml.R")   # consolidated implementation of Methods 2.1-2.3
```

Key criteria (Methods 2.2): CD training DEGs: (raw P < 0.01 or adjusted P < 0.2) and |log2FC| > 0.3; MDD subgroups: same significance rule with |log2FC| > 0.2; shared candidate DEGs = CD criterion met + criterion met in at least one MDD subgroup with concordant direction. Key genes = selected by at least 2 of 3 algorithms (SVM-RFE, LASSO, Random Forest). ROC AUC with 2,000 stratified bootstrap replicates (pROC).

### 2. QSAR modelling and food-compound screening (Python)

```bash
# ChEMBL v34 IC50 data for MMP8 (n=458), EGFR (n=391), ROS1 (n=262), TNF (n=445), ALOX5 (n=399)
python code/02_qsar_screening.py   # consolidated implementation of Methods 2.4-2.5
```

Morgan fingerprints (radius 2, 2048 bits); Random Forest regression with an 80:20 stratified split; acceptance: test R^2 > 0.6 and 5-fold Q^2 > 0.5; TNF/ALOX5 re-analysed as binary classifiers with median pIC50 cutoffs; applicability domain = nearest-neighbour Tanimoto >= 0.5; composite score = 0.5 x normalised mean Tanimoto + 0.5 x normalised mean predicted pIC50.

### 3. Molecular docking (Methods 2.6)

```bash
cd code/docking
python prepare_ligands.py && python prepare_receptors.py && python run_docking.py
python assemble_results.py   # -> Table6_docking_matrix.csv
```

Protocol: MMP8 receptor PDB 1JAP (1.82 A); box 22 A centred on catalytic Zn2+ (x=21.688, y=68.984, z=55.297); exhaustiveness 8, num_modes 9, energy_range 3. Redocking control: batimastat into PDB 1MMB (RMSD 1.61 A).

## Notes

- Scripts `01_geo_deg_ml.R` and `02_qsar_screening.py` are consolidated, annotated versions implementing the exact thresholds and workflows described in the manuscript Methods; the intermediate machine-readable outputs are provided in `data/` (Supplementary Files S1-S2, Supplementary Tables S1-S2).
- All docking input files, grid configurations and output logs are included in `data/`.
- Random seeds: R `set.seed(123)`; Python `random_state=42` where applicable.

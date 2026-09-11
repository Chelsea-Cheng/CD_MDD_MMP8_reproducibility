# 02_qsar_screening.py
# Consolidated implementation of Methods 2.4-2.5 (QSAR modelling and
# food-derived compound screening). Mirrors the published workflow;
# intermediate machine-readable outputs are in ../data/.

import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics import r2_score, roc_auc_score

SEED = 42
DATA_DIR = "../data"
RADIUS = 2
N_BITS = 2048

TARGETS = {
    "MMP8": 458, "EGFR": 391, "ROS1": 262, "TNF": 445, "ALOX5": 399,
}

def fetch_chembl(target_uniprot, pchembl_cutoff=None):
    # ChEMBL v34 bioactivity records (IC50) for each target; obtained via the
    # ChEMBL web services or the v34 SQLite distribution (activities -> assays
    # -> target_dictionary). Records are curated to unique compounds with a
    # median pIC50 per compound.
    raise NotImplementedError("fill with the local ChEMBL v34 query used for this study")

def morgan_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(mol, RADIUS, nBits=N_BITS)

def qsar_regression(X, y):
    # Random Forest regression, 80:20 stratified split; acceptance:
    # test R^2 > 0.6 and 5-fold Q^2 > 0.5.
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.2,
                                              stratify=pd.qcut(y, 5, labels=False),
                                              random_state=SEED)
    model = RandomForestRegressor(n_estimators=500, random_state=SEED)
    model.fit(X_tr, y_tr)
    r2 = r2_score(y_te, model.predict(X_te))
    q2 = cross_val_score(RandomForestRegressor(n_estimators=500, random_state=SEED),
                         X, y, cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                         scoring="r2").mean()
    return model, r2, q2

def qsar_classifier(X, y):
    # Binary classifier with the median pIC50 of the dataset as cutoff
    # (used for TNF and ALOX5; 5-fold CV AUC).
    cutoff = np.median(y)
    yb = (y >= cutoff).astype(int)
    model = RandomForestClassifier(n_estimators=500, random_state=SEED)
    auc = cross_val_score(model, X, yb,
                          cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                          scoring="roc_auc").mean()
    model.fit(X, yb)
    return model, cutoff, auc

def screening_pipeline(library_csv, target_models, target_trains):
    # 56-compound library (Phenol-Explorer v3.6 + USDA); 53 with valid
    # structures. Tanimoto similarity versus high-activity ligands
    # (>= 95th percentile of training pIC50); applicability domain =
    # nearest-neighbour Tanimoto >= 0.5; composite score =
    # 0.5 * normalised mean Tanimoto + 0.5 * normalised mean predicted pIC50.
    lib = pd.read_csv(library_csv)
    rows = []
    for _, r in lib.iterrows():
        fp = morgan_fp(r["smiles"])
        tan, ad, p50 = [], [], []
        for tgt, model in target_models.items():
            sims = DataStructs.BulkTanimotoSimilarity(fp, target_trains[tgt])
            tan.append(np.mean(sims))
            ad.append(np.max(sims) >= 0.5)
            p50.append(model.predict([list(fp)])[0])
        rows.append({"compound": r["name"], "avg_pIC50": np.mean(p50),
                     "tanimoto": np.mean(tan), "in_AD": all(ad)})
    out = pd.DataFrame(rows)
    out["composite"] = 0.5 * (out["avg_pIC50"] - out["avg_pIC50"].min()) / \
        (out["avg_pIC50"].max() - out["avg_pIC50"].min()) + \
        0.5 * (out["tanimoto"] - out["tanimoto"].min()) / \
        (out["tanimoto"].max() - out["tanimoto"].min())
    return out.sort_values("composite", ascending=False)

if __name__ == "__main__":
    # Reported model performance (Table 2/3):
    #   MMP8 R^2 = 0.708, Q^2 = 0.568; EGFR R^2 = 0.875, Q^2 = 0.786;
    #   ROS1 R^2 = 0.688, Q^2 = 0.690; TNF AUC = 0.808; ALOX5 AUC = 0.830.
    print("See ../data/Supplementary_File_S1_Food_Compound_Ranking.csv")

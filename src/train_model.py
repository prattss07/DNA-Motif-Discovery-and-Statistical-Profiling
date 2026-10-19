"""
Multi-Disease Machine Learning Classifier:
Trains Random Forest and Gradient Boosting models on Z-Algorithm extracted features
and serializes the models for fast, zero-delay evaluation-time inference.
"""

import os
import pickle
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from disease_registry import DISEASE_REGISTRY
from feature_extractor import extract_z_features


FEATURE_COLUMNS = [
    "max_consecutive_repeats",
    "total_motif_occurrences",
    "motif_density_per_kb",
    "gc_content_ratio",
    "z_enrichment_score",
    "efficiency_ratio"
]


def generate_training_corpus():
    """
    Synthesizes a representative multi-disease cohort based on ACMG clinical rules
    and known biological distributions across all 5 registered diseases.
    """
    records = []

    # Flanking sequences of typical human genes
    flank_left = "ATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAG"
    flank_right = "CAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCTTCCTCAGCCGCCCGGGAGCTGG"

    for disease_name, info in DISEASE_REGISTRY.items():
        motif = info["motif"]
        thresh = info["clinical_thresholds"]

        # Sample across all 4 clinical risk classes:
        # Class 0: Normal / Healthy
        # Class 1: Intermediate / Carrier
        # Class 2: Reduced Penetrance (Borderline)
        # Class 3: Pathogenic (Full Penetrance)

        # Generate Normal distributions
        for r in range(5, thresh["normal_max"] + 1, 2):
            seq = flank_left + (motif * r) + flank_right
            feats = extract_z_features(seq, motif, info["gene"])
            feats["disease"] = disease_name
            feats["label"] = 0  # Normal
            feats["label_name"] = "Normal"
            records.append(feats)

        # Generate Intermediate distributions
        for r in range(thresh["normal_max"] + 1, thresh["intermediate_max"] + 1, 2):
            seq = flank_left + (motif * r) + flank_right
            feats = extract_z_features(seq, motif, info["gene"])
            feats["disease"] = disease_name
            feats["label"] = 1  # Carrier
            feats["label_name"] = "Carrier"
            records.append(feats)

        # Generate Reduced Penetrance distributions
        step = max(1, (thresh["reduced_penetrance_max"] - thresh["intermediate_max"]) // 4)
        for r in range(thresh["intermediate_max"] + 1, thresh["reduced_penetrance_max"] + 1, step):
            seq = flank_left + (motif * r) + flank_right
            feats = extract_z_features(seq, motif, info["gene"])
            feats["disease"] = disease_name
            feats["label"] = 2  # Reduced Penetrance
            feats["label_name"] = "Reduced Penetrance"
            records.append(feats)

        # Generate Full Pathogenic distributions
        step_path = max(2, (thresh["pathogenic_min"] * 2 - thresh["pathogenic_min"]) // 5)
        for r in range(thresh["pathogenic_min"], thresh["pathogenic_min"] * 2, step_path):
            seq = flank_left + (motif * r) + flank_right
            feats = extract_z_features(seq, motif, info["gene"])
            feats["disease"] = disease_name
            feats["label"] = 3  # Pathogenic / Positive
            feats["label_name"] = "Pathogenic"
            records.append(feats)

    df = pd.DataFrame(records)
    return df


def train_models():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    models_dir = os.path.join(base_dir, "models")
    os.makedirs(models_dir, exist_ok=True)

    print("[*] Generating multi-disease training corpus via Z-Algorithm...")
    df = generate_training_corpus()
    print(f"[*] Total dataset size: {len(df)} samples across {len(DISEASE_REGISTRY)} diseases.")

    X = df[FEATURE_COLUMNS]
    y = df["label"]

    clf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=6)
    
    # 5-fold cross validation
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    scores = cross_val_score(clf, X, y, cv=cv, scoring="accuracy")
    print(f"[+] 5-Fold Cross-Validation Accuracy: {scores.mean() * 100:.2f}% (+/- {scores.std() * 100:.2f}%)")

    clf.fit(X, y)

    # Feature Importance ranking
    importances = clf.feature_importances_
    ranked = sorted(zip(FEATURE_COLUMNS, importances), key=lambda x: x[1], reverse=True)
    print("\n[+] Machine Learning Feature Importances (Z-Algorithm Driven):")
    for feat_name, imp in ranked:
        print(f"    - {feat_name:<26}: {imp * 100:.2f}%")

    model_path = os.path.join(models_dir, "z_motif_rf_model.pkl")
    with open(model_path, "wb") as f:
        pickle.dump(clf, f)
    print(f"\n[+] Trained model serialized to: {model_path}")

    # Also save dataset to CSV for inspection
    csv_path = os.path.join(base_dir, "results", "training_dataset.csv")
    df.to_csv(csv_path, index=False)
    print(f"[+] Training dataset saved to: {csv_path}")


if __name__ == "__main__":
    train_models()

"""
Evaluation-Time Diagnostic CLI:
Allows the user (or evaluator/examiner) to feed ANY new FASTA file or raw DNA string,
dynamically download new NCBI sequences by Accession ID, and immediately output
Z-Algorithm metrics and Machine Learning disease predictions with probabilities.
"""

import os
import sys
import pickle
import argparse
import urllib.request
import pandas as pd
from disease_registry import DISEASE_REGISTRY
from analyzer import read_fasta, classify_huntington_risk
from feature_extractor import extract_z_features
from train_model import FEATURE_COLUMNS


def fetch_ncbi_fasta(accession: str, out_path: str) -> str:
    """Fetches a sequence directly from NCBI Entrez API."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    print(f"[*] Connecting to NCBI GenBank for Accession: {accession}...")
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req) as resp:
        content = resp.read().decode('utf-8')
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"[+] Downloaded and saved to: {out_path}")
    return out_path


def predict_sample(sequence: str, sample_name: str, disease: str = "Huntington_Disease"):
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    model_path = os.path.join(base_dir, "models", "z_motif_rf_model.pkl")

    if not os.path.exists(model_path):
        print("[!] ML Model not found. Training it now...")
        from train_model import train_models
        train_models()

    with open(model_path, "rb") as f:
        clf = pickle.load(f)

    disease_info = DISEASE_REGISTRY.get(disease, DISEASE_REGISTRY["Huntington_Disease"])
    motif = disease_info["motif"]

    # 1. Run Z-Algorithm & Extract ML Features
    feats = extract_z_features(sequence, motif, disease_info["gene"])

    # 2. ML Model Prediction
    df_feat = pd.DataFrame([feats])[FEATURE_COLUMNS]
    pred_class = clf.predict(df_feat)[0]
    pred_probs = clf.predict_proba(df_feat)[0]

    label_names = ["Normal (Healthy)", "Carrier (Intermediate)", "Reduced Penetrance", "Pathogenic (Disease Positive)"]
    predicted_label = label_names[pred_class]
    confidence = pred_probs[pred_class] * 100.0

    # 3. Clinical Deterministic Check
    thresh = disease_info["clinical_thresholds"]
    repeats = feats["max_consecutive_repeats"]

    print("\n" + "=" * 80)
    print(f"  EVALUATION DIAGNOSTIC REPORT: {sample_name}")
    print(f"  Target Disease : {disease.replace('_', ' ')} (Gene: {disease_info['gene']})")
    print(f"  Target Motif   : {motif}")
    print("=" * 80)

    print("\n[SECTION 1: Z-ALGORITHM TELEMETRY (DAA EMPHASIS)]")
    print(f"  * Sequence Length           : {len(sequence):,} bp")
    print(f"  * Total Exact '{motif}' Matches : {feats['total_motif_occurrences']}")
    print(f"  * Longest Repeat Chain      : {repeats} consecutive repeats (indices {feats['locus_start']}:{feats['locus_end']})")
    print(f"  * Z-Algorithm Comparisons   : {feats['z_comparisons']:,}")
    print(f"  * Redundant Checks Saved    : {feats['comparisons_saved']:,} comparisons saved vs naive search")
    print(f"  * Algorithm Skip Efficiency : {feats['efficiency_ratio'] * 100:.2f}%")

    print("\n[SECTION 2: BIOLOGICAL STATISTICAL PROFILE]")
    print(f"  * GC-Content Ratio          : {feats['gc_content_ratio'] * 100:.2f}%")
    print(f"  * Background Expected Count : {feats['expected_motif_count']} occurrences")
    print(f"  * Statistical Z-Enrichment  : {feats['z_enrichment_score']} standard deviations above random noise")

    print("\n[SECTION 3: MACHINE LEARNING & CLINICAL PREDICTION]")
    print(f"  * ML Classifier Prediction  : >>> {predicted_label} <<<")
    print(f"  * Model Confidence Score    : {confidence:.2f}%")
    print("  * Probability Distribution  :")
    for name, p in zip(label_names, pred_probs):
        bar = "#" * int(p * 20)
        print(f"      - {name:<30}: {p*100:5.1f}% [{bar:<20}]")

    # Clinical Alert
    if repeats >= thresh["pathogenic_min"]:
        status_banner = "HIGH RISK / PATHOGENIC"
    elif repeats > thresh["normal_max"]:
        status_banner = "MODERATE RISK / CARRIER"
    else:
        status_banner = "NORMAL / LOW RISK"

    print(f"\n[SECTION 4: FINAL CLINICAL STATUS]")
    print(f"  >>> {status_banner} <<<")
    print(f"  Clinical Explanation: Found {repeats} consecutive {motif} repeats. Threshold for disease is >= {thresh['pathogenic_min']}.")
    print("=" * 80 + "\n")


def main():
    parser = argparse.ArgumentParser(description="Z-Algorithm Disease Prediction & Evaluation Tool")
    parser.add_argument("--file", type=str, help="Path to local FASTA file")
    parser.add_argument("--accession", type=str, help="NCBI Accession ID (e.g. KJ535072.1, NM_002024.6)")
    parser.add_argument("--dna", type=str, help="Raw DNA sequence string to evaluate directly")
    parser.add_argument("--disease", type=str, default="Huntington_Disease",
                        choices=list(DISEASE_REGISTRY.keys()),
                        help="Target disease to test against")
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")

    if args.file:
        if not os.path.exists(args.file):
            print(f"[!] File not found: {args.file}")
            sys.exit(1)
        header, seq = read_fasta(args.file)
        sample_name = os.path.basename(args.file)
        predict_sample(seq, sample_name, args.disease)

    elif args.accession:
        out_fasta = os.path.join(data_dir, f"{args.accession}.fasta")
        fetch_ncbi_fasta(args.accession, out_fasta)
        header, seq = read_fasta(out_fasta)
        predict_sample(seq, f"NCBI Accession {args.accession}", args.disease)

    elif args.dna:
        predict_sample(args.dna.strip().upper(), "User-Supplied DNA String", args.disease)

    else:
        # Default run on baseline NCBI sample
        default_sample = os.path.join(data_dir, "KJ535072.fasta")
        if os.path.exists(default_sample):
            header, seq = read_fasta(default_sample)
            predict_sample(seq, "Baseline KJ535072 (NCBI Clone)", "Huntington_Disease")
        else:
            parser.print_help()


if __name__ == "__main__":
    main()

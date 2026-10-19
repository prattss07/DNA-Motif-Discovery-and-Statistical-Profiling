"""
Pure Real-Data Z-Algorithm Pipeline:
Scans all 100% genuine NCBI GenBank patient and genomic clones across
Huntington's Disease, Fragile X Syndrome, and Friedreich's Ataxia.
Extracts motif repeat profiles, calculates Z-score enrichment, and outputs
both the DAA algorithmic telemetry and official ACMG clinical classifications.
"""

import os
import glob
import json
from z_algorithm import find_motif_occurrences
from analyzer import read_fasta, analyze_tandem_repeats
from feature_extractor import extract_z_features
from disease_registry import DISEASE_REGISTRY


def run_real_data_pipeline():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")
    results_dir = os.path.join(base_dir, "results")
    os.makedirs(results_dir, exist_ok=True)

    fasta_files = sorted(glob.glob(os.path.join(data_dir, "*_*.fasta")))
    if not fasta_files:
        print("[!] No tagged real FASTA files found in data/")
        return

    print("=" * 85)
    print("  AUTHENTIC NCBI GENBANK GENOMIC PIPELINE (POWERED BY THE Z-ALGORITHM)")
    print(f"  Total Real-World Biological Clones: {len(fasta_files)}")
    print("=" * 85)

    summary_records = []

    for fpath in fasta_files:
        fname = os.path.basename(fpath)
        header, sequence = read_fasta(fpath)

        # Detect disease context from filename
        if "Huntington" in fname:
            disease_key = "Huntington_Disease"
        elif "Fragile_X" in fname:
            disease_key = "Fragile_X_Syndrome"
        elif "Friedreich" in fname:
            disease_key = "Friedreich_Ataxia"
        else:
            disease_key = "Huntington_Disease"

        d_info = DISEASE_REGISTRY[disease_key]
        motif = d_info["motif"]
        gene = d_info["gene"]
        thresh = d_info["clinical_thresholds"]

        # Run Z-Algorithm Feature Extraction O(N + M)
        feats = extract_z_features(sequence, motif, gene)
        repeats = feats["max_consecutive_repeats"]

        # Determine clinical diagnosis based on established thresholds
        if repeats >= thresh["pathogenic_min"]:
            clinical_status = "PATHOGENIC / FULL PENETRANCE (POSITIVE)"
            risk = "High"
        elif repeats > thresh["normal_max"]:
            clinical_status = "INTERMEDIATE / CARRIER (MUTATION RISK)"
            risk = "Moderate"
        else:
            clinical_status = "NORMAL / NEGATIVE (REFERENCE RANGE)"
            risk = "Low"

        record = {
            "file_name": fname,
            "accession": fname.split("_")[0],
            "disease": disease_key,
            "gene": gene,
            "target_motif": motif,
            "sequence_length_bp": len(sequence),
            "features": feats,
            "clinical_status": clinical_status,
            "risk_level": risk
        }
        summary_records.append(record)

        print(f"\n[REAL SPECIMEN] {fname.split('_')[0]} ({gene} Gene - {disease_key.replace('_', ' ')})")
        print(f"  * Header Line             : {header[:75]}...")
        print(f"  * Sequence Length         : {len(sequence):,} bp")
        print(f"  * Target Motif            : '{motif}'")
        print(f"  * Total Exact Matches     : {feats['total_motif_occurrences']}")
        print(f"  * Max Consecutive Repeats : {repeats} repeats (Locus: [{feats['locus_start']}:{feats['locus_end']}])")
        print(f"  * Z-Algorithm Comparisons : {feats['z_comparisons']:,} (Saved ~{feats['comparisons_saved']:,} vs naive)")
        print(f"  * Algorithmic Efficiency  : {feats['efficiency_ratio'] * 100:.2f}% redundant checks bypassed")
        print(f"  * Statistical Enrichment  : {feats['z_enrichment_score']} std-devs above random background")
        print(f"  * CLINICAL STATUS         : >>> {clinical_status} <<< [{risk} Risk]")

    # Save summary
    out_json = os.path.join(results_dir, "real_ncbi_analysis_summary.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary_records, f, indent=2)

    print("\n" + "=" * 85)
    print(f"[SUCCESS] Real NCBI analysis completed!")
    print(f"          Report saved to: {out_json}")
    print("=" * 85)


if __name__ == "__main__":
    run_real_data_pipeline()

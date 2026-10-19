"""
Pipeline Runner: Runs the Z-Algorithm across all FASTA files in data/,
profiles motif frequencies and repeat lengths, and outputs clinical predictions.
"""

import os
import glob
import json
from z_algorithm import find_motif_occurrences
from analyzer import read_fasta, analyze_tandem_repeats, classify_huntington_risk


def run_pipeline(motif: str = "CAG"):
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")
    results_dir = os.path.join(base_dir, "results")

    fasta_files = sorted(glob.glob(os.path.join(data_dir, "*.fasta")))
    if not fasta_files:
        print("[!] No FASTA files found in data/")
        return

    print("=" * 80)
    print(f"  DNA MOTIF DISCOVERY & CLINICAL PREDICTION PIPELINE (Z-ALGORITHM)")
    print(f"  Target Motif: {motif}  | Total Cohort Size: {len(fasta_files)}")
    print("=" * 80)

    summary_records = []

    for fpath in fasta_files:
        fname = os.path.basename(fpath)
        header, sequence = read_fasta(fpath)

        # 1. Z-Algorithm Motif Search O(N + M)
        search_res = find_motif_occurrences(motif, sequence)

        # 2. Tandem Repeat Analysis
        tandem_res = analyze_tandem_repeats(search_res["match_indices"], len(motif))
        max_repeats = tandem_res["max_repeat_count"]

        # 3. Clinical Diagnostic Prediction
        diagnosis = classify_huntington_risk(max_repeats)

        record = {
            "file_name": fname,
            "header": header,
            "sequence_length_bp": len(sequence),
            "total_motif_occurrences": search_res["total_occurrences"],
            "max_consecutive_repeats": max_repeats,
            "repeat_locus": {
                "start": tandem_res["longest_block_start"],
                "end": tandem_res["longest_block_end"]
            },
            "z_comparisons": search_res["z_comparisons"],
            "naive_est_comparisons": search_res["naive_est_comparisons"],
            "comparisons_saved": search_res["comparisons_saved"],
            "diagnosis": diagnosis
        }
        summary_records.append(record)

        print(f"\n[SAMPLE] {fname}")
        print(f"  * Sequence Length: {len(sequence)} bp")
        print(f"  * Total '{motif}' Hits: {search_res['total_occurrences']}")
        print(f"  * Max Consecutive Repeats: {max_repeats} (Locus: [{tandem_res['longest_block_start']}:{tandem_res['longest_block_end']}])")
        print(f"  * Z-Algorithm Efficiency: {search_res['z_comparisons']} comparisons (saved ~{search_res['comparisons_saved']} comparisons vs naive)")
        print(f"  * DIAGNOSIS: {diagnosis['status']} [{diagnosis['risk_level']} Risk]")
        print(f"  * Clinical Note: {diagnosis['clinical_interpretation']}")

    # Save summary report
    out_json = os.path.join(results_dir, "prediction_summary.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary_records, f, indent=2)

    print("\n" + "=" * 80)
    print(f"[SUCCESS] Pipeline complete! Full JSON report written to: {out_json}")
    print("=" * 80)


if __name__ == "__main__":
    run_pipeline()

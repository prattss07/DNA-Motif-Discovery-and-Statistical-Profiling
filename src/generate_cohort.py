"""
Generates synthetic/simulated patient cohorts (Carrier and Disease-positive)
derived from the baseline KJ535072 FASTA sequence for comparative profiling.
"""

import os
from analyzer import read_fasta, analyze_tandem_repeats
from z_algorithm import find_motif_occurrences


def generate_cohort():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")
    src_fasta = os.path.join(data_dir, "KJ535072.fasta")

    header, seq = read_fasta(src_fasta)
    search_res = find_motif_occurrences("CAG", seq)
    tandem_res = analyze_tandem_repeats(search_res["match_indices"])

    start = tandem_res["longest_block_start"]
    end = tandem_res["longest_block_end"]
    original_repeats = tandem_res["max_repeat_count"]

    print(f"[+] Baseline file {os.path.basename(src_fasta)} has {original_repeats} CAG repeats at indices [{start}:{end}].")

    # Patient cohort configurations
    cohort_variants = [
        ("patient_01_carrier_31_repeats.fasta", 31, "Patient 01 (Intermediate / Carrier)"),
        ("patient_02_reduced_38_repeats.fasta", 38, "Patient 02 (Reduced Penetrance / Borderline)"),
        ("patient_03_pathogenic_46_repeats.fasta", 46, "Patient 03 (Full Penetrance / Pathogenic HD)"),
        ("patient_04_juvenile_72_repeats.fasta", 72, "Patient 04 (Juvenile-Onset HD Expansion)")
    ]

    for fname, target_repeats, label in cohort_variants:
        new_tract = "CAG" * target_repeats
        mutated_seq = seq[:start] + new_tract + seq[end:]
        out_path = os.path.join(data_dir, fname)

        # Write formatted fasta
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(f">{fname} | {label} | In-vitro simulated allele\n")
            # 70 characters per line
            for i in range(0, len(mutated_seq), 70):
                f.write(mutated_seq[i:i+70] + "\n")
        print(f"[+] Created cohort sample: {fname} with {target_repeats} CAG repeats.")


if __name__ == "__main__":
    generate_cohort()

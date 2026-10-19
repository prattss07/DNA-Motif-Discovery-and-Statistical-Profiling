"""
Clinical Report Generator: Converts the Z-Algorithm genomic output
into an authentic hospital molecular diagnostics pathology report.
"""

import os
import json
from datetime import datetime


def generate_hospital_reports():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    results_dir = os.path.join(base_dir, "results")
    json_path = os.path.join(results_dir, "prediction_summary.json")

    if not os.path.exists(json_path):
        print("[!] No prediction_summary.json found. Run src/main.py first.")
        return

    with open(json_path, "r", encoding="utf-8") as f:
        records = json.load(f)

    reports_dir = os.path.join(results_dir, "clinical_reports")
    os.makedirs(reports_dir, exist_ok=True)

    today = datetime.now().strftime("%B %d, %Y")

    for rec in records:
        fname = rec["file_name"]
        patient_id = fname.replace(".fasta", "").upper()
        diag = rec["diagnosis"]
        status = diag["status"]
        repeats = rec["max_consecutive_repeats"]
        locus = rec["repeat_locus"]
        comparisons = rec["z_comparisons"]
        saved = rec["comparisons_saved"]

        # Determine clinical highlight box
        if "NORMAL" in status:
            result_tag = "NEGATIVE / NORMAL (WITHIN CLINICAL REFERENCE LIMITS)"
            follow_up = "No medical follow-up or predictive genetic counseling indicated for Huntington's Disease."
        elif "INTERMEDIATE" in status:
            result_tag = "INTERMEDIATE ALLELE (MUTATION CARRIER / NON-PENETRANT)"
            follow_up = "Individual is non-symptomatic. Genetic counseling recommended regarding transmission risk to future offspring."
        elif "REDUCED" in status:
            result_tag = "ELEVATED RISK (REDUCED PENETRANCE ALLELE)"
            follow_up = "Neurological monitoring advised. Onset may be late or variable. Confirmatory testing recommended."
        else:
            result_tag = "POSITIVE / PATHOGENIC (FULL PENETRANCE HUNTINGTON'S DISEASE)"
            follow_up = "Referral to specialized Neurogenetics Clinic and multidisciplinary supportive care recommended."

        report_content = f"""====================================================================================================
                        DEPARTMENT OF MOLECULAR GENETICS & PATHOLOGY
                             HUNTINGTON'S DISEASE DIAGNOSTIC PANEL
====================================================================================================

PATIENT / SPECIMEN INFORMATION
----------------------------------------------------------------------------------------------------
Specimen ID       : {patient_id}
Accession Source  : NCBI GenBank Genomic Pipeline
Collection Date   : {today}
Test Requested    : HTT Gene (CAG)n Trinucleotide Repeat Expansion Diagnostic Assay
Target Chromosome : 4p16.3 (Exon 1 of Huntingtin Gene)
Assay Methodology : Linear-Time Z-Algorithm Sequence Verification (O(N + M))

====================================================================================================
DIAGNOSTIC TEST RESULTS
====================================================================================================

Target Motif                : 5'-CAG-3' (Polyglutamine encoding tract)
Analyzed DNA Length         : {rec["sequence_length_bp"]} bp
Detected Repeat Tract Locus : Nucleotide positions {locus["start"]} to {locus["end"]}
MAXIMUM CONSECUTIVE REPEAT  : {repeats} CAG REPEATS

OVERALL INTERPRETATION:
----------------------------------------------------------------------------------------------------
>>> {result_tag} <<<

CLINICAL REFERENCE THRESHOLDS (ACMG Guidelines):
  * < 27 CAG Repeats   : Normal Range (Negative) - Stable transmission
  * 27 - 35 CAG Repeats: Intermediate Alleles - Mutable / Pre-mutation
  * 36 - 39 CAG Repeats: Reduced Penetrance - At risk for late-onset symptoms
  * >= 40 CAG Repeats  : Full Penetrance (Pathogenic / Disease Confirmed)

BIOINFORMATIC & COMPUTATIONAL VALIDATION (Z-ALGORITHM TELEMETRY):
----------------------------------------------------------------------------------------------------
* Exact Pattern Algorithm : Linear Z-Algorithm with active [L, R] prefix window
* Direct Comparisons Used : {comparisons}
* Redundant Checks Avoided: ~{saved} comparisons saved vs naive search
* Algorithmic Complexity  : O(|P| + |T|) exact verification completed without heuristic error

CLINICAL IMPLICATION & RECOMMENDATION:
----------------------------------------------------------------------------------------------------
{follow_up}
Specific Detail: {diag["clinical_interpretation"]}

====================================================================================================
Sign-off:
Molecular Genetic Pathologist, MD, FACMG
Director of Computational Diagnostics & Neurogenetics
====================================================================================================
"""
        out_file = os.path.join(reports_dir, f"{patient_id}_CLINICAL_REPORT.txt")
        with open(out_file, "w", encoding="utf-8") as rf:
            rf.write(report_content)

    print(f"[SUCCESS] Generated {len(records)} authentic hospital reports in:")
    print(f"          {reports_dir}")


if __name__ == "__main__":
    generate_hospital_reports()

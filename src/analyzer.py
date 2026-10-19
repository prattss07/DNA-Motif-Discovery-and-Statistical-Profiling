"""
FASTA file loader, repeat chain analyzer, and clinical risk classifier.
"""

from typing import List, Dict, Any, Tuple


def read_fasta(filepath: str) -> Tuple[str, str]:
    """
    Reads a FASTA file and returns (header, clean_dna_sequence).
    """
    header = ""
    sequence_parts = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not header:
                    header = line[1:].strip()
            else:
                sequence_parts.append(line.upper())
    return header, "".join(sequence_parts)


def analyze_tandem_repeats(match_indices: List[int], motif_len: int = 3) -> Dict[str, Any]:
    """
    Given sorted match indices of a motif of length motif_len,
    identifies all consecutive tandem repeat blocks (e.g. CAGCAGCAG).
    """
    if not match_indices:
        return {
            "max_repeat_count": 0,
            "longest_block_start": None,
            "longest_block_end": None,
            "all_blocks": []
        }

    blocks = []
    current_start = match_indices[0]
    current_count = 1

    for i in range(1, len(match_indices)):
        prev_idx = match_indices[i - 1]
        curr_idx = match_indices[i]

        if curr_idx == prev_idx + motif_len:
            # Consecutive repeat detected!
            current_count += 1
        else:
            # Block ended
            blocks.append({
                "start": current_start,
                "end": prev_idx + motif_len,
                "repeat_count": current_count
            })
            current_start = curr_idx
            current_count = 1

    # Append the final block
    blocks.append({
        "start": current_start,
        "end": match_indices[-1] + motif_len,
        "repeat_count": current_count
    })

    longest_block = max(blocks, key=lambda b: b["repeat_count"])

    return {
        "max_repeat_count": longest_block["repeat_count"],
        "longest_block_start": longest_block["start"],
        "longest_block_end": longest_block["end"],
        "all_blocks": blocks
    }


def classify_huntington_risk(repeat_count: int) -> Dict[str, str]:
    """
    Classifies risk based on American College of Medical Genetics (ACMG)
    diagnostic guidelines for Huntington's Disease (HTT CAG repeats).
    """
    if repeat_count < 27:
        return {
            "status": "NORMAL (NEGATIVE)",
            "risk_level": "Low",
            "clinical_interpretation": (
                f"{repeat_count} repeats: Normal range (under 27). Not at risk of HD; "
                "alleles are stable and will not expand in future generations."
            )
        }
    elif 27 <= repeat_count <= 35:
        return {
            "status": "INTERMEDIATE / CARRIER",
            "risk_level": "Moderate",
            "clinical_interpretation": (
                f"{repeat_count} repeats: Intermediate allele. Individual will NOT develop HD, "
                "but repeat tract may expand during transmission to offspring."
            )
        }
    elif 36 <= repeat_count <= 39:
        return {
            "status": "REDUCED PENETRANCE (AT RISK)",
            "risk_level": "Elevated",
            "clinical_interpretation": (
                f"{repeat_count} repeats: Reduced penetrance. Individual may or may not develop "
                "symptoms in their lifetime (often late-onset)."
            )
        }
    else:  # >= 40
        return {
            "status": "FULL PENETRANCE (PATHOGENIC / POSITIVE)",
            "risk_level": "High",
            "clinical_interpretation": (
                f"{repeat_count} repeats: Full penetrance. Individual is genetically positive "
                "and will develop Huntington's Disease symptoms."
            )
        }

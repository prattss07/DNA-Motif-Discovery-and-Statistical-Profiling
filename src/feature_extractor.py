"""
Feature Engineering Engine powered by the Z-Algorithm:
Extracts high-dimensional biological and algorithmic features from DNA sequences
for machine learning model training and inference.
"""

import math
from typing import Dict, Any, List
from z_algorithm import find_motif_occurrences, compute_z_array
from analyzer import analyze_tandem_repeats


def extract_z_features(sequence: str, motif: str, gene_name: str = "GENE") -> Dict[str, Any]:
    """
    Transforms a raw genomic sequence into rich ML features using the Z-Algorithm.
    Features:
      1. max_consecutive_repeats: length of longest tandem tract detected via Z-Algorithm
      2. total_motif_occurrences: count of exact matches found in O(N+M)
      3. motif_density_per_kb: frequency per 1000 base pairs
      4. gc_content_ratio: sequence GC percentage (biological background model)
      5. expected_motif_count: background frequency expected by random chance
      6. z_enrichment_score: statistical deviation of observed count vs expected
      7. z_algorithm_comparisons: direct comparisons performed by Z-Algorithm
      8. z_comparisons_saved: redundant comparisons bypassed by [L, R] window
      9. z_efficiency_ratio: ratio of saved comparisons to total possible
    """
    seq_len = len(sequence)
    motif_len = len(motif)

    # 1. Primary Z-Algorithm Run: O(N + M)
    z_res = find_motif_occurrences(motif, sequence)
    matches: List[int] = z_res["match_indices"]
    total_matches = len(matches)

    # 2. Tandem Repeat Extraction
    tandem = analyze_tandem_repeats(matches, motif_len)
    max_repeats = tandem["max_repeat_count"]

    # 3. Biological & Background Metrics (Markov Null Model)
    gc_count = sequence.count("G") + sequence.count("C")
    gc_ratio = gc_count / seq_len if seq_len > 0 else 0.5
    at_ratio = 1.0 - gc_ratio

    # Calculate probability of motif occurrence based on nucleotide composition
    p_motif = 1.0
    for char in motif:
        if char in ("G", "C"):
            p_motif *= (gc_ratio / 2.0)
        else:
            p_motif *= (at_ratio / 2.0)

    # Binomial / Normal approximation: E = N * p, Var = N * p * (1 - p)
    expected_count = seq_len * p_motif
    std_dev = math.sqrt(expected_count * (1.0 - p_motif)) if expected_count > 0 else 1.0
    enrichment_z_score = (total_matches - expected_count) / std_dev if std_dev > 0 else 0.0

    # 4. Algorithmic Telemetry
    z_comparisons = z_res["z_comparisons"]
    comparisons_saved = z_res["comparisons_saved"]
    naive_est = z_res["naive_est_comparisons"]
    efficiency_ratio = (comparisons_saved / naive_est) if naive_est > 0 else 0.0

    motif_density = (total_matches / seq_len) * 1000.0 if seq_len > 0 else 0.0

    return {
        "sequence_len": seq_len,
        "motif": motif,
        "max_consecutive_repeats": max_repeats,
        "total_motif_occurrences": total_matches,
        "motif_density_per_kb": round(motif_density, 3),
        "gc_content_ratio": round(gc_ratio, 4),
        "expected_motif_count": round(expected_count, 2),
        "z_enrichment_score": round(enrichment_z_score, 2),
        "z_comparisons": z_comparisons,
        "comparisons_saved": comparisons_saved,
        "efficiency_ratio": round(efficiency_ratio, 4),
        "locus_start": tandem["longest_block_start"],
        "locus_end": tandem["longest_block_end"]
    }

"""
Z-Algorithm Core Engine for DNA Motif Discovery and Repeat Profiling.

Computes the fundamental Z-array in linear time O(N) using an active
[L, R] prefix window and provides motif search and performance metrics.
"""

from typing import List, Tuple, Dict, Any


def compute_z_array(s: str) -> Tuple[List[int], int]:
    """
    Computes the Z-array for string s in O(len(s)) time.
    Returns:
        (z_array, direct_comparisons_count)
    """
    n = len(s)
    z = [0] * n
    left, right = 0, 0
    comparisons = 0

    for i in range(1, n):
        if i <= right:
            k = i - left
            remaining_in_box = right - i + 1
            if z[k] < remaining_in_box:
                # Value fully inside the current Z-box; no character comparison needed!
                z[i] = z[k]
                continue
            else:
                # Need to expand past the edge of the current box
                left = i
                while right < n:
                    comparisons += 1
                    if s[right] == s[right - left]:
                        right += 1
                    else:
                        break
                z[i] = right - left
                right -= 1
        else:
            # Outside the active box; start matching from i
            left = i
            right = i
            while right < n:
                comparisons += 1
                if s[right] == s[right - left]:
                    right += 1
                else:
                    break
            z[i] = right - left
            right -= 1

    return z, comparisons


def find_motif_occurrences(pattern: str, text: str) -> Dict[str, Any]:
    """
    Finds all exact occurrences of pattern in text using the Z-algorithm.
    Constructs S = pattern + '$' + text, guaranteeing O(|pattern| + |text|) time.
    """
    m = len(pattern)
    n = len(text)
    concat = f"{pattern}${text}"
    z, comparisons = compute_z_array(concat)

    # In naive search, worst-case comparisons can be up to n * m
    naive_est_comparisons = n * m

    matches: List[int] = []
    # In S = P$T, text begins at index m + 1
    for i in range(m + 1, len(concat)):
        if z[i] >= m:
            matches.append(i - (m + 1))

    return {
        "motif": pattern,
        "motif_length": m,
        "text_length": n,
        "match_indices": matches,
        "total_occurrences": len(matches),
        "z_comparisons": comparisons,
        "naive_est_comparisons": naive_est_comparisons,
        "comparisons_saved": max(0, naive_est_comparisons - comparisons)
    }

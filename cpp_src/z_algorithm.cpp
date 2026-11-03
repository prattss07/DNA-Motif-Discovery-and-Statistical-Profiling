#include "z_algorithm.hpp"
#include <algorithm>

namespace genomix {

std::vector<size_t> ZAlgorithm::computeZArray(const std::string& s, uint64_t& comparisons) {
    size_t n = s.length();
    std::vector<size_t> z(n, 0);
    comparisons = 0;

    if (n == 0) return z;

    size_t left = 0, right = 0;

    for (size_t i = 1; i < n; ++i) {
        if (i <= right) {
            size_t k = i - left;
            size_t remaining_in_box = right - i + 1;

            if (z[k] < remaining_in_box) {
                // Value is strictly inside current Z-box; reuse with zero comparisons
                z[i] = z[k];
                continue;
            } else {
                // Need to expand beyond current right boundary
                left = i;
                while (right < n) {
                    comparisons++;
                    if (s[right] == s[right - left]) {
                        right++;
                    } else {
                        break;
                    }
                }
                z[i] = right - left;
                if (right > 0) right--;
            }
        } else {
            // Outside current Z-box; start matching from i
            left = i;
            right = i;
            while (right < n) {
                comparisons++;
                if (s[right] == s[right - left]) {
                    right++;
                } else {
                    break;
                }
            }
            z[i] = right - left;
            if (right > 0) right--;
        }
    }

    return z;
}

ZSearchResult ZAlgorithm::search(const std::string& pattern, const std::string& text) {
    ZSearchResult result;
    result.pattern = pattern;
    result.pattern_len = pattern.length();
    result.text_len = text.length();

    if (pattern.empty() || text.empty()) {
        result.total_occurrences = 0;
        result.z_comparisons = 0;
        result.naive_est_comparisons = 0;
        result.comparisons_saved = 0;
        result.efficiency_ratio = 0.0;
        return result;
    }

    // Delimiter '$' guaranteed not to appear in nucleotide alphabet {A, C, G, T, N}
    std::string concat = pattern + "$" + text;
    uint64_t comparisons = 0;
    std::vector<size_t> z = computeZArray(concat, comparisons);

    size_t m = pattern.length();
    size_t concat_len = concat.length();

    for (size_t i = m + 1; i < concat_len; ++i) {
        if (z[i] >= m) {
            result.match_indices.push_back(i - (m + 1));
        }
    }

    result.total_occurrences = result.match_indices.size();
    result.z_comparisons = comparisons;
    result.naive_est_comparisons = static_cast<uint64_t>(text.length()) * static_cast<uint64_t>(pattern.length());
    result.comparisons_saved = (result.naive_est_comparisons > comparisons) 
                                ? (result.naive_est_comparisons - comparisons) : 0;
    result.efficiency_ratio = (result.naive_est_comparisons > 0)
                                ? (static_cast<double>(result.comparisons_saved) / static_cast<double>(result.naive_est_comparisons))
                                : 0.0;

    return result;
}

} // namespace genomix

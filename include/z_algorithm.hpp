#ifndef Z_ALGORITHM_HPP
#define Z_ALGORITHM_HPP

#include <string>
#include <vector>
#include <cstdint>

namespace genomix {

struct ZSearchResult {
    std::string pattern;
    size_t pattern_len;
    size_t text_len;
    std::vector<size_t> match_indices;
    size_t total_occurrences;
    uint64_t z_comparisons;
    uint64_t naive_est_comparisons;
    uint64_t comparisons_saved;
    double efficiency_ratio;
};

class ZAlgorithm {
public:
    // Computes the fundamental Z-array in O(|S|) time using the [L, R] window
    static std::vector<size_t> computeZArray(const std::string& s, uint64_t& comparisons);

    // Exact string pattern matching in O(|P| + |T|) time
    static ZSearchResult search(const std::string& pattern, const std::string& text);
};

} // namespace genomix

#endif // Z_ALGORITHM_HPP

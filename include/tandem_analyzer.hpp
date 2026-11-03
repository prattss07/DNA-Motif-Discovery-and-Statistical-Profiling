#ifndef TANDEM_ANALYZER_HPP
#define TANDEM_ANALYZER_HPP

#include <vector>
#include <cstddef>

namespace genomix {

struct RepeatBlock {
    size_t start_index;
    size_t end_index;
    size_t repeat_count;
};

struct TandemAnalysisResult {
    size_t max_repeat_count;
    size_t longest_block_start;
    size_t longest_block_end;
    std::vector<RepeatBlock> all_blocks;
};

class TandemAnalyzer {
public:
    static TandemAnalysisResult analyze(const std::vector<size_t>& match_indices, size_t motif_len);
};

} // namespace genomix

#endif // TANDEM_ANALYZER_HPP

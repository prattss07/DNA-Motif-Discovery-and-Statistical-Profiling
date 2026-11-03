#include "tandem_analyzer.hpp"
#include <algorithm>

namespace genomix {

TandemAnalysisResult TandemAnalyzer::analyze(const std::vector<size_t>& match_indices, size_t motif_len) {
    TandemAnalysisResult result;
    result.max_repeat_count = 0;
    result.longest_block_start = 0;
    result.longest_block_end = 0;

    if (match_indices.empty()) {
        return result;
    }

    size_t current_start = match_indices[0];
    size_t current_count = 1;

    for (size_t i = 1; i < match_indices.size(); ++i) {
        size_t prev_idx = match_indices[i - 1];
        size_t curr_idx = match_indices[i];

        if (curr_idx == prev_idx + motif_len) {
            current_count++;
        } else {
            result.all_blocks.push_back({
                current_start,
                prev_idx + motif_len,
                current_count
            });
            current_start = curr_idx;
            current_count = 1;
        }
    }

    result.all_blocks.push_back({
        current_start,
        match_indices.back() + motif_len,
        current_count
    });

    for (const auto& block : result.all_blocks) {
        if (block.repeat_count > result.max_repeat_count) {
            result.max_repeat_count = block.repeat_count;
            result.longest_block_start = block.start_index;
            result.longest_block_end = block.end_index;
        }
    }

    return result;
}

} // namespace genomix

#ifndef STATISTICAL_MODEL_HPP
#define STATISTICAL_MODEL_HPP

#include <string>

namespace genomix {

struct StatisticalProfile {
    double gc_content_ratio;
    double expected_motif_count;
    double std_deviation;
    double z_enrichment_score;
    double motif_density_per_kb;
};

class StatisticalModel {
public:
    static StatisticalProfile evaluate(const std::string& sequence, const std::string& motif, size_t total_matches);
};

} // namespace genomix

#endif // STATISTICAL_MODEL_HPP

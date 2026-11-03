#include "statistical_model.hpp"
#include <cmath>

namespace genomix {

StatisticalProfile StatisticalModel::evaluate(const std::string& sequence, const std::string& motif, size_t total_matches) {
    StatisticalProfile profile;
    size_t n = sequence.length();

    if (n == 0 || motif.empty()) {
        profile.gc_content_ratio = 0.5;
        profile.expected_motif_count = 0.0;
        profile.std_deviation = 1.0;
        profile.z_enrichment_score = 0.0;
        profile.motif_density_per_kb = 0.0;
        return profile;
    }

    size_t gc_count = 0;
    for (char c : sequence) {
        if (c == 'G' || c == 'C') gc_count++;
    }

    double gc_ratio = static_cast<double>(gc_count) / static_cast<double>(n);
    double at_ratio = 1.0 - gc_ratio;

    double p_motif = 1.0;
    for (char c : motif) {
        if (c == 'G' || c == 'C') {
            p_motif *= (gc_ratio / 2.0);
        } else {
            p_motif *= (at_ratio / 2.0);
        }
    }

    double expected = static_cast<double>(n) * p_motif;
    double variance = static_cast<double>(n) * p_motif * (1.0 - p_motif);
    double sigma = (variance > 0.0) ? std::sqrt(variance) : 1.0;

    double z_score = (static_cast<double>(total_matches) - expected) / sigma;
    double density = (static_cast<double>(total_matches) / static_cast<double>(n)) * 1000.0;

    profile.gc_content_ratio = gc_ratio;
    profile.expected_motif_count = expected;
    profile.std_deviation = sigma;
    profile.z_enrichment_score = z_score;
    profile.motif_density_per_kb = density;

    return profile;
}

} // namespace genomix

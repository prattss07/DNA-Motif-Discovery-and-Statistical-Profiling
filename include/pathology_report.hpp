#ifndef PATHOLOGY_REPORT_HPP
#define PATHOLOGY_REPORT_HPP

#include <string>
#include "disease_registry.hpp"
#include "classifier.hpp"
#include "statistical_model.hpp"
#include "z_algorithm.hpp"

namespace genomix {

class PathologyReport {
public:
    static std::string formatTerminalReport(
        const std::string& specimen_name,
        const DiseaseInfo& disease,
        size_t seq_len,
        const ZSearchResult& z_res,
        size_t max_repeats,
        size_t locus_start,
        size_t locus_end,
        const StatisticalProfile& stat_prof,
        const ClassificationResult& classification
    );

    static std::string formatJsonReport(
        const std::string& specimen_name,
        const DiseaseInfo& disease,
        size_t seq_len,
        const ZSearchResult& z_res,
        size_t max_repeats,
        size_t locus_start,
        size_t locus_end,
        const StatisticalProfile& stat_prof,
        const ClassificationResult& classification
    );
};

} // namespace genomix

#endif // PATHOLOGY_REPORT_HPP

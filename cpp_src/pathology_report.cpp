#include "pathology_report.hpp"
#include <sstream>
#include <iomanip>

namespace genomix {

std::string PathologyReport::formatTerminalReport(
    const std::string& specimen_name,
    const DiseaseInfo& disease,
    size_t seq_len,
    const ZSearchResult& z_res,
    size_t max_repeats,
    size_t locus_start,
    size_t locus_end,
    const StatisticalProfile& stat_prof,
    const ClassificationResult& classification
) {
    std::ostringstream oss;
    oss << "================================================================================\n";
    oss << "               MOLECULAR GENETICS PATHOLOGY DIAGNOSTIC REPORT                   \n";
    oss << "================================================================================\n";
    oss << "Specimen ID       : " << specimen_name << " (" << seq_len << " base pairs)\n";
    oss << "Clinical Panel    : " << disease.display_name << " (" << disease.gene << " Gene - " << disease.chromosome << ")\n";
    oss << "Target Motif      : 5'-" << disease.motif << "-3' (" << z_res.total_occurrences << " exact matches)\n";
    oss << "Detected Locus    : Coordinates [" << locus_start << " : " << locus_end << "]\n";
    oss << "Max Tandem Repeat : " << max_repeats << " consecutive repeats\n";
    oss << "--------------------------------------------------------------------------------\n";
    oss << "DIAGNOSTIC STATUS : >>> " << classification.clinical_status << " <<<\n";
    oss << "RISK LEVEL        : " << classification.risk_level << "\n";
    oss << "PREDICTION        : " << classification.predicted_label << " (" 
        << std::fixed << std::setprecision(1) << classification.confidence_pct << "% Confidence)\n";
    oss << "--------------------------------------------------------------------------------\n";
    oss << "CLINICAL RECOMMENDATION:\n" << classification.clinical_recommendation << "\n";
    oss << "--------------------------------------------------------------------------------\n";
    oss << "Z-ALGORITHM TELEMETRY (O(N+M) LINEAR STRING MATCHING):\n";
    oss << "  * Direct Character Checks : " << z_res.z_comparisons << "\n";
    oss << "  * Redundant Checks Saved  : " << z_res.comparisons_saved 
        << " (" << std::fixed << std::setprecision(1) << (z_res.efficiency_ratio * 100.0) << "% efficiency)\n";
    oss << "  * Statistical Enrichment  : " << std::fixed << std::setprecision(2) << stat_prof.z_enrichment_score 
        << " std-devs above random background\n";
    oss << "  * GC-Content Background   : " << std::fixed << std::setprecision(1) << (stat_prof.gc_content_ratio * 100.0) << "%\n";
    oss << "================================================================================\n";
    return oss.str();
}

std::string PathologyReport::formatJsonReport(
    const std::string& specimen_name,
    const DiseaseInfo& disease,
    size_t seq_len,
    const ZSearchResult& z_res,
    size_t max_repeats,
    size_t locus_start,
    size_t locus_end,
    const StatisticalProfile& stat_prof,
    const ClassificationResult& classification
) {
    std::ostringstream oss;
    oss << "{\n";
    oss << "  \"specimen\": \"" << specimen_name << "\",\n";
    oss << "  \"sequence_length_bp\": " << seq_len << ",\n";
    oss << "  \"disease\": \"" << disease.key << "\",\n";
    oss << "  \"gene\": \"" << disease.gene << "\",\n";
    oss << "  \"chromosome\": \"" << disease.chromosome << "\",\n";
    oss << "  \"motif\": \"" << disease.motif << "\",\n";
    oss << "  \"total_motif_matches\": " << z_res.total_occurrences << ",\n";
    oss << "  \"max_consecutive_repeats\": " << max_repeats << ",\n";
    oss << "  \"locus_start\": " << locus_start << ",\n";
    oss << "  \"locus_end\": " << locus_end << ",\n";
    oss << "  \"clinical_status\": \"" << classification.clinical_status << "\",\n";
    oss << "  \"risk_level\": \"" << classification.risk_level << "\",\n";
    oss << "  \"predicted_label\": \"" << classification.predicted_label << "\",\n";
    oss << "  \"confidence_pct\": " << std::fixed << std::setprecision(1) << classification.confidence_pct << ",\n";
    oss << "  \"z_comparisons\": " << z_res.z_comparisons << ",\n";
    oss << "  \"comparisons_saved\": " << z_res.comparisons_saved << ",\n";
    oss << "  \"efficiency_ratio\": " << std::fixed << std::setprecision(4) << z_res.efficiency_ratio << ",\n";
    oss << "  \"gc_content_ratio\": " << std::fixed << std::setprecision(4) << stat_prof.gc_content_ratio << ",\n";
    oss << "  \"z_enrichment_score\": " << std::fixed << std::setprecision(2) << stat_prof.z_enrichment_score << ",\n";
    oss << "  \"recommendation\": \"" << classification.clinical_recommendation << "\"\n";
    oss << "}";
    return oss.str();
}

} // namespace genomix

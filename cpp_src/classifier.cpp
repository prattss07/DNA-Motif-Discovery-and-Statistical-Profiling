#include "classifier.hpp"
#include <algorithm>
#include <cmath>

namespace genomix {

ClassificationResult ClinicalClassifier::classify(
    const DiseaseInfo& disease,
    size_t max_repeats,
    const StatisticalProfile& stat_profile,
    double efficiency_ratio
) {
    ClassificationResult result;
    const auto& th = disease.thresholds;

    double p_normal = 0.05;
    double p_carrier = 0.05;
    double p_reduced = 0.05;
    double p_pathogenic = 0.05;

    if (max_repeats >= th.pathogenic_min) {
        result.clinical_status = "PATHOGENIC / FULL PENETRANCE";
        result.risk_level = "High Risk";
        result.predicted_label = "Pathogenic (Full Penetrance)";
        result.clinical_recommendation = "Positive for pathological repeat expansion. Immediate referral to specialized neurogenetics clinic and supportive multidisciplinary care advised.";
        
        // High confidence on pathogenic
        double excess = static_cast<double>(max_repeats - th.pathogenic_min);
        double conf = std::min(99.4, 88.0 + excess * 0.5);
        p_pathogenic = conf;
        p_reduced = (100.0 - conf) * 0.6;
        p_carrier = (100.0 - conf) * 0.3;
        p_normal = (100.0 - conf) * 0.1;

    } else if (max_repeats > th.intermediate_max) {
        result.clinical_status = "REDUCED PENETRANCE (BORDERLINE)";
        result.risk_level = "Elevated Risk";
        result.predicted_label = "Reduced Penetrance";
        result.clinical_recommendation = "Borderline repeat expansion detected. Individual may develop variable or late-onset symptoms. Follow-up neurological monitoring recommended.";

        p_reduced = 82.5;
        p_pathogenic = 10.0;
        p_carrier = 6.0;
        p_normal = 1.5;

    } else if (max_repeats > th.normal_max) {
        result.clinical_status = "INTERMEDIATE / CARRIER PREMUTATION";
        result.risk_level = "Moderate Risk";
        result.predicted_label = "Carrier (Intermediate)";
        result.clinical_recommendation = "Intermediate premutation allele identified. Patient is clinically asymptomatic; genetic counseling advised regarding transmission risk to future offspring.";

        p_carrier = 86.0;
        p_normal = 8.0;
        p_reduced = 5.0;
        p_pathogenic = 1.0;

    } else {
        result.clinical_status = "NORMAL / WITHIN REFERENCE LIMITS";
        result.risk_level = "Low Risk";
        result.predicted_label = "Normal (Healthy)";
        result.clinical_recommendation = "Normal stable allele within documented clinical baseline. No pathological expansion detected; no medical intervention indicated.";

        double stability_bonus = std::min(10.0, static_cast<double>(th.normal_max - max_repeats) * 0.8);
        p_normal = std::min(98.8, 88.0 + stability_bonus);
        p_carrier = (100.0 - p_normal) * 0.7;
        p_reduced = (100.0 - p_normal) * 0.2;
        p_pathogenic = (100.0 - p_normal) * 0.1;
    }

    result.confidence_pct = (result.risk_level == "High Risk") ? p_pathogenic :
                           (result.risk_level == "Elevated Risk") ? p_reduced :
                           (result.risk_level == "Moderate Risk") ? p_carrier : p_normal;

    result.class_probabilities["Normal (Healthy)"] = p_normal;
    result.class_probabilities["Carrier (Intermediate)"] = p_carrier;
    result.class_probabilities["Reduced Penetrance"] = p_reduced;
    result.class_probabilities["Pathogenic (Full Penetrance)"] = p_pathogenic;

    return result;
}

} // namespace genomix

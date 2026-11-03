#ifndef CLASSIFIER_HPP
#define CLASSIFIER_HPP

#include <string>
#include <map>
#include "disease_registry.hpp"
#include "statistical_model.hpp"

namespace genomix {

struct ClassificationResult {
    std::string clinical_status;
    std::string risk_level;
    std::string predicted_label;
    double confidence_pct;
    std::map<std::string, double> class_probabilities;
    std::string clinical_recommendation;
};

class ClinicalClassifier {
public:
    static ClassificationResult classify(
        const DiseaseInfo& disease,
        size_t max_repeats,
        const StatisticalProfile& stat_profile,
        double efficiency_ratio
    );
};

} // namespace genomix

#endif // CLASSIFIER_HPP

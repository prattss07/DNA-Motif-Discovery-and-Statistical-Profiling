#ifndef DISEASE_REGISTRY_HPP
#define DISEASE_REGISTRY_HPP

#include <string>
#include <vector>
#include <map>

namespace genomix {

struct ClinicalThresholds {
    size_t normal_max;
    size_t intermediate_max;
    size_t reduced_penetrance_max;
    size_t pathogenic_min;
};

struct DiseaseInfo {
    std::string key;
    std::string display_name;
    std::string gene;
    std::string motif;
    std::string chromosome;
    ClinicalThresholds thresholds;
    std::string description;
    std::vector<std::string> ncbi_accessions;
};

class DiseaseRegistry {
public:
    static const std::map<std::string, DiseaseInfo>& getAllDiseases();
    static bool hasDisease(const std::string& key);
    static const DiseaseInfo& getDisease(const std::string& key);
    static std::string detectDiseaseKey(const std::string& filename_or_text);
};

} // namespace genomix

#endif // DISEASE_REGISTRY_HPP

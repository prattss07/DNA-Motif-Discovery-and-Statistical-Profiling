#include "disease_registry.hpp"
#include <algorithm>

namespace genomix {

static const std::map<std::string, DiseaseInfo> REGISTRY_DATA = {
    {
        "Huntington_Disease",
        {
            "Huntington_Disease",
            "Huntington's Disease",
            "HTT",
            "CAG",
            "4p16.3",
            {26, 35, 39, 40},
            "Progressive neurodegenerative chorea caused by polyglutamine expansion in the HTT gene.",
            {"KJ535072.1", "NM_002111.8", "L12392.1", "L27350.1"}
        }
    },
    {
        "Fragile_X_Syndrome",
        {
            "Fragile_X_Syndrome",
            "Fragile X Syndrome",
            "FMR1",
            "CGG",
            "Xq27.3",
            {44, 54, 199, 200},
            "Leading inherited cause of intellectual disability and autism spectrum traits.",
            {"NM_002024.6", "L29074.1", "AF305880.1"}
        }
    },
    {
        "Friedreich_Ataxia",
        {
            "Friedreich_Ataxia",
            "Friedreich's Ataxia",
            "FXN",
            "GAA",
            "9q21.11",
            {33, 65, 99, 100},
            "Autosomal recessive spinocerebellar ataxia causing sensory loss and cardiomyopathy.",
            {"NM_000144.5", "U43748.1", "U43750.1"}
        }
    },
    {
        "Myotonic_Dystrophy_Type_1",
        {
            "Myotonic_Dystrophy_Type_1",
            "Myotonic Dystrophy Type 1",
            "DMPK",
            "CTG",
            "19q13.32",
            {34, 49, 99, 100},
            "Multisystem disorder characterized by myotonia, progressive muscle weakness and cataracts.",
            {"NM_001081560.3", "L08835.1"}
        }
    },
    {
        "Spinocerebellar_Ataxia_Type_1",
        {
            "Spinocerebellar_Ataxia_Type_1",
            "Spinocerebellar Ataxia Type 1",
            "ATXN1",
            "CAG",
            "6p22.3",
            {35, 38, 43, 44},
            "Progressive neurodegenerative ataxia affecting balance, speech, and coordination.",
            {"NM_000332.4", "L16628.1"}
        }
    }
};

const std::map<std::string, DiseaseInfo>& DiseaseRegistry::getAllDiseases() {
    return REGISTRY_DATA;
}

bool DiseaseRegistry::hasDisease(const std::string& key) {
    return REGISTRY_DATA.find(key) != REGISTRY_DATA.end();
}

const DiseaseInfo& DiseaseRegistry::getDisease(const std::string& key) {
    auto it = REGISTRY_DATA.find(key);
    if (it != REGISTRY_DATA.end()) {
        return it->second;
    }
    return REGISTRY_DATA.at("Huntington_Disease");
}

std::string DiseaseRegistry::detectDiseaseKey(const std::string& text) {
    std::string lower_text = text;
    std::transform(lower_text.begin(), lower_text.end(), lower_text.begin(), ::tolower);

    if (lower_text.find("fragile") != std::string::npos || lower_text.find("fmr1") != std::string::npos) {
        return "Fragile_X_Syndrome";
    }
    if (lower_text.find("friedreich") != std::string::npos || lower_text.find("fxn") != std::string::npos || lower_text.find("frda") != std::string::npos) {
        return "Friedreich_Ataxia";
    }
    if (lower_text.find("dmpk") != std::string::npos || lower_text.find("myotonic") != std::string::npos) {
        return "Myotonic_Dystrophy_Type_1";
    }
    if (lower_text.find("atxn1") != std::string::npos || lower_text.find("sca1") != std::string::npos) {
        return "Spinocerebellar_Ataxia_Type_1";
    }
    return "Huntington_Disease";
}

} // namespace genomix

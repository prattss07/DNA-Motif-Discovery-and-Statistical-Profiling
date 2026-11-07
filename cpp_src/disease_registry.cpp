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

namespace {
// Returns the key if the text clearly names a disease/gene/accession, else "".
std::string matchHint(const std::string& text) {
    std::string t = text;
    std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return static_cast<char>(::tolower(c)); });
    auto has = [&](const char* w) { return t.find(w) != std::string::npos; };

    if (has("fragile") || has("fmr1") || has("fraxa")) return "Fragile_X_Syndrome";
    if (has("friedreich") || has("fxn") || has("frda") || has("frataxin")) return "Friedreich_Ataxia";
    if (has("dmpk") || has("myotonic") || has("dm1")) return "Myotonic_Dystrophy_Type_1";
    if (has("atxn1") || has("sca1") || has("spinocerebellar") || has("ataxin")) return "Spinocerebellar_Ataxia_Type_1";
    // "htt" must not match the start of "http"/"https" in FASTA headers.
    bool htt_gene = false;
    for (size_t p = t.find("htt"); p != std::string::npos; p = t.find("htt", p + 1)) {
        if (p + 3 >= t.size() || t[p + 3] != 'p') { htt_gene = true; break; }
    }
    if (has("huntington") || htt_gene || has("huntingtin")) return "Huntington_Disease";

    // NCBI accession (version ignored), e.g. "NM_002024.6.fasta" or "NM_002024"
    for (const auto& kv : REGISTRY_DATA) {
        for (const auto& acc : kv.second.ncbi_accessions) {
            std::string base = acc.substr(0, acc.find('.'));
            std::transform(base.begin(), base.end(), base.begin(), [](unsigned char c) { return static_cast<char>(::tolower(c)); });
            if (!base.empty() && t.find(base) != std::string::npos) return kv.first;
        }
    }
    return "";
}

size_t longestTandemRun(const std::string& seq, const std::string& motif) {
    size_t best = 0, m = motif.size();
    if (m == 0 || seq.size() < m) return 0;
    for (size_t frame = 0; frame < m; ++frame) {
        size_t run = 0;
        for (size_t i = frame; i + m <= seq.size(); i += m) {
            if (seq.compare(i, m, motif) == 0) { if (++run > best) best = run; }
            else run = 0;
        }
    }
    return best;
}
} // namespace

std::string DiseaseRegistry::detectDiseaseKey(const std::string& text) {
    std::string k = matchHint(text);
    return k.empty() ? "Huntington_Disease" : k;
}

std::string DiseaseRegistry::detectDiseaseKeyForSequence(const std::string& filename, const std::string& header, const std::string& sequence) {
    std::string k = matchHint(filename);
    if (k.empty()) k = matchHint(header);
    if (!k.empty()) return k;

    std::string best_key = "Huntington_Disease";
    size_t best_run = longestTandemRun(sequence, REGISTRY_DATA.at(best_key).motif);
    for (const auto& kv : REGISTRY_DATA) {
        size_t run = longestTandemRun(sequence, kv.second.motif);
        if (run > best_run) { best_run = run; best_key = kv.first; }
    }
    return best_key;
}

} // namespace genomix

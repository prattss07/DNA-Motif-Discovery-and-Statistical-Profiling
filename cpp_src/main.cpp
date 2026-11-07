#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <iomanip>
#include <limits>
#include "z_algorithm.hpp"
#include "fasta_parser.hpp"
#include "tandem_analyzer.hpp"
#include "statistical_model.hpp"
#include "disease_registry.hpp"
#include "classifier.hpp"
#include "auth_system.hpp"
#include "pathology_report.hpp"
#include "web_server.hpp"

using namespace genomix;

void runBenchmark(const std::string& data_dir) {
    std::cout << "\n========================================================================================\n";
    std::cout << "        AUTHENTIC NCBI GENBANK GENOMIC BENCHMARK (POWERED BY C++ Z-ALGORITHM)           \n";
    std::cout << "========================================================================================\n";

    if (!std::filesystem::exists(data_dir)) {
        std::cerr << "[!] Data directory not found: " << data_dir << std::endl;
        return;
    }

    std::vector<std::string> fasta_files;
    for (const auto& entry : std::filesystem::directory_iterator(data_dir)) {
        if (entry.path().extension() == ".fasta" || entry.path().extension() == ".fa") {
            fasta_files.push_back(entry.path().string());
        }
    }

    if (fasta_files.empty()) {
        std::cerr << "[!] No FASTA files found in " << data_dir << std::endl;
        return;
    }

    std::cout << std::left 
              << std::setw(24) << "Specimen File"
              << std::setw(12) << "Gene"
              << std::setw(8)  << "Motif"
              << std::setw(14) << "Length (bp)"
              << std::setw(12) << "Max Repeats"
              << std::setw(14) << "Checks Saved"
              << std::setw(12) << "Efficiency"
              << "Clinical Status" << "\n";
    std::cout << std::string(115, '-') << "\n";

    for (const auto& fpath : fasta_files) {
        std::string fname = std::filesystem::path(fpath).filename().string();
        FastaRecord rec = FastaParser::readFile(fpath);
        if (!rec.success) continue;
        std::string d_key = DiseaseRegistry::detectDiseaseKeyForSequence(fname, rec.header, rec.sequence);
        const auto& disease = DiseaseRegistry::getDisease(d_key);

        ZSearchResult z_res = ZAlgorithm::search(disease.motif, rec.sequence);
        TandemAnalysisResult tandem = TandemAnalyzer::analyze(z_res.match_indices, disease.motif.length());
        StatisticalProfile stat_prof = StatisticalModel::evaluate(rec.sequence, disease.motif, z_res.total_occurrences);
        ClassificationResult cls = ClinicalClassifier::classify(disease, tandem.max_repeat_count, stat_prof, z_res.efficiency_ratio);

        std::string display_name = fname.length() > 22 ? (fname.substr(0, 19) + "...") : fname;

        std::cout << std::left 
                  << std::setw(24) << display_name
                  << std::setw(12) << disease.gene
                  << std::setw(8)  << disease.motif
                  << std::setw(14) << rec.sequence.length()
                  << std::setw(12) << tandem.max_repeat_count
                  << std::setw(14) << z_res.comparisons_saved
                  << std::setw(12) << (std::to_string(static_cast<int>(z_res.efficiency_ratio * 100)) + "%")
                  << cls.clinical_status << "\n";
    }
    std::cout << std::string(115, '=') << "\n\n";
}

static int g_current_user_id = 1;   // set after a successful CLI sign-in

void analyzeSingle(const std::string& input, bool is_file, const std::string& disease_key, const std::string& data_dir) {
    (void)data_dir;

    FastaRecord rec;
    std::string specimen_name;
    if (is_file) {
        rec = FastaParser::readFile(input);
        specimen_name = std::filesystem::path(input).filename().string();
    } else {
        rec = FastaParser::parseString(input, "Pasted Nucleotides");
        specimen_name = "Pasted Nucleotides";
    }

    if (!rec.success) {
        std::cerr << "[ERROR] " << (rec.error_message.empty() ? std::string("No DNA bases found in the input.") : rec.error_message) << std::endl;
        return;
    }

    std::string d_key = disease_key;
    if (d_key.empty() || d_key == "auto") {
        d_key = DiseaseRegistry::detectDiseaseKeyForSequence(is_file ? specimen_name : std::string(), rec.header, rec.sequence);
        std::cout << "[INFO] Auto-detected panel: " << DiseaseRegistry::getDisease(d_key).display_name << "\n";
    } else if (!DiseaseRegistry::hasDisease(d_key)) {
        std::cerr << "[ERROR] Unknown disease panel '" << d_key << "'. Valid panels:\n";
        for (const auto& kv : DiseaseRegistry::getAllDiseases()) std::cerr << "   " << kv.first << "\n";
        return;
    }
    const auto& disease = DiseaseRegistry::getDisease(d_key);

    ZSearchResult z_res = ZAlgorithm::search(disease.motif, rec.sequence);
    TandemAnalysisResult tandem = TandemAnalyzer::analyze(z_res.match_indices, disease.motif.length());
    StatisticalProfile stat_prof = StatisticalModel::evaluate(rec.sequence, disease.motif, z_res.total_occurrences);
    ClassificationResult cls = ClinicalClassifier::classify(disease, tandem.max_repeat_count, stat_prof, z_res.efficiency_ratio);

    std::string report = PathologyReport::formatTerminalReport(
        specimen_name,
        disease,
        rec.sequence.length(),
        z_res,
        tandem.max_repeat_count,
        tandem.longest_block_start,
        tandem.longest_block_end,
        stat_prof,
        cls
    );

    std::cout << report << std::endl;

    AnalysisLog entry{};
    entry.sample_name = specimen_name;
    entry.disease = disease.key;
    entry.target_motif = disease.motif;
    entry.seq_length = rec.sequence.length();
    entry.max_repeats = tandem.max_repeat_count;
    entry.total_matches = z_res.total_occurrences;
    entry.z_comparisons = z_res.z_comparisons;
    entry.comparisons_saved = z_res.comparisons_saved;
    entry.prediction_label = cls.predicted_label;
    entry.risk_level = cls.risk_level;
    entry.confidence_pct = cls.confidence_pct;
    AuthSystem::recordAnalysis(g_current_user_id, entry);
}

// Reads a pasted sequence (raw bases or FASTA, any number of lines).
// Finish with an empty line, a line containing END, or Ctrl+Z / Ctrl+D.
std::string readPastedSequence() {
    std::string all, line;
    while (std::getline(std::cin, line)) {
        std::string t = line;
        while (!t.empty() && (t.back() == '\r' || t.back() == ' ' || t.back() == '\t')) t.pop_back();
        if (t == "END" || t == "end") break;
        if (t.empty()) { if (all.empty()) continue; else break; }
        all += t + "\n";
    }
    if (std::cin.eof()) std::cin.clear();
    return all;
}

std::string trimQuotes(std::string s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == ' ' || s.back() == '\t')) s.pop_back();
    size_t start = s.find_first_not_of(" \t");
    s = (start == std::string::npos) ? "" : s.substr(start);
    if (s.size() >= 2 && ((s.front() == '"' && s.back() == '"') || (s.front() == '\'' && s.back() == '\''))) s = s.substr(1, s.size() - 2);
    return s;
}

void printHelp() {
    std::cout << "GenomiX C++ Clinical Platform\n";
    std::cout << "Usage:\n";
    std::cout << "  genomix.exe --web [port]             Launch embedded Web UI (Default port 8080)\n";
    std::cout << "  genomix.exe --benchmark              Run batch benchmark across all NCBI GenBank samples\n";
    std::cout << "  genomix.exe --file <path> [disease]  Analyze a local FASTA file\n";
    std::cout << "  genomix.exe --dna <seq> [disease]    Analyze a raw nucleotide string\n";
    std::cout << "  genomix.exe                          Interactive Terminal Console\n";
}

int main(int argc, char* argv[]) {
    // Resolve project directories relative to executable or current working directory
    // Works from the project folder, from inside bin\, or wherever the exe lives (folder may be moved/renamed).
    std::string base_dir = ".";
    if (!std::filesystem::exists("./data")) {
        std::error_code ec;
        std::filesystem::path exe_root = std::filesystem::absolute(argv[0], ec).parent_path().parent_path();
        if (!ec && std::filesystem::exists(exe_root / "data")) base_dir = exe_root.string();
        else if (std::filesystem::exists("../data")) base_dir = "..";
    }

    std::string data_dir = base_dir + "/data";
    std::string results_dir = base_dir + "/results";
    std::string web_dir = base_dir + "/web";

    AuthSystem::init(results_dir);

    if (argc > 1) {
        std::string arg1 = argv[1];
        if (arg1 == "--web") {
            int port = 8080;
            if (argc > 2) {
                try { port = std::stoi(argv[2]); }
                catch (...) { std::cerr << "[ERROR] Invalid port '" << argv[2] << "'. Using 8080.\n"; }
                if (port < 1 || port > 65535) { std::cerr << "[ERROR] Port out of range. Using 8080.\n"; port = 8080; }
            }
            WebServer::start(port, web_dir, data_dir);
            return 0;
        } else if (arg1 == "--benchmark") {
            runBenchmark(data_dir);
            return 0;
        } else if (arg1 == "--file" && argc > 2) {
            std::string fpath = argv[2];
            std::string d_key = (argc > 3) ? argv[3] : "";
            analyzeSingle(fpath, true, d_key, data_dir);
            return 0;
        } else if (arg1 == "--dna" && argc > 2) {
            std::string dna = argv[2];
            std::string d_key = (argc > 3) ? argv[3] : "";
            analyzeSingle(dna, false, d_key, data_dir);
            return 0;
        } else if (arg1 == "--help" || arg1 == "-h") {
            printHelp();
            return 0;
        }
    }

    // Default Interactive Console Menu
    std::cout << "================================================================================\n";
    std::cout << "      GENOMIX | DNA MOTIF DISCOVERY & CLINICAL SUITE (C++ EDITION)              \n";
    std::cout << "      Design and Analysis of Algorithms - Pure Linear-Time Z-Algorithm Engine   \n";
    std::cout << "================================================================================\n";

    while (true) {
        std::cout << "\n[MENU OPTIONS]\n";
        std::cout << "  1. Launch Web Dashboard UI (http://localhost:8080)\n";
        std::cout << "  2. Run Batch NCBI Benchmark Across All Cohort Samples\n";
        std::cout << "  3. Screen Local FASTA File\n";
        std::cout << "  4. Analyze Raw DNA Nucleotide Sequence\n";
        std::cout << "  5. View Clinical Disease Thresholds (ACMG Registry)\n";
        std::cout << "  6. Practitioner Sign In / Sign Up\n";
        std::cout << "  7. View User Analysis History & Audit Logs\n";
        std::cout << "  0. Exit\n";
        std::cout << "Select Option [0-7]: ";

        int choice = -1;
        if (!(std::cin >> choice)) {
            std::cin.clear();
            std::string dump;
            std::getline(std::cin, dump);
            continue;
        }

        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if (choice == 0) {
            std::cout << "Exiting GenomiX. Goodbye.\n";
            break;
        } else if (choice == 1) {
            WebServer::start(8080, web_dir, data_dir);
        } else if (choice == 2) {
            runBenchmark(data_dir);
        } else if (choice == 3) {
            std::cout << "Enter path to FASTA file: ";
            std::string fpath;
            std::getline(std::cin, fpath);
            fpath = trimQuotes(fpath);
            std::cout << "Enter Disease Panel (or press Enter for auto-detect): ";
            std::string d_key;
            std::getline(std::cin, d_key);
            analyzeSingle(fpath, true, trimQuotes(d_key), data_dir);
        } else if (choice == 4) {
            std::cout << "\nPaste your DNA sequence (raw bases or FASTA, multiple lines are fine).\n";
            std::cout << "When finished, press Enter on an empty line (or type END):\n";
            std::string seq = readPastedSequence();
            if (seq.empty()) {
                std::cout << "[ERROR] Nothing was pasted.\n";
            } else {
                std::cout << "Enter Disease Panel (e.g. Huntington_Disease, Fragile_X_Syndrome) or press Enter for auto-detect: ";
                std::string d_key;
                std::getline(std::cin, d_key);
                analyzeSingle(seq, false, trimQuotes(d_key), data_dir);
            }
        } else if (choice == 5) {
            std::cout << "\nACMG Clinical Thresholds Registry:\n";
            for (const auto& [k, d] : DiseaseRegistry::getAllDiseases()) {
                std::cout << "  * " << d.display_name << " (" << d.gene << " / " << d.motif << " on " << d.chromosome << "):\n"
                          << "      Normal Limit: <= " << d.thresholds.normal_max << " repeats\n"
                          << "      Intermediate: <= " << d.thresholds.intermediate_max << " repeats\n"
                          << "      Pathogenic  : >= " << d.thresholds.pathogenic_min << " repeats\n";
            }
        } else if (choice == 6) {
            std::cout << "1. Sign In\n2. Create Account\nChoice: ";
            int auth_choice = 1;
            if (!(std::cin >> auth_choice)) { std::cin.clear(); std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); continue; }
            if (auth_choice == 1) {
                std::string u, p;
                std::cout << "Username/Email: ";
                std::cin >> u;
                std::cout << "Password: ";
                std::cin >> p;
                auto user = AuthSystem::authenticate(u, p);
                if (user.has_value()) {
                    g_current_user_id = user->id;
                    std::cout << "[SUCCESS] Welcome, " << user->full_name << " (" << user->institution << ")\n";
                } else {
                    std::cout << "[ERROR] Invalid credentials.\n";
                }
            } else {
                std::string u, email, pass, name, inst;
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                std::cout << "Full Name: "; std::getline(std::cin, name);
                std::cout << "Institution: "; std::getline(std::cin, inst);
                std::cout << "Username: "; std::cin >> u;
                std::cout << "Email: "; std::cin >> email;
                std::cout << "Password: "; std::cin >> pass;
                auto [ok, msg] = AuthSystem::registerUser(u, email, pass, name, inst);
                std::cout << (ok ? "[SUCCESS] " : "[ERROR] ") << msg << "\n";
            }
        } else if (choice == 7) {
            std::cout << "\nRecent Diagnostic Audit Logs:\n";
            auto logs = AuthSystem::getUserAnalyses(g_current_user_id);
            if (logs.empty()) {
                std::cout << "  No audit logs recorded yet.\n";
            } else {
                for (const auto& l : logs) {
                    std::cout << "  [" << l.created_at << "] " << l.sample_name << " | " << l.disease 
                              << " | " << l.max_repeats << " repeats | " << l.prediction_label 
                              << " (" << l.risk_level << ")\n";
                }
            }
        }
    }

    return 0;
}

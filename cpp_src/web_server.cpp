#define WIN32_LEAN_AND_MEAN
#include <winsock2.h>
#include <ws2tcpip.h>
#include "web_server.hpp"
#include "z_algorithm.hpp"
#include "fasta_parser.hpp"
#include "tandem_analyzer.hpp"
#include "statistical_model.hpp"
#include "disease_registry.hpp"
#include "classifier.hpp"
#include "auth_system.hpp"
#include "pathology_report.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include <filesystem>
#include <algorithm>

#pragma comment(lib, "ws2_32.lib")

namespace genomix {

namespace {

std::string extractJsonField(const std::string& json, const std::string& field) {
    std::string key = "\"" + field + "\"";
    size_t pos = json.find(key);
    if (pos == std::string::npos) return "";

    pos = json.find(':', pos);
    if (pos == std::string::npos) return "";

    pos = json.find_first_not_of(" \t\r\n", pos + 1);
    if (pos == std::string::npos) return "";

    if (json[pos] == '\"') {
        size_t end_pos = json.find('\"', pos + 1);
        if (end_pos == std::string::npos) return "";
        return json.substr(pos + 1, end_pos - (pos + 1));
    } else {
        size_t end_pos = json.find_first_of(",}\r\n", pos);
        if (end_pos == std::string::npos) end_pos = json.length();
        return json.substr(pos, end_pos - pos);
    }
}

std::string urlDecode(const std::string& str) {
    std::string ret;
    char ch;
    int i, ii;
    for (i = 0; i < static_cast<int>(str.length()); i++) {
        if (str[i] == '%') {
            sscanf(str.substr(i + 1, 2).c_str(), "%x", &ii);
            ch = static_cast<char>(ii);
            ret += ch;
            i = i + 2;
        } else if (str[i] == '+') {
            ret += ' ';
        } else {
            ret += str[i];
        }
    }
    return ret;
}

std::string buildHttpResponse(const std::string& content, const std::string& content_type = "application/json", int status = 200) {
    std::ostringstream oss;
    oss << "HTTP/1.1 " << status << " OK\r\n";
    oss << "Content-Type: " << content_type << "; charset=utf-8\r\n";
    oss << "Content-Length: " << content.length() << "\r\n";
    oss << "Access-Control-Allow-Origin: *\r\n";
    oss << "Access-Control-Allow-Methods: GET, POST, OPTIONS\r\n";
    oss << "Access-Control-Allow-Headers: Content-Type\r\n";
    oss << "Connection: close\r\n\r\n";
    oss << content;
    return oss.str();
}

void handleClient(SOCKET client_socket, std::string web_dir, std::string data_dir) {
    char buffer[65536];
    int bytes_received = recv(client_socket, buffer, sizeof(buffer) - 1, 0);

    if (bytes_received <= 0) {
        closesocket(client_socket);
        return;
    }

    buffer[bytes_received] = '\0';
    std::string request(buffer);

    std::istringstream req_stream(request);
    std::string method, path, proto;
    req_stream >> method >> path >> proto;

    // Handle OPTIONS for CORS
    if (method == "OPTIONS") {
        std::string resp = buildHttpResponse("", "text/plain", 204);
        send(client_socket, resp.c_str(), resp.length(), 0);
        closesocket(client_socket);
        return;
    }

    // Split headers and body
    size_t body_pos = request.find("\r\n\r\n");
    std::string body = (body_pos != std::string::npos) ? request.substr(body_pos + 4) : "";

    // 1. GET or HEAD / -> Serve Web UI (index.html)
    if ((method == "GET" || method == "HEAD") && (path == "/" || path == "/index.html")) {
        std::string html_path = web_dir + "/index.html";
        std::ifstream f(html_path, std::ios::binary);
        if (f.is_open()) {
            std::ostringstream ss;
            ss << f.rdbuf();
            std::string content = (method == "HEAD") ? "" : ss.str();
            std::string resp = buildHttpResponse(content, "text/html");
            send(client_socket, resp.c_str(), resp.length(), 0);
        } else {
            std::string resp = buildHttpResponse("<h1>404 UI Not Found</h1>", "text/html", 404);
            send(client_socket, resp.c_str(), resp.length(), 0);
        }
        closesocket(client_socket);
        return;
    }

    // 2. POST /api/login
    if (method == "POST" && path == "/api/login") {
        std::string username = extractJsonField(body, "username");
        std::string password = extractJsonField(body, "password");

        auto user_opt = AuthSystem::authenticate(username, password);
        std::ostringstream oss;
        if (user_opt.has_value()) {
            const auto& u = user_opt.value();
            oss << "{\"success\":true,\"user\":{\"id\":" << u.id 
                << ",\"username\":\"" << u.username << "\""
                << ",\"full_name\":\"" << u.full_name << "\""
                << ",\"email\":\"" << u.email << "\""
                << ",\"institution\":\"" << u.institution << "\"}}";
        } else {
            oss << "{\"success\":false,\"message\":\"Invalid credentials.\"}";
        }
        std::string resp = buildHttpResponse(oss.str());
        send(client_socket, resp.c_str(), resp.length(), 0);
        closesocket(client_socket);
        return;
    }

    // 3. POST /api/register
    if (method == "POST" && path == "/api/register") {
        std::string username = extractJsonField(body, "username");
        std::string email = extractJsonField(body, "email");
        std::string password = extractJsonField(body, "password");
        std::string full_name = extractJsonField(body, "full_name");
        std::string institution = extractJsonField(body, "institution");

        auto [ok, msg] = AuthSystem::registerUser(username, email, password, full_name, institution);
        std::ostringstream oss;
        oss << "{\"success\":" << (ok ? "true" : "false") 
            << ",\"message\":\"" << msg << "\"}";

        std::string resp = buildHttpResponse(oss.str());
        send(client_socket, resp.c_str(), resp.length(), 0);
        closesocket(client_socket);
        return;
    }

    // 4. GET /api/files -> List pre-loaded NCBI files in data directory
    if (method == "GET" && path == "/api/files") {
        std::ostringstream oss;
        oss << "{\"files\":[";
        bool first = true;
        if (std::filesystem::exists(data_dir)) {
            for (const auto& entry : std::filesystem::directory_iterator(data_dir)) {
                if (entry.path().extension() == ".fasta" || entry.path().extension() == ".fa") {
                    if (!first) oss << ",";
                    oss << "\"" << entry.path().filename().string() << "\"";
                    first = false;
                }
            }
        }
        oss << "]}";
        std::string resp = buildHttpResponse(oss.str());
        send(client_socket, resp.c_str(), resp.length(), 0);
        closesocket(client_socket);
        return;
    }

    // 5. POST /api/analyze
    if (method == "POST" && path == "/api/analyze") {
        std::string disease_key = extractJsonField(body, "disease");
        std::string input_mode = extractJsonField(body, "mode"); // "file" or "dna"
        std::string payload = extractJsonField(body, "payload");
        std::string user_id_str = extractJsonField(body, "user_id");
        int user_id = user_id_str.empty() ? 1 : std::stoi(user_id_str);

        if (disease_key.empty()) disease_key = "Huntington_Disease";
        const auto& disease = DiseaseRegistry::getDisease(disease_key);

        FastaRecord record;
        std::string specimen_name;

        if (input_mode == "file") {
            std::string file_path = data_dir + "/" + payload;
            record = FastaParser::readFile(file_path);
            specimen_name = payload;
        } else {
            record = FastaParser::parseString(payload, "Pasted Nucleotides");
            specimen_name = "Pasted Nucleotides";
        }

        std::string resp;
        if (!record.success) {
            std::ostringstream err_oss;
            err_oss << "{\"success\":false,\"message\":\"" << record.error_message << "\"}";
            resp = buildHttpResponse(err_oss.str(), "application/json", 400);
        } else {
            // Run Z-Algorithm Engine
            ZSearchResult z_res = ZAlgorithm::search(disease.motif, record.sequence);

            // Run Tandem Repeat Analyzer
            TandemAnalysisResult tandem = TandemAnalyzer::analyze(z_res.match_indices, disease.motif.length());

            // Run Statistical Background Model
            StatisticalProfile stat_prof = StatisticalModel::evaluate(record.sequence, disease.motif, z_res.total_occurrences);

            // Run Clinical Classifier
            ClassificationResult classification = ClinicalClassifier::classify(
                disease,
                tandem.max_repeat_count,
                stat_prof,
                z_res.efficiency_ratio
            );

            // Record to User Analytics History
            AnalysisLog log_entry;
            log_entry.sample_name = specimen_name;
            log_entry.disease = disease.key;
            log_entry.target_motif = disease.motif;
            log_entry.seq_length = record.sequence.length();
            log_entry.max_repeats = tandem.max_repeat_count;
            log_entry.total_matches = z_res.total_occurrences;
            log_entry.z_comparisons = z_res.z_comparisons;
            log_entry.comparisons_saved = z_res.comparisons_saved;
            log_entry.prediction_label = classification.predicted_label;
            log_entry.risk_level = classification.risk_level;
            log_entry.confidence_pct = classification.confidence_pct;
            AuthSystem::recordAnalysis(user_id, log_entry);

            // Build Json Response
            std::ostringstream oss;
            oss << "{\"success\":true,\"data\":";
            oss << PathologyReport::formatJsonReport(
                specimen_name,
                disease,
                record.sequence.length(),
                z_res,
                tandem.max_repeat_count,
                tandem.longest_block_start,
                tandem.longest_block_end,
                stat_prof,
                classification
            );
            oss << "}";
            resp = buildHttpResponse(oss.str());
        }

        send(client_socket, resp.c_str(), resp.length(), 0);
        closesocket(client_socket);
        return;
    }

    // 6. GET /api/history?user_id=X
    if (method == "GET" && path.rfind("/api/history", 0) == 0) {
        int user_id = 1;
        size_t q_pos = path.find("user_id=");
        if (q_pos != std::string::npos) {
            user_id = std::stoi(path.substr(q_pos + 8));
        }

        auto logs = AuthSystem::getUserAnalyses(user_id);
        std::ostringstream oss;
        oss << "{\"history\":[";
        for (size_t i = 0; i < logs.size(); ++i) {
            if (i > 0) oss << ",";
            const auto& l = logs[i];
            oss << "{\"id\":" << l.id
                << ",\"sample_name\":\"" << l.sample_name << "\""
                << ",\"disease\":\"" << l.disease << "\""
                << ",\"motif\":\"" << l.target_motif << "\""
                << ",\"seq_length\":" << l.seq_length
                << ",\"max_repeats\":" << l.max_repeats
                << ",\"total_matches\":" << l.total_matches
                << ",\"z_comparisons\":" << l.z_comparisons
                << ",\"comparisons_saved\":" << l.comparisons_saved
                << ",\"prediction\":\"" << l.prediction_label << "\""
                << ",\"risk_level\":\"" << l.risk_level << "\""
                << ",\"confidence_pct\":" << l.confidence_pct
                << ",\"created_at\":\"" << l.created_at << "\"}";
        }
        oss << "]}";
        std::string resp = buildHttpResponse(oss.str());
        send(client_socket, resp.c_str(), resp.length(), 0);
        closesocket(client_socket);
        return;
    }

    // 7. GET /api/cohort -> Benchmark summary of all data/ files
    if (method == "GET" && path == "/api/cohort") {
        std::ostringstream oss;
        oss << "{\"cohort\":[";
        bool first = true;

        if (std::filesystem::exists(data_dir)) {
            for (const auto& entry : std::filesystem::directory_iterator(data_dir)) {
                if (entry.path().extension() == ".fasta" || entry.path().extension() == ".fa") {
                    std::string fname = entry.path().filename().string();
                    std::string d_key = DiseaseRegistry::detectDiseaseKey(fname);
                    const auto& disease = DiseaseRegistry::getDisease(d_key);

                    FastaRecord rec = FastaParser::readFile(entry.path().string());
                    if (rec.success) {
                        ZSearchResult z_res = ZAlgorithm::search(disease.motif, rec.sequence);
                        TandemAnalysisResult tandem = TandemAnalyzer::analyze(z_res.match_indices, disease.motif.length());
                        StatisticalProfile stat_prof = StatisticalModel::evaluate(rec.sequence, disease.motif, z_res.total_occurrences);
                        ClassificationResult cls = ClinicalClassifier::classify(disease, tandem.max_repeat_count, stat_prof, z_res.efficiency_ratio);

                        if (!first) oss << ",";
                        oss << "{\"filename\":\"" << fname << "\""
                            << ",\"accession\":\"" << fname.substr(0, fname.find('_')) << "\""
                            << ",\"disease\":\"" << disease.display_name << "\""
                            << ",\"gene\":\"" << disease.gene << "\""
                            << ",\"motif\":\"" << disease.motif << "\""
                            << ",\"length\":" << rec.sequence.length()
                            << ",\"max_repeats\":" << tandem.max_repeat_count
                            << ",\"status\":\"" << cls.clinical_status << "\"}";
                        first = false;
                    }
                }
            }
        }
        oss << "]}";
        std::string resp = buildHttpResponse(oss.str());
        send(client_socket, resp.c_str(), resp.length(), 0);
        closesocket(client_socket);
        return;
    }

    // Fallback 404
    std::string resp = buildHttpResponse("{\"error\":\"Not Found\"}", "application/json", 404);
    send(client_socket, resp.c_str(), resp.length(), 0);
    closesocket(client_socket);
}

} // anonymous namespace

void WebServer::start(int port, const std::string& web_dir, const std::string& data_dir) {
    WSADATA wsaData;
    int wsaerr = WSAStartup(MAKEWORD(2, 2), &wsaData);
    if (wsaerr != 0) {
        std::cerr << "[ERROR] WSAStartup failed with error: " << wsaerr << std::endl;
        return;
    }

    SOCKET server_socket = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if (server_socket == INVALID_SOCKET) {
        std::cerr << "[ERROR] Socket creation failed: " << WSAGetLastError() << std::endl;
        WSACleanup();
        return;
    }

    sockaddr_in server_addr;
    server_addr.sin_family = AF_INET;
    server_addr.sin_addr.s_addr = INADDR_ANY;
    server_addr.sin_port = htons(port);

    if (bind(server_socket, (sockaddr*)&server_addr, sizeof(server_addr)) == SOCKET_ERROR) {
        std::cerr << "[ERROR] Bind failed with error: " << WSAGetLastError() << std::endl;
        closesocket(server_socket);
        WSACleanup();
        return;
    }

    if (listen(server_socket, SOMAXCONN) == SOCKET_ERROR) {
        std::cerr << "[ERROR] Listen failed with error: " << WSAGetLastError() << std::endl;
        closesocket(server_socket);
        WSACleanup();
        return;
    }

    std::cout << "================================================================================" << std::endl;
    std::cout << "     GENOMIX C++ CLINICAL SUITE - EMBEDDED HIGH-PERFORMANCE HTTP SERVER         " << std::endl;
    std::cout << "================================================================================" << std::endl;
    std::cout << "  Web UI Address   : http://localhost:" << port << std::endl;
    std::cout << "  Static Web Root  : " << web_dir << std::endl;
    std::cout << "  Data Repository  : " << data_dir << std::endl;
    std::cout << "  Engine           : Pure C++17 Linear-Time Z-Algorithm & Analytical Model      " << std::endl;
    std::cout << "================================================================================" << std::endl;
    std::cout << "[INFO] Server listening on port " << port << "... Press Ctrl+C to terminate." << std::endl;

    while (true) {
        SOCKET client_socket = accept(server_socket, NULL, NULL);
        if (client_socket == INVALID_SOCKET) {
            continue;
        }

        std::thread([client_socket, web_dir, data_dir]() {
            handleClient(client_socket, web_dir, data_dir);
        }).detach();
    }

    closesocket(server_socket);
    WSACleanup();
}

} // namespace genomix

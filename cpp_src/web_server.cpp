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
#include <cstring>
#include <cctype>
#include <cstdio>

#pragma comment(lib, "ws2_32.lib")

namespace genomix {

namespace {

// Decodes a JSON string starting right after the opening quote at `pos`.
// Handles \" \\ \/ \b \f \n \r \t and \uXXXX (encoded as UTF-8).
std::string decodeJsonString(const std::string& json, size_t pos) {
    std::string out;
    for (size_t i = pos; i < json.size(); ++i) {
        char c = json[i];
        if (c == '"') return out;
        if (c != '\\') { out.push_back(c); continue; }
        if (++i >= json.size()) break;
        switch (json[i]) {
            case 'n': out.push_back('\n'); break;
            case 'r': out.push_back('\r'); break;
            case 't': out.push_back('\t'); break;
            case 'b': out.push_back('\b'); break;
            case 'f': out.push_back('\f'); break;
            case 'u': {
                if (i + 4 >= json.size()) return out;
                unsigned cp = 0;
                try { cp = static_cast<unsigned>(std::stoul(json.substr(i + 1, 4), nullptr, 16)); } catch (...) { cp = '?'; }
                i += 4;
                if (cp < 0x80) out.push_back(static_cast<char>(cp));
                else if (cp < 0x800) { out.push_back(static_cast<char>(0xC0 | (cp >> 6))); out.push_back(static_cast<char>(0x80 | (cp & 0x3F))); }
                else { out.push_back(static_cast<char>(0xE0 | (cp >> 12))); out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F))); out.push_back(static_cast<char>(0x80 | (cp & 0x3F))); }
                break;
            }
            default: out.push_back(json[i]); break;   // \" \\ \/
        }
    }
    return out;
}

std::string extractJsonField(const std::string& json, const std::string& field) {
    std::string key = "\"" + field + "\"";
    size_t pos = json.find(key);
    if (pos == std::string::npos) return "";

    pos = json.find(':', pos + key.size());
    if (pos == std::string::npos) return "";

    pos = json.find_first_not_of(" \t\r\n", pos + 1);
    if (pos == std::string::npos) return "";

    if (json[pos] == '"') {
        return decodeJsonString(json, pos + 1);
    }
    size_t end_pos = json.find_first_of(",}\r\n", pos);
    if (end_pos == std::string::npos) end_pos = json.length();
    std::string v = json.substr(pos, end_pos - pos);
    while (!v.empty() && (v.back() == ' ' || v.back() == '\t')) v.pop_back();
    return v;
}

std::string jsonEscape(const std::string& in) {
    std::string out;
    for (unsigned char c : in) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) { char buf[8]; snprintf(buf, sizeof(buf), "\\u%04x", c); out += buf; }
                else out.push_back(static_cast<char>(c));
        }
    }
    return out;
}

int parseIntOr(const std::string& s, int fallback) {
    try { size_t used = 0; int v = std::stoi(s, &used); return used ? v : fallback; }
    catch (...) { return fallback; }
}

// File names from the browser may only name a FASTA file inside data_dir.
bool isSafeFastaName(const std::string& name) {
    if (name.empty() || name.size() > 200) return false;
    if (name.find_first_of("/\\:") != std::string::npos) return false;
    if (name.find("..") != std::string::npos) return false;
    std::filesystem::path ext = std::filesystem::path(name).extension();
    return ext == ".fasta" || ext == ".fa";
}

const char* statusText(int code) {
    switch (code) {
        case 200: return "OK";
        case 204: return "No Content";
        case 400: return "Bad Request";
        case 404: return "Not Found";
        case 413: return "Payload Too Large";
        case 500: return "Internal Server Error";
        default:  return "OK";
    }
}

std::string buildHttpResponse(const std::string& content, const std::string& content_type = "application/json", int status = 200) {
    std::ostringstream oss;
    oss << "HTTP/1.1 " << status << " " << statusText(status) << "\r\n";
    oss << "Content-Type: " << content_type << "; charset=utf-8\r\n";
    oss << "Content-Length: " << content.length() << "\r\n";
    oss << "Access-Control-Allow-Origin: *\r\n";
    oss << "Access-Control-Allow-Methods: GET, POST, OPTIONS\r\n";
    oss << "Access-Control-Allow-Headers: Content-Type\r\n";
    oss << "Connection: close\r\n\r\n";
    oss << content;
    return oss.str();
}

bool readFullRequest(SOCKET sock, std::string& request, bool& too_large) {
    const size_t kMaxRequest = 16u * 1024u * 1024u;   // 16 MB cap
    too_large = false;
    char buf[16384];
    size_t header_end = std::string::npos;
    size_t content_length = 0;
    while (true) {
        int n = recv(sock, buf, sizeof(buf), 0);
        if (n <= 0) break;
        request.append(buf, static_cast<size_t>(n));
        if (request.size() > kMaxRequest) { too_large = true; return false; }
        if (header_end == std::string::npos) {
            header_end = request.find("\r\n\r\n");
            if (header_end != std::string::npos) {
                std::string head = request.substr(0, header_end);
                std::transform(head.begin(), head.end(), head.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                size_t cl = head.find("content-length:");
                if (cl != std::string::npos) {
                    try { content_length = static_cast<size_t>(std::stoull(head.substr(cl + 15))); } catch (...) { content_length = 0; }
                }
                if (content_length > kMaxRequest) { too_large = true; return false; }
            }
        }
        if (header_end != std::string::npos && request.size() >= header_end + 4 + content_length) break;
    }
    return !request.empty();
}

void handleClientImpl(SOCKET client_socket, const std::string& web_dir, const std::string& data_dir);

void handleClient(SOCKET client_socket, std::string web_dir, std::string data_dir) {
    // An uncaught exception inside a detached thread would terminate the whole server.
    try {
        handleClientImpl(client_socket, web_dir, data_dir);
    } catch (const std::exception& e) {
        std::cerr << "[WARN] Request failed: " << e.what() << std::endl;
        std::string resp = buildHttpResponse("{\"success\":false,\"message\":\"Internal server error.\"}", "application/json", 500);
        send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
        closesocket(client_socket);
    } catch (...) {
        closesocket(client_socket);
    }
}

void handleClientImpl(SOCKET client_socket, const std::string& web_dir, const std::string& data_dir) {
    std::string request;
    bool too_large = false;
    if (!readFullRequest(client_socket, request, too_large)) {
        if (too_large) {
            std::string resp = buildHttpResponse("{\"success\":false,\"message\":\"Request too large.\"}", "application/json", 413);
            send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
        }
        closesocket(client_socket);
        return;
    }

    std::istringstream req_stream(request);
    std::string method, path, proto;
    req_stream >> method >> path >> proto;

    // Handle OPTIONS for CORS
    if (method == "OPTIONS") {
        std::string resp = buildHttpResponse("", "text/plain", 204);
        send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
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
            send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
        } else {
            std::string resp = buildHttpResponse("<h1>404 UI Not Found</h1>", "text/html", 404);
            send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
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
                << ",\"username\":\"" << jsonEscape(u.username) << "\""
                << ",\"full_name\":\"" << jsonEscape(u.full_name) << "\""
                << ",\"email\":\"" << jsonEscape(u.email) << "\""
                << ",\"institution\":\"" << jsonEscape(u.institution) << "\"}}";
        } else {
            oss << "{\"success\":false,\"message\":\"Invalid credentials.\"}";
        }
        std::string resp = buildHttpResponse(oss.str());
        send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
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
            << ",\"message\":\"" << jsonEscape(msg) << "\"}";

        std::string resp = buildHttpResponse(oss.str());
        send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
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
        send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
        closesocket(client_socket);
        return;
    }

    // 5. POST /api/analyze
    if (method == "POST" && path == "/api/analyze") {
        std::string disease_key = extractJsonField(body, "disease");
        std::string input_mode = extractJsonField(body, "mode"); // "file" or "dna"
        std::string payload = extractJsonField(body, "payload");
        std::string user_id_str = extractJsonField(body, "user_id");
        int user_id = parseIntOr(user_id_str, 1);

        FastaRecord record;
        std::string specimen_name;

        if (input_mode == "file") {
            if (!isSafeFastaName(payload)) {
                std::string resp = buildHttpResponse("{\"success\":false,\"message\":\"Invalid specimen file name.\"}", "application/json", 400);
                send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
                closesocket(client_socket);
                return;
            }
            record = FastaParser::readFile(data_dir + "/" + payload);
            specimen_name = payload;
        } else {
            record = FastaParser::parseString(payload, "Pasted Nucleotides");
            specimen_name = "Pasted Nucleotides";
        }

        // "auto" (or empty): work out the panel from the specimen instead of silently assuming Huntington's,
        // which produced false "Normal" results for Fragile X / Friedreich / DM1 files.
        if (disease_key.empty() || disease_key == "auto") {
            disease_key = record.success ? DiseaseRegistry::detectDiseaseKeyForSequence(specimen_name, record.header, record.sequence)
                                         : "Huntington_Disease";
        }
        const auto& disease = DiseaseRegistry::getDisease(disease_key);

        std::string resp;
        if (!record.success) {
            std::ostringstream err_oss;
            err_oss << "{\"success\":false,\"message\":\"" << jsonEscape(record.error_message.empty() ? std::string("Could not read a DNA sequence from the input.") : record.error_message) << "\"}";
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
            oss << ",\"text_report\":\"" << jsonEscape(PathologyReport::formatTerminalReport(
                specimen_name,
                disease,
                record.sequence.length(),
                z_res,
                tandem.max_repeat_count,
                tandem.longest_block_start,
                tandem.longest_block_end,
                stat_prof,
                classification
            )) << "\"";
            oss << "}";
            resp = buildHttpResponse(oss.str());
        }

        send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
        closesocket(client_socket);
        return;
    }

    // 6. GET /api/history?user_id=X
    if (method == "GET" && path.rfind("/api/history", 0) == 0) {
        int user_id = 1;
        size_t q_pos = path.find("user_id=");
        if (q_pos != std::string::npos) {
            user_id = parseIntOr(path.substr(q_pos + 8), 1);
        }

        auto logs = AuthSystem::getUserAnalyses(user_id);
        std::ostringstream oss;
        oss << "{\"history\":[";
        for (size_t i = 0; i < logs.size(); ++i) {
            if (i > 0) oss << ",";
            const auto& l = logs[i];
            oss << "{\"id\":" << l.id
                << ",\"sample_name\":\"" << jsonEscape(l.sample_name) << "\""
                << ",\"disease\":\"" << jsonEscape(l.disease) << "\""
                << ",\"motif\":\"" << jsonEscape(l.target_motif) << "\""
                << ",\"seq_length\":" << l.seq_length
                << ",\"max_repeats\":" << l.max_repeats
                << ",\"total_matches\":" << l.total_matches
                << ",\"z_comparisons\":" << l.z_comparisons
                << ",\"comparisons_saved\":" << l.comparisons_saved
                << ",\"prediction\":\"" << jsonEscape(l.prediction_label) << "\""
                << ",\"risk_level\":\"" << jsonEscape(l.risk_level) << "\""
                << ",\"confidence_pct\":" << l.confidence_pct
                << ",\"created_at\":\"" << l.created_at << "\"}";
        }
        oss << "]}";
        std::string resp = buildHttpResponse(oss.str());
        send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
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
                    FastaRecord rec = FastaParser::readFile(entry.path().string());
                    std::string d_key = rec.success
                        ? DiseaseRegistry::detectDiseaseKeyForSequence(fname, rec.header, rec.sequence)
                        : DiseaseRegistry::detectDiseaseKey(fname);
                    const auto& disease = DiseaseRegistry::getDisease(d_key);

                    if (rec.success) {
                        ZSearchResult z_res = ZAlgorithm::search(disease.motif, rec.sequence);
                        TandemAnalysisResult tandem = TandemAnalyzer::analyze(z_res.match_indices, disease.motif.length());
                        StatisticalProfile stat_prof = StatisticalModel::evaluate(rec.sequence, disease.motif, z_res.total_occurrences);
                        ClassificationResult cls = ClinicalClassifier::classify(disease, tandem.max_repeat_count, stat_prof, z_res.efficiency_ratio);

                        if (!first) oss << ",";
                        oss << "{\"filename\":\"" << jsonEscape(fname) << "\""
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
        send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
        closesocket(client_socket);
        return;
    }

    // Fallback 404
    std::string resp = buildHttpResponse("{\"error\":\"Not Found\"}", "application/json", 404);
    send(client_socket, resp.c_str(), static_cast<int>(resp.length()), 0);
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

    // Prevent another process from hijacking the port (SO_REUSEADDR is unsafe on Windows).
    int opt = 1;
    setsockopt(server_socket, SOL_SOCKET, SO_EXCLUSIVEADDRUSE, (const char*)&opt, sizeof(opt));

    // Bind to loopback only: no firewall prompt, and the clinical dashboard is not exposed to the LAN.
    // If the port is blocked (WSAEACCES 10013: reserved by Hyper-V/WinNAT, or 10048: already in use),
    // fall back to the next free ports.
    const int kMaxAttempts = 20;
    int current_port = port;
    int last_err = 0;
    bool bound = false;
    for (int attempt = 0; attempt < kMaxAttempts && current_port <= 65535; ++attempt, ++current_port) {
        sockaddr_in server_addr;
        std::memset(&server_addr, 0, sizeof(server_addr));
        server_addr.sin_family = AF_INET;
        server_addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
        server_addr.sin_port = htons(static_cast<u_short>(current_port));
        if (bind(server_socket, (sockaddr*)&server_addr, sizeof(server_addr)) == 0) {
            bound = true;
            break;
        }
        last_err = WSAGetLastError();
    }

    if (!bound) {
        std::cerr << "[ERROR] Bind failed on ports " << port << "-" << (current_port - 1)
                  << " (Winsock error " << last_err << ")." << std::endl;
        std::cerr << "[HINT] Try another port: genomix.exe --web 9000   (check reserved ranges with: netsh int ipv4 show excludedportrange protocol=tcp)" << std::endl;
        closesocket(server_socket);
        WSACleanup();
        return;
    }

    if (current_port != port) {
        std::cout << "[INFO] Port " << port << " unavailable; using port " << current_port << " instead." << std::endl;
    }
    port = current_port;

    if (listen(server_socket, SOMAXCONN) == SOCKET_ERROR) {
        int err = WSAGetLastError();
        std::cerr << "[ERROR] Listen failed with error: " << err << std::endl;
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
            handleClient(client_socket, web_dir, data_dir);   // has its own try/catch
        }).detach();
    }

    closesocket(server_socket);
    WSACleanup();
}

} // namespace genomix

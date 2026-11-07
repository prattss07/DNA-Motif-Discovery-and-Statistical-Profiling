#include "auth_system.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <cstring>
#include <filesystem>
#include <mutex>

namespace {
std::mutex g_auth_mutex;   // request threads share users/logs

// Fields are stored tab-separated, one record per line: strip characters that would corrupt the file.
std::string cleanField(std::string v) {
    for (char& c : v) if (c == '\t' || c == '\n' || c == '\r') c = ' ';
    return v;
}
}

namespace genomix {

// -------------------------------------------------------------
// Pure C++ Standard SHA-256 Implementation (Zero Dependency)
// -------------------------------------------------------------
namespace {
    inline uint32_t rotr(uint32_t x, uint32_t n) {
        return (x >> n) | (x << (32 - n));
    }
    inline uint32_t ch(uint32_t x, uint32_t y, uint32_t z) {
        return (x & y) ^ (~x & z);
    }
    inline uint32_t maj(uint32_t x, uint32_t y, uint32_t z) {
        return (x & y) ^ (x & z) ^ (y & z);
    }
    inline uint32_t sig0(uint32_t x) {
        return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22);
    }
    inline uint32_t sig1(uint32_t x) {
        return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25);
    }
    inline uint32_t theta0(uint32_t x) {
        return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3);
    }
    inline uint32_t theta1(uint32_t x) {
        return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10);
    }

    const uint32_t K[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    };
}

std::string AuthSystem::sha256(const std::string& input) {
    std::string salted = input + "genomiX_cpp_salt_2026";
    std::vector<uint8_t> msg(salted.begin(), salted.end());

    uint64_t bit_len = msg.size() * 8;
    msg.push_back(0x80);
    while ((msg.size() * 8) % 512 != 448) {
        msg.push_back(0x00);
    }
    for (int i = 7; i >= 0; --i) {
        msg.push_back(static_cast<uint8_t>((bit_len >> (i * 8)) & 0xff));
    }

    uint32_t H[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    };

    for (size_t chunk = 0; chunk < msg.size(); chunk += 64) {
        uint32_t W[64];
        for (int i = 0; i < 16; ++i) {
            W[i] = (msg[chunk + i * 4] << 24) |
                   (msg[chunk + i * 4 + 1] << 16) |
                   (msg[chunk + i * 4 + 2] << 8) |
                   (msg[chunk + i * 4 + 3]);
        }
        for (int i = 16; i < 64; ++i) {
            W[i] = theta1(W[i - 2]) + W[i - 7] + theta0(W[i - 15]) + W[i - 16];
        }

        uint32_t a = H[0], b = H[1], c = H[2], d = H[3];
        uint32_t e = H[4], f = H[5], g = H[6], h = H[7];

        for (int i = 0; i < 64; ++i) {
            uint32_t T1 = h + sig1(e) + ch(e, f, g) + K[i] + W[i];
            uint32_t T2 = sig0(a) + maj(a, b, c);
            h = g;
            g = f;
            f = e;
            e = d + T1;
            d = c;
            c = b;
            b = a;
            a = T1 + T2;
        }

        H[0] += a; H[1] += b; H[2] += c; H[3] += d;
        H[4] += e; H[5] += f; H[6] += g; H[7] += h;
    }

    std::ostringstream oss;
    for (int i = 0; i < 8; ++i) {
        oss << std::hex << std::setw(8) << std::setfill('0') << H[i];
    }
    return oss.str();
}

std::string AuthSystem::users_file_path = "";
std::string AuthSystem::logs_file_path = "";
std::vector<UserProfile> AuthSystem::users;
std::vector<AnalysisLog> AuthSystem::logs;

std::string AuthSystem::getCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;
#if defined(_WIN32)
    localtime_s(&tm_buf, &tt);
#else
    localtime_r(&tt, &tm_buf);
#endif
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tm_buf);
    return std::string(buf);
}

bool AuthSystem::init(const std::string& data_dir) {
    std::filesystem::create_directories(data_dir);
    users_file_path = data_dir + "/users.csv";
    logs_file_path = data_dir + "/analysis_logs.csv";

    loadData();

    // Ensure default demo user exists
    if (users.empty()) {
        registerUser("prattss07", "prathampsg266@gmail.com", "SecurePass123!", "Prathamesh P", "Neurogenetics Laboratory");
    }

    return true;
}

void AuthSystem::loadData() {
    users.clear();
    logs.clear();

    // Read Users
    std::ifstream uf(users_file_path);
    if (uf.is_open()) {
        std::string line;
        while (std::getline(uf, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string id_s, user, email, hash, name, inst, created, last_l;
            if (std::getline(ss, id_s, '\t') &&
                std::getline(ss, user, '\t') &&
                std::getline(ss, email, '\t') &&
                std::getline(ss, hash, '\t') &&
                std::getline(ss, name, '\t') &&
                std::getline(ss, inst, '\t') &&
                std::getline(ss, created, '\t') &&
                std::getline(ss, last_l)) {
                try {
                    users.push_back({
                        std::stoi(id_s), user, email, hash, name, inst, created, last_l
                    });
                } catch (...) { /* skip corrupt line */ }
            }
        }
        uf.close();
    }

    // Read Logs
    std::ifstream lf(logs_file_path);
    if (lf.is_open()) {
        std::string line;
        while (std::getline(lf, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string id_s, uid_s, sample, disease, motif, len_s, rep_s, hits_s, comps_s, saved_s, pred, risk, conf_s, created;
            if (std::getline(ss, id_s, '\t') &&
                std::getline(ss, uid_s, '\t') &&
                std::getline(ss, sample, '\t') &&
                std::getline(ss, disease, '\t') &&
                std::getline(ss, motif, '\t') &&
                std::getline(ss, len_s, '\t') &&
                std::getline(ss, rep_s, '\t') &&
                std::getline(ss, hits_s, '\t') &&
                std::getline(ss, comps_s, '\t') &&
                std::getline(ss, saved_s, '\t') &&
                std::getline(ss, pred, '\t') &&
                std::getline(ss, risk, '\t') &&
                std::getline(ss, conf_s, '\t') &&
                std::getline(ss, created)) {
                try {
                    logs.push_back({
                        std::stoi(id_s),
                        std::stoi(uid_s),
                        sample, disease, motif,
                        static_cast<size_t>(std::stoull(len_s)),
                        static_cast<size_t>(std::stoull(rep_s)),
                        static_cast<size_t>(std::stoull(hits_s)),
                        static_cast<uint64_t>(std::stoull(comps_s)),
                        static_cast<uint64_t>(std::stoull(saved_s)),
                        pred, risk,
                        std::stod(conf_s),
                        created
                    });
                } catch (...) { /* skip corrupt line */ }
            }
        }
        lf.close();
    }
}

void AuthSystem::saveData() {
    std::ofstream uf(users_file_path, std::ios::trunc);
    if (uf.is_open()) {
        for (const auto& u : users) {
            uf << u.id << '\t' << u.username << '\t' << u.email << '\t' << u.password_hash << '\t'
               << u.full_name << '\t' << u.institution << '\t' << u.created_at << '\t' << u.last_login << '\n';
        }
        uf.close();
    }

    std::ofstream lf(logs_file_path, std::ios::trunc);
    if (lf.is_open()) {
        for (const auto& l : logs) {
            lf << l.id << '\t' << l.user_id << '\t' << l.sample_name << '\t' << l.disease << '\t'
               << l.target_motif << '\t' << l.seq_length << '\t' << l.max_repeats << '\t'
               << l.total_matches << '\t' << l.z_comparisons << '\t' << l.comparisons_saved << '\t'
               << l.prediction_label << '\t' << l.risk_level << '\t' << l.confidence_pct << '\t'
               << l.created_at << '\n';
        }
        lf.close();
    }
}

std::pair<bool, std::string> AuthSystem::registerUser(
    const std::string& username,
    const std::string& email,
    const std::string& password,
    const std::string& full_name,
    const std::string& institution
) {
    std::lock_guard<std::mutex> lock(g_auth_mutex);
    if (username.empty() || email.empty() || password.empty() || full_name.empty()) {
        return {false, "All required fields must be populated."};
    }

    std::string lower_user = username;
    std::string lower_email = email;
    std::transform(lower_user.begin(), lower_user.end(), lower_user.begin(), ::tolower);
    std::transform(lower_email.begin(), lower_email.end(), lower_email.begin(), ::tolower);

    for (const auto& u : users) {
        std::string existing_user = u.username;
        std::string existing_email = u.email;
        std::transform(existing_user.begin(), existing_user.end(), existing_user.begin(), ::tolower);
        std::transform(existing_email.begin(), existing_email.end(), existing_email.begin(), ::tolower);

        if (existing_user == lower_user) {
            return {false, "Username is already registered."};
        }
        if (existing_email == lower_email) {
            return {false, "Email address is already registered."};
        }
    }

    int next_id = users.empty() ? 1 : (users.back().id + 1);
    std::string hash = sha256(password);
    std::string now_ts = getCurrentTimestamp();

    users.push_back({
        next_id,
        cleanField(username),
        cleanField(email),
        hash,
        cleanField(full_name),
        cleanField(institution),
        now_ts,
        now_ts
    });

    saveData();
    return {true, "Account registered successfully."};
}

std::optional<UserProfile> AuthSystem::authenticate(
    const std::string& username_or_email,
    const std::string& password
) {
    std::lock_guard<std::mutex> lock(g_auth_mutex);
    std::string lower_id = username_or_email;
    std::transform(lower_id.begin(), lower_id.end(), lower_id.begin(), ::tolower);
    std::string hash = sha256(password);

    for (auto& u : users) {
        std::string u_name = u.username;
        std::string u_mail = u.email;
        std::transform(u_name.begin(), u_name.end(), u_name.begin(), ::tolower);
        std::transform(u_mail.begin(), u_mail.end(), u_mail.begin(), ::tolower);

        if ((u_name == lower_id || u_mail == lower_id) && u.password_hash == hash) {
            u.last_login = getCurrentTimestamp();
            saveData();
            return u;
        }
    }
    return std::nullopt;
}

void AuthSystem::recordAnalysis(int user_id, const AnalysisLog& log) {
    std::lock_guard<std::mutex> lock(g_auth_mutex);
    AnalysisLog entry = log;
    entry.sample_name = cleanField(entry.sample_name);
    entry.id = logs.empty() ? 1 : (logs.back().id + 1);
    entry.user_id = user_id;
    entry.created_at = getCurrentTimestamp();
    logs.push_back(entry);
    saveData();
}

std::vector<AnalysisLog> AuthSystem::getUserAnalyses(int user_id) {
    std::lock_guard<std::mutex> lock(g_auth_mutex);
    std::vector<AnalysisLog> res;
    for (const auto& l : logs) {
        if (l.user_id == user_id) {
            res.push_back(l);
        }
    }
    std::reverse(res.begin(), res.end());
    return res;
}

std::vector<UserProfile> AuthSystem::getAllUsers() {
    std::lock_guard<std::mutex> lock(g_auth_mutex);
    return users;
}

} // namespace genomix

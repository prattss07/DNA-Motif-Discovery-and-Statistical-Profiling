#ifndef AUTH_SYSTEM_HPP
#define AUTH_SYSTEM_HPP

#include <string>
#include <vector>
#include <optional>
#include <cstdint>

namespace genomix {

struct UserProfile {
    int id;
    std::string username;
    std::string email;
    std::string password_hash;
    std::string full_name;
    std::string institution;
    std::string created_at;
    std::string last_login;
};

struct AnalysisLog {
    int id;
    int user_id;
    std::string sample_name;
    std::string disease;
    std::string target_motif;
    size_t seq_length;
    size_t max_repeats;
    size_t total_matches;
    uint64_t z_comparisons;
    uint64_t comparisons_saved;
    std::string prediction_label;
    std::string risk_level;
    double confidence_pct;
    std::string created_at;
};

class AuthSystem {
public:
    static std::string sha256(const std::string& input);

    static bool init(const std::string& data_dir);

    static std::pair<bool, std::string> registerUser(
        const std::string& username,
        const std::string& email,
        const std::string& password,
        const std::string& full_name,
        const std::string& institution
    );

    static std::optional<UserProfile> authenticate(
        const std::string& username_or_email,
        const std::string& password
    );

    static void recordAnalysis(int user_id, const AnalysisLog& log);
    static std::vector<AnalysisLog> getUserAnalyses(int user_id);
    static std::vector<UserProfile> getAllUsers();

private:
    static std::string users_file_path;
    static std::string logs_file_path;
    static std::vector<UserProfile> users;
    static std::vector<AnalysisLog> logs;

    static void loadData();
    static void saveData();
    static std::string getCurrentTimestamp();
};

} // namespace genomix

#endif // AUTH_SYSTEM_HPP

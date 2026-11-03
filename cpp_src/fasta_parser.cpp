#include "fasta_parser.hpp"
#include <fstream>
#include <cctype>
#include <filesystem>

namespace genomix {

FastaRecord FastaParser::readFile(const std::string& filepath) {
    FastaRecord record;
    record.success = false;
    record.filename = std::filesystem::path(filepath).filename().string();

    std::ifstream file(filepath);
    if (!file.is_open()) {
        record.error_message = "Failed to open file: " + filepath;
        return record;
    }

    std::string line;
    bool found_header = false;

    while (std::getline(file, line)) {
        // Strip carriage returns and trailing whitespace
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' ' || line.back() == '\t')) {
            line.pop_back();
        }
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!found_header) {
                record.header = line.substr(1);
                found_header = true;
            }
        } else {
            for (char ch : line) {
                if (std::isalpha(static_cast<unsigned char>(ch))) {
                    record.sequence.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
                }
            }
        }
    }

    file.close();
    record.success = !record.sequence.empty();
    if (!record.success) {
        record.error_message = "File contained no valid sequence data.";
    }
    return record;
}

FastaRecord FastaParser::parseString(const std::string& raw_content, const std::string& default_header) {
    FastaRecord record;
    record.filename = "User_Stream";
    record.header = default_header;
    record.success = false;

    size_t pos = 0;
    size_t len = raw_content.length();
    bool in_header = false;

    while (pos < len) {
        size_t next_line = raw_content.find('\n', pos);
        if (next_line == std::string::npos) next_line = len;

        std::string line = raw_content.substr(pos, next_line - pos);
        pos = next_line + 1;

        while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t')) {
            line.pop_back();
        }
        if (line.empty()) continue;

        if (line[0] == '>') {
            record.header = line.substr(1);
        } else {
            for (char ch : line) {
                if (std::isalpha(static_cast<unsigned char>(ch))) {
                    record.sequence.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
                }
            }
        }
    }

    record.success = !record.sequence.empty();
    return record;
}

} // namespace genomix

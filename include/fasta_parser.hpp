#ifndef FASTA_PARSER_HPP
#define FASTA_PARSER_HPP

#include <string>

namespace genomix {

struct FastaRecord {
    std::string header;
    std::string sequence;
    std::string filename;
    bool success;
    std::string error_message;
};

class FastaParser {
public:
    static FastaRecord readFile(const std::string& filepath);
    static FastaRecord parseString(const std::string& raw_content, const std::string& default_header = "Raw Input");
};

} // namespace genomix

#endif // FASTA_PARSER_HPP

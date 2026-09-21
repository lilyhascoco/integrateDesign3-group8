#include "parser/parser_common.h"

#include <algorithm>
#include <cctype>
#include <sstream>

namespace plc::parser_util {

std::string trim(const std::string& s) {
    std::size_t begin = 0;
    while (begin < s.size() && std::isspace(static_cast<unsigned char>(s[begin])) != 0) {
        ++begin;
    }
    std::size_t end = s.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(s[end - 1])) != 0) {
        --end;
    }
    return s.substr(begin, end - begin);
}

std::vector<std::string> split_ws(const std::string& s) {
    std::istringstream iss(s);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

bool starts_with(const std::string& s, const std::string& prefix) {
    return s.rfind(prefix, 0) == 0;
}

bool is_comment_or_empty(const std::string& s) {
    const std::string t = trim(s);
    return t.empty() || starts_with(t, "#") || starts_with(t, "UCLA");
}

plc::Orient parse_orient(const std::string& token) {
    if (token == "N") return plc::Orient::kN;
    if (token == "S") return plc::Orient::kS;
    if (token == "E") return plc::Orient::kE;
    if (token == "W") return plc::Orient::kW;
    if (token == "FN") return plc::Orient::kFN;
    if (token == "FS") return plc::Orient::kFS;
    if (token == "FE") return plc::Orient::kFE;
    if (token == "FW") return plc::Orient::kFW;
    return plc::Orient::kUnknown;
}

}  // namespace plc::parser_util

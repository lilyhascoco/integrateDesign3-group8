#pragma once

#include <string>
#include <vector>

#include "common/enums.h"

namespace plc::parser_util {

std::string trim(const std::string& s);
std::vector<std::string> split_ws(const std::string& s);
bool starts_with(const std::string& s, const std::string& prefix);
bool is_comment_or_empty(const std::string& s);
plc::Orient parse_orient(const std::string& token);

}  // namespace plc::parser_util

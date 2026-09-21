#include "parser/wts_parser.h"

#include <fstream>

#include "parser/parser_common.h"

namespace plc {

bool WtsParser::parse(const std::string& wts_path, PlacementDB* db, ValidationReport* report) const {
    if (db == nullptr || report == nullptr) {
        return false;
    }

    if (wts_path.empty()) {
        return true;
    }

    std::ifstream fin(wts_path);
    if (!fin) {
        report->add_warning("Failed to open wts file: " + wts_path);
        return true;
    }

    std::string line;
    while (std::getline(fin, line)) {
        if (parser_util::is_comment_or_empty(line)) {
            continue;
        }
        const auto tokens = parser_util::split_ws(line);
        if (tokens.size() < 2) {
            continue;
        }
        const double weight = std::stod(tokens[1]);
        const auto cell_iter = db->design.cell_name_to_id.find(tokens[0]);
        if (cell_iter != db->design.cell_name_to_id.end()) {
            db->cell(cell_iter->second).weight = weight;
            continue;
        }
        const auto net_iter = db->design.net_name_to_id.find(tokens[0]);
        if (net_iter != db->design.net_name_to_id.end()) {
            db->net(net_iter->second).weight = weight;
        }
    }

    return report->ok();
}

}  // namespace plc

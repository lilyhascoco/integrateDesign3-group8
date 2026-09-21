#include "parser/pl_parser.h"

#include <fstream>

#include "parser/parser_common.h"

namespace plc {

bool PlParser::parse(const std::string& pl_path, PlacementDB* db, ValidationReport* report) const {
    if (db == nullptr || report == nullptr) {
        return false;
    }

    if (pl_path.empty()) {
        report->add_error("Failed to parse .pl because file path is empty");
        return false;
    }

    std::ifstream fin(pl_path);
    if (!fin) {
        report->add_error("Failed to open pl file: " + pl_path);
        return false;
    }

    std::string line;
    while (std::getline(fin, line)) {
        if (parser_util::is_comment_or_empty(line)) {
            continue;
        }

        const auto tokens = parser_util::split_ws(line);
        if (tokens.size() < 3) {
            continue;
        }

        const auto iter = db->design.cell_name_to_id.find(tokens[0]);
        if (iter == db->design.cell_name_to_id.end()) {
            report->add_error("Unknown node in pl file: " + tokens[0]);
            continue;
        }

        Cell& cell = db->cell(iter->second);
        cell.init_x = std::stod(tokens[1]);
        cell.init_y = std::stod(tokens[2]);

        for (std::size_t i = 3; i < tokens.size(); ++i) {
            const std::string& token = tokens[i];
            if (token == ":") {
                continue;
            }
            const Orient orient = parser_util::parse_orient(token);
            if (orient != Orient::kUnknown) {
                cell.orient = orient;
                continue;
            }
            if (token == "/FIXED") {
                cell.place_status = PlaceStatus::kFixed;
            } else if (token == "/FIXED_NI") {
                cell.place_status = PlaceStatus::kFixedNI;
            } else if (token == "/PLACED" && !cell.is_fixed()) {
                cell.place_status = PlaceStatus::kMovable;
            }
        }
    }

    return report->ok();
}

}  // namespace plc

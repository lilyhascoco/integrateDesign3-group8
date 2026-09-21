#include "parser/scl_parser.h"

#include <fstream>
#include <limits>

#include "parser/parser_common.h"

namespace plc {

bool SclParser::parse(const std::string& scl_path, PlacementDB* db, ValidationReport* report) const {
    if (db == nullptr || report == nullptr) {
        return false;
    }

    std::ifstream fin(scl_path);
    if (!fin) {
        report->add_error("Failed to open scl file: " + scl_path);
        return false;
    }

    std::string line;
    bool in_row = false;
    Row row;
    while (std::getline(fin, line)) {
        const std::string trimmed = parser_util::trim(line);
        if (trimmed.empty() || parser_util::starts_with(trimmed, "#") || parser_util::starts_with(trimmed, "UCLA")) {
            continue;
        }
        if (trimmed == "CoreRow Horizontal") {
            in_row = true;
            row = Row{};
            row.id = static_cast<RowId>(db->design.rows.size());
            continue;
        }
        if (trimmed == "End") {
            if (in_row) {
                db->design.rows.push_back(row);
                in_row = false;
            }
            continue;
        }
        if (!in_row) {
            continue;
        }

        const auto pos = trimmed.find(':');
        if (pos == std::string::npos) {
            continue;
        }
        const std::string key = parser_util::trim(trimmed.substr(0, pos));
        const std::string value = parser_util::trim(trimmed.substr(pos + 1));
        const auto tokens = parser_util::split_ws(value);
        if (key == "Coordinate" && !tokens.empty()) {
            row.origin_y = std::stod(tokens[0]);
        } else if (key == "Height" && !tokens.empty()) {
            row.height = std::stod(tokens[0]);
        } else if (key == "Sitewidth" && !tokens.empty()) {
            row.site_width = std::stod(tokens[0]);
        } else if (key == "Sitespacing" && !tokens.empty()) {
            row.site_spacing = std::stod(tokens[0]);
        } else if (key == "Siteorient" && !tokens.empty()) {
            row.orient = parser_util::parse_orient(tokens[0]);
        } else if (key == "Sitesymmetry" && !tokens.empty()) {
            row.symmetry = tokens[0];
        } else if (key == "SubrowOrigin" && tokens.size() >= 4) {
            row.origin_x = std::stod(tokens[0]);
            row.num_sites = std::stoi(tokens.back());
        }
    }

    update_core_from_rows(db);
    if (db->design.rows.empty()) {
        report->add_warning("No rows parsed from scl file");
    }
    return report->ok();
}

void SclParser::update_core_from_rows(PlacementDB* db) const {
    if (db->design.rows.empty()) {
        return;
    }

    Coord min_x = std::numeric_limits<Coord>::max();
    Coord min_y = std::numeric_limits<Coord>::max();
    Coord max_x = std::numeric_limits<Coord>::lowest();
    Coord max_y = std::numeric_limits<Coord>::lowest();

    for (const auto& row : db->design.rows) {
        min_x = std::min(min_x, row.x_min());
        min_y = std::min(min_y, row.origin_y);
        max_x = std::max(max_x, row.x_max());
        max_y = std::max(max_y, row.origin_y + row.height);
    }

    db->design.core.lx = min_x;
    db->design.core.ly = min_y;
    db->design.core.ux = max_x;
    db->design.core.uy = max_y;
    db->design.core.default_row_height = db->design.rows.front().height;
    db->design.core.default_site_width = db->design.rows.front().site_width;
}

}  // namespace plc

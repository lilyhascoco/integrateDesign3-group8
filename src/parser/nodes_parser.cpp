#include "parser/nodes_parser.h"

#include <fstream>

#include "parser/parser_common.h"

namespace plc {

bool NodesParser::parse(const std::string& nodes_path, PlacementDB* db, ValidationReport* report) const {
    if (db == nullptr || report == nullptr) {
        return false;
    }

    std::ifstream fin(nodes_path);
    if (!fin) {
        report->add_error("Failed to open nodes file: " + nodes_path);
        return false;
    }

    int header_num_nodes = -1;
    int header_num_terminals = -1;
    int parsed_terminals = 0;
    std::string line;
    while (std::getline(fin, line)) {
        if (parser_util::is_comment_or_empty(line)) {
            continue;
        }
        const auto tokens = parser_util::split_ws(line);
        if (tokens.empty()) {
            continue;
        }
        if (tokens.size() >= 2 && tokens[0] == "NumNodes") {
            header_num_nodes = std::stoi(tokens.back());
            continue;
        }
        if (tokens.size() >= 2 && tokens[0] == "NumTerminals") {
            header_num_terminals = std::stoi(tokens.back());
            continue;
        }
        if (tokens.size() < 3) {
            continue;
        }

        const bool is_terminal = (tokens.size() >= 4 && tokens[3] == "terminal");
        const bool is_terminal_ni = (tokens.size() >= 4 && tokens[3] == "terminal_NI");

        Cell cell;
        cell.id = static_cast<CellId>(db->design.cells.size());
        cell.name = tokens[0];
        cell.width = std::stod(tokens[1]);
        cell.height = std::stod(tokens[2]);
        if (is_terminal_ni) {
            cell.type = CellType::kTerminalNI;
        } else if (is_terminal) {
            cell.type = CellType::kTerminal;
        } else {
            cell.type = infer_cell_type(cell.width, cell.height, db->design.core.default_row_height);
        }
        cell.place_status = is_terminal ? PlaceStatus::kFixed
                                        : (is_terminal_ni ? PlaceStatus::kFixedNI : PlaceStatus::kMovable);

        if (is_terminal || is_terminal_ni) {
            ++parsed_terminals;
        }

        db->design.cell_name_to_id[cell.name] = cell.id;
        db->design.cells.push_back(cell);
    }

    if (header_num_nodes >= 0 && static_cast<int>(db->design.cells.size()) != header_num_nodes) {
        report->add_warning("NumNodes header does not match parsed node count");
    }
    if (header_num_terminals >= 0 && parsed_terminals != header_num_terminals) {
        report->add_warning("NumTerminals header does not match parsed terminal count");
    }

    return report->ok();
}

CellType NodesParser::infer_cell_type(Coord width, Coord height, Coord row_height) const {
    if (row_height > 0.0 && height > row_height * 2.0) {
        return CellType::kMacro;
    }
    return CellType::kStdCell;
}

}  // namespace plc

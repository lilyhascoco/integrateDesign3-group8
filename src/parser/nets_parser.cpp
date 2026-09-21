#include "parser/nets_parser.h"

#include <fstream>

#include "parser/parser_common.h"

namespace plc {

bool NetsParser::parse(const std::string& nets_path, PlacementDB* db, ValidationReport* report) const {
    if (db == nullptr || report == nullptr) {
        return false;
    }

    std::ifstream fin(nets_path);
    if (!fin) {
        report->add_error("Failed to open nets file: " + nets_path);
        return false;
    }

    int header_num_nets = -1;
    int header_num_pins = -1;
    std::string line;
    while (std::getline(fin, line)) {
        if (parser_util::is_comment_or_empty(line)) {
            continue;
        }

        const auto tokens = parser_util::split_ws(line);
        if (tokens.empty()) {
            continue;
        }

        if (tokens.size() >= 2 && tokens[0] == "NumNets") {
            header_num_nets = std::stoi(tokens.back());
            continue;
        }
        if (tokens.size() >= 2 && tokens[0] == "NumPins") {
            header_num_pins = std::stoi(tokens.back());
            continue;
        }

        if (tokens.size() >= 3 && tokens[0] == "NetDegree") {
            const int degree = std::stoi(tokens[2]);
            Net net;
            net.id = static_cast<NetId>(db->design.nets.size());
            net.name = (tokens.size() >= 4) ? tokens[3] : ("net_" + std::to_string(net.id));
            db->design.net_name_to_id[net.name] = net.id;

            for (int pin_index = 0; pin_index < degree; ++pin_index) {
                std::string pin_line;
                if (!std::getline(fin, pin_line)) {
                    report->add_error("Unexpected end of file while reading pins for " + net.name);
                    break;
                }
                if (parser_util::is_comment_or_empty(pin_line)) {
                    --pin_index;
                    continue;
                }

                const auto pin_tokens = parser_util::split_ws(pin_line);
                if (pin_tokens.size() < 2) {
                    report->add_warning("Skip malformed pin line in net " + net.name);
                    continue;
                }
                const auto cell_iter = db->design.cell_name_to_id.find(pin_tokens[0]);
                if (cell_iter == db->design.cell_name_to_id.end()) {
                    report->add_error("Unknown node in nets file: " + pin_tokens[0]);
                    continue;
                }

                Pin pin;
                pin.id = static_cast<PinId>(db->design.pins.size());
                pin.cell_id = cell_iter->second;
                pin.net_id = net.id;
                pin.name = net.name + "_pin_" + std::to_string(pin.id);
                pin.direction = parse_pin_direction(pin_tokens[1]);

                // Common BookShelf pin format: name dir : offset_x offset_y
                if (pin_tokens.size() >= 5 && pin_tokens[2] == ":") {
                    pin.offset_x = std::stod(pin_tokens[3]);
                    pin.offset_y = std::stod(pin_tokens[4]);
                } else if (pin_tokens.size() >= 4) {
                    pin.offset_x = std::stod(pin_tokens[pin_tokens.size() - 2]);
                    pin.offset_y = std::stod(pin_tokens[pin_tokens.size() - 1]);
                }

                db->design.pins.push_back(pin);
                db->cell(pin.cell_id).pin_ids.push_back(pin.id);
                net.pin_ids.push_back(pin.id);
            }

            db->design.nets.push_back(net);
        }
    }

    if (header_num_nets >= 0 && static_cast<int>(db->design.nets.size()) != header_num_nets) {
        report->add_warning("NumNets header does not match parsed net count");
    }
    if (header_num_pins >= 0 && static_cast<int>(db->design.pins.size()) != header_num_pins) {
        report->add_warning("NumPins header does not match parsed pin count");
    }
    return report->ok();
}

PinDirection NetsParser::parse_pin_direction(const std::string& token) const {
    if (token == "I") return PinDirection::kInput;
    if (token == "O") return PinDirection::kOutput;
    if (token == "B") return PinDirection::kInout;
    return PinDirection::kUnknown;
}

}  // namespace plc

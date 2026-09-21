#include "parser/aux_parser.h"

#include <filesystem>
#include <fstream>
#include <sstream>

#include "parser/parser_common.h"

namespace fs = std::filesystem;

namespace plc {

bool AuxParser::parse(const std::string& aux_path, BookshelfInputSpec* spec, ValidationReport* report) const {
    if (spec == nullptr || report == nullptr) {
        return false;
    }

    std::ifstream fin(aux_path);
    if (!fin) {
        report->add_error("Failed to open aux file: " + aux_path);
        return false;
    }

    fs::path aux_fs_path(aux_path);
    spec->aux_path = aux_path;
    spec->benchmark_name = aux_fs_path.stem().string();
    spec->benchmark_dir = aux_fs_path.parent_path().string();

    std::string line;
    while (std::getline(fin, line)) {
        const std::string trimmed = parser_util::trim(line);
        if (trimmed.empty() || parser_util::starts_with(trimmed, "#")) {
            continue;
        }
        const auto pos = trimmed.find(':');
        const std::string rhs = (pos == std::string::npos) ? trimmed : trimmed.substr(pos + 1);
        std::istringstream iss(rhs);
        std::string token;
        while (iss >> token) {
            fs::path file_path = fs::path(spec->benchmark_dir) / token;
            const std::string ext = fs::path(token).extension().string();
            if (ext == ".nodes") {
                spec->nodes_path = file_path.string();
            } else if (ext == ".nets") {
                spec->nets_path = file_path.string();
            } else if (ext == ".pl") {
                spec->pl_path = file_path.string();
            } else if (ext == ".scl") {
                spec->scl_path = file_path.string();
            } else if (ext == ".wts") {
                spec->wts_path = file_path.string();
            }
        }
    }

    if (spec->nodes_path.empty()) report->add_error("Missing .nodes file in aux");
    if (spec->nets_path.empty()) report->add_error("Missing .nets file in aux");
    if (spec->pl_path.empty()) report->add_error("Missing .pl file in aux");
    if (spec->scl_path.empty()) report->add_error("Missing .scl file in aux");

    return report->ok();
}

}  // namespace plc

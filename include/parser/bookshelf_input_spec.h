#pragma once

#include <string>

namespace plc {

struct BookshelfInputSpec {
    std::string aux_path;
    std::string nodes_path;
    std::string nets_path;
    std::string pl_path;
    std::string scl_path;
    std::string wts_path;
    std::string benchmark_name;
    std::string benchmark_dir;
};

}  // namespace plc

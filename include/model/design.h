#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "model/cell.h"
#include "model/net.h"
#include "model/pin.h"
#include "model/row.h"

namespace plc {

struct DesignInfo {
    std::string design_name;
    std::string benchmark_dir;
    std::string aux_path;
    std::string scl_path;
    std::string nodes_path;
    std::string nets_path;
    std::string pl_path;
    std::string wts_path;
};

struct DesignStorage {
    DesignInfo info;
    CoreRegion core;
    std::vector<Row> rows;
    std::vector<Cell> cells;
    std::vector<Pin> pins;
    std::vector<Net> nets;
    std::unordered_map<std::string, CellId> cell_name_to_id;
    std::unordered_map<std::string, NetId> net_name_to_id;
};

}  // namespace plc

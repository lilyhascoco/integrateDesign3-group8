#pragma once

#include <string>

#include "database/placement_db.h"
#include "database/validation_report.h"

namespace plc {

class NodesParser {
public:
    bool parse(const std::string& nodes_path, PlacementDB* db, ValidationReport* report) const;

private:
    CellType infer_cell_type(Coord width, Coord height, Coord row_height) const;
};

}  // namespace plc

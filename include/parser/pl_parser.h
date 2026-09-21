#pragma once

#include <string>

#include "database/placement_db.h"
#include "database/validation_report.h"

namespace plc {

class PlParser {
public:
    bool parse(const std::string& pl_path, PlacementDB* db, ValidationReport* report) const;
};

}  // namespace plc

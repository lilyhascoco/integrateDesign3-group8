#pragma once

#include <string>

#include "database/placement_db.h"
#include "database/validation_report.h"

namespace plc {

class WtsParser {
public:
    bool parse(const std::string& wts_path, PlacementDB* db, ValidationReport* report) const;
};

}  // namespace plc

#pragma once

#include <string>

#include "database/placement_db.h"
#include "database/validation_report.h"

namespace plc {

class SclParser {
public:
    bool parse(const std::string& scl_path, PlacementDB* db, ValidationReport* report) const;

private:
    void update_core_from_rows(PlacementDB* db) const;
};

}  // namespace plc

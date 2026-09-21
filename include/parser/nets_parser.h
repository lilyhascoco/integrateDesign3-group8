#pragma once

#include <string>

#include "database/placement_db.h"
#include "database/validation_report.h"

namespace plc {

class NetsParser {
public:
    bool parse(const std::string& nets_path, PlacementDB* db, ValidationReport* report) const;

private:
    PinDirection parse_pin_direction(const std::string& token) const;
};

}  // namespace plc

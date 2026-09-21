#pragma once

#include "database/placement_db.h"

namespace plc {

struct AreaStats {
    double total_cell_area = 0.0;
    double movable_area = 0.0;
    double fixed_area = 0.0;
    double core_area = 0.0;
    double utilization = 0.0;
};

class AreaStatsCalculator {
public:
    AreaStats compute(const PlacementDB& db) const;
};

}  // namespace plc

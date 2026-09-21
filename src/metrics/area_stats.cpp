#include "metrics/area_stats.h"

namespace plc {

AreaStats AreaStatsCalculator::compute(const PlacementDB& db) const {
    AreaStats stats;
    stats.total_cell_area = db.total_cell_area();
    stats.movable_area = db.movable_cell_area();
    stats.fixed_area = stats.total_cell_area - stats.movable_area;
    stats.core_area = (db.design.core.ux - db.design.core.lx) * (db.design.core.uy - db.design.core.ly);
    if (stats.core_area > 0.0) {
        stats.utilization = stats.movable_area / stats.core_area;
    }
    return stats;
}

}  // namespace plc

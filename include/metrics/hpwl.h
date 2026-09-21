#pragma once

#include "database/placement_db.h"

namespace plc {

class HpwlEvaluator {
public:
    double compute_total_hpwl(const PlacementDB& db, const PlacementState& state) const;
    double compute_net_hpwl(const PlacementDB& db, const PlacementState& state, NetId net_id) const;
};

}  // namespace plc

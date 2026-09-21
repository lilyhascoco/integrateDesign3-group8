#pragma once

#include "database/placement_db.h"

namespace plc {

class GlobalPlacement {
public:
    virtual ~GlobalPlacement() = default;
    virtual void initialize(const PlacementDB& db, PlacementState* state) = 0;
    virtual void run_iteration() = 0;
    virtual bool converged() const = 0;
};

}  // namespace plc

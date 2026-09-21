#pragma once

#include "database/placement_db.h"

namespace plc {

class InitialPlacement {
public:
    virtual ~InitialPlacement() = default;
    virtual void run(const PlacementDB& db, PlacementState* state) = 0;
};

}  // namespace plc

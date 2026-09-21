#pragma once

#include <string>

#include "common/enums.h"
#include "common/types.h"

namespace plc {

struct Pin {
    PinId id = 0;
    std::string name;
    CellId cell_id = 0;
    NetId net_id = 0;
    Coord offset_x = 0.0;
    Coord offset_y = 0.0;
    PinDirection direction = PinDirection::kUnknown;
};

}  // namespace plc

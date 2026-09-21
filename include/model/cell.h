#pragma once

#include <string>
#include <vector>

#include "common/enums.h"
#include "common/types.h"

namespace plc {

struct Cell {
    CellId id = 0;
    std::string name;

    Coord width = 0.0;
    Coord height = 0.0;
    Coord init_x = 0.0;
    Coord init_y = 0.0;

    CellType type = CellType::kUnknown;
    PlaceStatus place_status = PlaceStatus::kUnplaced;
    Orient orient = Orient::kUnknown;
    double weight = 1.0;

    std::vector<PinId> pin_ids;

    Coord area() const { return width * height; }
    bool is_movable() const { return place_status == PlaceStatus::kMovable; }
    bool is_fixed() const {
        return place_status == PlaceStatus::kFixed || place_status == PlaceStatus::kFixedNI;
    }
};

}  // namespace plc

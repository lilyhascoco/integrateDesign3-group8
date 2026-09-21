#pragma once

#include "common/types.h"

namespace plc {

struct BBox {
    Coord lx = 0.0;
    Coord ly = 0.0;
    Coord ux = 0.0;
    Coord uy = 0.0;

    Coord width() const { return ux - lx; }
    Coord height() const { return uy - ly; }
    Coord area() const { return width() * height(); }
};

}  // namespace plc

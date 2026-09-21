#pragma once

#include <string>

#include "common/enums.h"
#include "common/types.h"

namespace plc {

struct Row {
    RowId id = 0;
    Coord origin_x = 0.0;
    Coord origin_y = 0.0;
    Coord height = 0.0;
    Coord site_width = 0.0;
    Coord site_spacing = 0.0;
    int num_sites = 0;
    Orient orient = Orient::kN;
    std::string symmetry;

    Coord x_min() const { return origin_x; }
    Coord x_max() const { return origin_x + static_cast<Coord>(num_sites) * site_spacing; }
};

struct CoreRegion {
    Coord lx = 0.0;
    Coord ly = 0.0;
    Coord ux = 0.0;
    Coord uy = 0.0;
    Coord default_row_height = 0.0;
    Coord default_site_width = 0.0;
};

}  // namespace plc

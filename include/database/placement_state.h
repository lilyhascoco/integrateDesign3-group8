#pragma once

#include <cstddef>
#include <vector>

#include "common/enums.h"
#include "common/types.h"

namespace plc {

class PlacementState {
public:
    void resize(std::size_t num_cells);
    std::size_t size() const { return x_.size(); }

    Coord x(CellId id) const;
    Coord y(CellId id) const;
    void set_pos(CellId id, Coord x, Coord y);

    Orient orient(CellId id) const;
    void set_orient(CellId id, Orient orient);

private:
    std::vector<Coord> x_;
    std::vector<Coord> y_;
    std::vector<Orient> orient_;
};

}  // namespace plc

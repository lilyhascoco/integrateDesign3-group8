#pragma once

#include <cstddef>

#include "database/placement_state.h"
#include "model/design.h"

namespace plc {

class PlacementDB {
public:
    DesignStorage design;
    PlacementState init_state;

    const Cell& cell(CellId id) const;
    Cell& cell(CellId id);
    const Net& net(NetId id) const;
    Net& net(NetId id);
    const Pin& pin(PinId id) const;
    const Row& row(RowId id) const;

    std::size_t num_cells() const;
    std::size_t num_nets() const;
    std::size_t num_pins() const;
    std::size_t num_rows() const;

    std::size_t num_movable_cells() const;
    std::size_t num_fixed_cells() const;
    double total_cell_area() const;
    double movable_cell_area() const;
};

}  // namespace plc

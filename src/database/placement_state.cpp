#include "database/placement_state.h"

#include <stdexcept>

namespace plc {

void PlacementState::resize(std::size_t num_cells) {
    x_.assign(num_cells, 0.0);
    y_.assign(num_cells, 0.0);
    orient_.assign(num_cells, Orient::kUnknown);
}

Coord PlacementState::x(CellId id) const {
    return x_.at(id);
}

Coord PlacementState::y(CellId id) const {
    return y_.at(id);
}

void PlacementState::set_pos(CellId id, Coord x, Coord y) {
    if (id >= x_.size() || id >= y_.size()) {
        throw std::out_of_range("cell id out of range in PlacementState::set_pos");
    }
    x_[id] = x;
    y_[id] = y;
}

Orient PlacementState::orient(CellId id) const {
    return orient_.at(id);
}

void PlacementState::set_orient(CellId id, Orient orient) {
    if (id >= orient_.size()) {
        throw std::out_of_range("cell id out of range in PlacementState::set_orient");
    }
    orient_[id] = orient;
}

}  // namespace plc

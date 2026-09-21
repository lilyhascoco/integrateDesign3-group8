#include "database/placement_db.h"

#include <numeric>

namespace plc {

const Cell& PlacementDB::cell(CellId id) const { return design.cells.at(id); }
Cell& PlacementDB::cell(CellId id) { return design.cells.at(id); }
const Net& PlacementDB::net(NetId id) const { return design.nets.at(id); }
Net& PlacementDB::net(NetId id) { return design.nets.at(id); }
const Pin& PlacementDB::pin(PinId id) const { return design.pins.at(id); }
const Row& PlacementDB::row(RowId id) const { return design.rows.at(id); }

std::size_t PlacementDB::num_cells() const { return design.cells.size(); }
std::size_t PlacementDB::num_nets() const { return design.nets.size(); }
std::size_t PlacementDB::num_pins() const { return design.pins.size(); }
std::size_t PlacementDB::num_rows() const { return design.rows.size(); }

std::size_t PlacementDB::num_movable_cells() const {
    std::size_t count = 0;
    for (const auto& c : design.cells) {
        if (c.is_movable()) {
            ++count;
        }
    }
    return count;
}

std::size_t PlacementDB::num_fixed_cells() const {
    std::size_t count = 0;
    for (const auto& c : design.cells) {
        if (c.is_fixed()) {
            ++count;
        }
    }
    return count;
}

double PlacementDB::total_cell_area() const {
    double area = 0.0;
    for (const auto& c : design.cells) {
        area += c.area();
    }
    return area;
}

double PlacementDB::movable_cell_area() const {
    double area = 0.0;
    for (const auto& c : design.cells) {
        if (c.is_movable()) {
            area += c.area();
        }
    }
    return area;
}

}  // namespace plc

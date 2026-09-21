#include "metrics/hpwl.h"

#include <limits>

namespace plc {

double HpwlEvaluator::compute_total_hpwl(const PlacementDB& db, const PlacementState& state) const {
    double total = 0.0;
    for (NetId net_id = 0; net_id < db.num_nets(); ++net_id) {
        total += compute_net_hpwl(db, state, net_id);
    }
    return total;
}

double HpwlEvaluator::compute_net_hpwl(const PlacementDB& db, const PlacementState& state, NetId net_id) const {
    const auto& net = db.net(net_id);
    if (net.pin_ids.size() <= 1) {
        return 0.0;
    }

    Coord min_x = std::numeric_limits<Coord>::max();
    Coord min_y = std::numeric_limits<Coord>::max();
    Coord max_x = std::numeric_limits<Coord>::lowest();
    Coord max_y = std::numeric_limits<Coord>::lowest();

    for (PinId pin_id : net.pin_ids) {
        const auto& pin = db.pin(pin_id);
        const Coord x = state.x(pin.cell_id) + pin.offset_x;
        const Coord y = state.y(pin.cell_id) + pin.offset_y;
        min_x = std::min(min_x, x);
        min_y = std::min(min_y, y);
        max_x = std::max(max_x, x);
        max_y = std::max(max_y, y);
    }

    return static_cast<double>((max_x - min_x) + (max_y - min_y)) * net.weight;
}

}  // namespace plc

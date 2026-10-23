#include "algorithms/local_sa_placer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>

#include "metrics/hpwl.h"

namespace plc {

namespace {

std::size_t bin_index(double center_x,
                      double center_y,
                      const CoreRegion& core,
                      std::size_t bins_x,
                      std::size_t bins_y) {
    const double bin_width = (core.ux - core.lx) / static_cast<double>(bins_x);
    const double bin_height = (core.uy - core.ly) / static_cast<double>(bins_y);

    const auto x = static_cast<std::size_t>(
        std::floor((center_x - core.lx) / bin_width));
    const auto y = static_cast<std::size_t>(
        std::floor((center_y - core.ly) / bin_height));
    const std::size_t clamped_x = std::min(x, bins_x - 1);
    const std::size_t clamped_y = std::min(y, bins_y - 1);
    return clamped_y * bins_x + clamped_x;
}

bool center_is_legal(const Cell& cell, double center_x, double center_y, const CoreRegion& core) {
    return center_x - cell.width * 0.5 >= core.lx &&
           center_x + cell.width * 0.5 <= core.ux &&
           center_y - cell.height * 0.5 >= core.ly &&
           center_y + cell.height * 0.5 <= core.uy;
}

}  // namespace

LocalSAPlacer::LocalSAPlacer(LocalSAConfig config, unsigned seed)
    : config_(config), rng_(seed) {
    config_.bins_x = std::max<std::size_t>(1, config_.bins_x);
    config_.bins_y = std::max<std::size_t>(1, config_.bins_y);
    config_.cooling_interval = std::max<std::size_t>(1, config_.cooling_interval);
    config_.initial_temperature = std::max(0.0, config_.initial_temperature);
    config_.cooling_factor = std::clamp(config_.cooling_factor, 0.0, 1.0);
    config_.minimum_temperature = std::max(0.0, config_.minimum_temperature);
}

void LocalSAPlacer::run(const PlacementDB& db, PlacementState* state) {
    stats_ = LocalSAStats{};
    if (state == nullptr || db.num_cells() == 0) {
        return;
    }

    const auto& core = db.design.core;
    const double core_width = core.ux - core.lx;
    const double core_height = core.uy - core.ly;
    if (core_width <= 0.0 || core_height <= 0.0) {
        return;
    }

    // Build bins from the existing placement. PlacementState stores lower-left
    // coordinates, so the center used by the binning rule is origin + size / 2.
    std::vector<std::vector<CellId>> bins(config_.bins_x * config_.bins_y);
    for (const auto& cell : db.design.cells) {
        if (!cell.is_movable()) {
            continue;
        }

        const double center_x = state->x(cell.id) + cell.width * 0.5;
        const double center_y = state->y(cell.id) + cell.height * 0.5;
        if (!center_is_legal(cell, center_x, center_y, core)) {
            ++stats_.skipped_moves;
            continue;
        }
        bins[bin_index(center_x, center_y, core, config_.bins_x, config_.bins_y)].push_back(cell.id);
    }

    std::vector<std::size_t> selectable_bins;
    for (std::size_t index = 0; index < bins.size(); ++index) {
        if (bins[index].size() >= 2) {
            selectable_bins.push_back(index);
        }
    }
    if (selectable_bins.empty()) {
        stats_.final_hpwl = HpwlEvaluator{}.compute_total_hpwl(db, *state);
        stats_.initial_hpwl = stats_.final_hpwl;
        return;
    }

    HpwlEvaluator hpwl_evaluator;
    stats_.initial_hpwl = hpwl_evaluator.compute_total_hpwl(db, *state);
    double total_hpwl = stats_.initial_hpwl;
    double temperature = config_.initial_temperature;
    std::uniform_real_distribution<double> probability(0.0, 1.0);
    std::uniform_int_distribution<std::size_t> bin_picker(0, selectable_bins.size() - 1);

    for (std::size_t iteration = 0;
         iteration < config_.max_iterations && temperature >= config_.minimum_temperature;
         ++iteration) {
        const auto& bin = bins[selectable_bins[bin_picker(rng_)]];
        std::uniform_int_distribution<std::size_t> cell_picker(0, bin.size() - 1);
        const CellId first_id = bin[cell_picker(rng_)];

        CellId second_id = first_id;
        while (second_id == first_id) {
            second_id = bin[cell_picker(rng_)];
        }

        const Cell& first = db.cell(first_id);
        const Cell& second = db.cell(second_id);
        const double first_x = state->x(first_id);
        const double first_y = state->y(first_id);
        const double second_x = state->x(second_id);
        const double second_y = state->y(second_id);

        const double first_center_x = first_x + first.width * 0.5;
        const double first_center_y = first_y + first.height * 0.5;
        const double second_center_x = second_x + second.width * 0.5;
        const double second_center_y = second_y + second.height * 0.5;

        // Swapping centers preserves locality. Different cell sizes can make a
        // swapped center illegal near a core boundary, so validate both cells.
        if (!center_is_legal(first, second_center_x, second_center_y, core) ||
            !center_is_legal(second, first_center_x, first_center_y, core)) {
            ++stats_.skipped_moves;
            continue;
        }

        // Only nets incident to either moved cell can change HPWL. Deduplication
        // avoids evaluating a shared net twice in the incremental delta.
        std::unordered_set<NetId> affected_nets;
        for (PinId pin_id : first.pin_ids) {
            affected_nets.insert(db.pin(pin_id).net_id);
        }
        for (PinId pin_id : second.pin_ids) {
            affected_nets.insert(db.pin(pin_id).net_id);
        }

        double old_local_hpwl = 0.0;
        for (NetId net_id : affected_nets) {
            old_local_hpwl += compute_net_hpwl(db, *state, net_id);
        }

        // PlacementState uses lower-left positions. Convert the destination
        // center back to an origin before evaluating the new pin coordinates.
        state->set_pos(first_id, second_center_x - first.width * 0.5,
                       second_center_y - first.height * 0.5);
        state->set_pos(second_id, first_center_x - second.width * 0.5,
                       first_center_y - second.height * 0.5);

        double new_local_hpwl = 0.0;
        for (NetId net_id : affected_nets) {
            new_local_hpwl += compute_net_hpwl(db, *state, net_id);
        }

        const double delta_cost = new_local_hpwl - old_local_hpwl;
        const bool accept = delta_cost <= 0.0 ||
                            (temperature > 0.0 &&
                             probability(rng_) < std::exp(-delta_cost / temperature));
        ++stats_.attempted_moves;

        if (accept) {
            total_hpwl += delta_cost;
            ++stats_.accepted_moves;
        } else {
            state->set_pos(first_id, first_x, first_y);
            state->set_pos(second_id, second_x, second_y);
            ++stats_.rejected_moves;
        }

        if ((iteration + 1) % config_.cooling_interval == 0) {
            temperature *= config_.cooling_factor;
        }
    }

    // Re-evaluate once after all incremental updates to avoid accumulated
    // floating-point error in the reported global result.
    stats_.final_hpwl = hpwl_evaluator.compute_total_hpwl(db, *state);
    stats_.final_temperature = temperature;
    (void)total_hpwl;
}

double LocalSAPlacer::compute_net_hpwl(const PlacementDB& db,
                                       const PlacementState& state,
                                       NetId net_id) const {
    return HpwlEvaluator{}.compute_net_hpwl(db, state, net_id);
}

}  // namespace plc

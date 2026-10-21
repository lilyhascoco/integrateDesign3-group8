#include "algorithms/random_placement.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace plc {

namespace {

double uniform_coord(std::mt19937* rng, double lower, double upper) {
    if (upper <= lower) {
        return lower;
    }
    std::uniform_real_distribution<double> dist(lower, upper);
    return dist(*rng);
}

}  // namespace

RandomPlacement::RandomPlacement(unsigned seed,
                                 std::size_t overlap_bins_x,
                                 std::size_t overlap_bins_y)
    : rng_(seed),
      overlap_bins_x_(std::max<std::size_t>(1, overlap_bins_x)),
      overlap_bins_y_(std::max<std::size_t>(1, overlap_bins_y)) {}

void RandomPlacement::run(const PlacementDB& db, PlacementState* state) {
    if (state == nullptr) {
        return;
    }

    stats_ = RandomPlacementStats{};
    state->resize(db.num_cells());

    const auto& core = db.design.core;
    const double core_width = core.ux - core.lx;
    const double core_height = core.uy - core.ly;

    for (const auto& cell : db.design.cells) {
        state->set_pos(cell.id, cell.init_x, cell.init_y);
        state->set_orient(cell.id, cell.orient);

        if (!cell.is_movable()) {
            continue;
        }

        ++stats_.movable_cells;
        stats_.movable_area += cell.area();

        if (core_width <= 0.0 || core_height <= 0.0 ||
            cell.width > core_width || cell.height > core_height) {
            ++stats_.skipped_cells;
            continue;
        }

        // Randomize the cell center inside the legal center range:
        // [core.lx + width / 2, core.ux - width / 2] and
        // [core.ly + height / 2, core.uy - height / 2].
        const double min_cx = core.lx + cell.width * 0.5;
        const double max_cx = core.ux - cell.width * 0.5;
        const double min_cy = core.ly + cell.height * 0.5;
        const double max_cy = core.uy - cell.height * 0.5;

        const double cx = uniform_coord(&rng_, min_cx, max_cx);
        const double cy = uniform_coord(&rng_, min_cy, max_cy);

        // The current PlacementState stores lower-left coordinates because
        // BookShelf .pl positions are parsed as lower-left positions. Convert
        // the sampled center back to lower-left before updating the state.
        const double llx = cx - cell.width * 0.5;
        const double lly = cy - cell.height * 0.5;
        state->set_pos(cell.id, llx, lly);

        ++stats_.placed_cells;
    }

    stats_.estimated_overlap_area = estimate_overlap_area(db, *state);
    if (stats_.movable_area > 0.0) {
        stats_.estimated_overlap_ratio = stats_.estimated_overlap_area / stats_.movable_area;
    }
}

double RandomPlacement::estimate_overlap_area(const PlacementDB& db, const PlacementState& state) const {
    const auto& core = db.design.core;
    const double core_width = core.ux - core.lx;
    const double core_height = core.uy - core.ly;
    if (core_width <= 0.0 || core_height <= 0.0) {
        return 0.0;
    }

    const double bin_width = core_width / static_cast<double>(overlap_bins_x_);
    const double bin_height = core_height / static_cast<double>(overlap_bins_y_);
    const double bin_capacity = bin_width * bin_height;
    std::vector<double> bin_demand(overlap_bins_x_ * overlap_bins_y_, 0.0);

    for (const auto& cell : db.design.cells) {
        const double cell_lx = std::max(state.x(cell.id), core.lx);
        const double cell_ly = std::max(state.y(cell.id), core.ly);
        const double cell_ux = std::min(state.x(cell.id) + cell.width, core.ux);
        const double cell_uy = std::min(state.y(cell.id) + cell.height, core.uy);
        if (cell_ux <= cell_lx || cell_uy <= cell_ly) {
            continue;
        }

        const std::size_t begin_x = static_cast<std::size_t>(
            std::floor((cell_lx - core.lx) / bin_width));
        const std::size_t begin_y = static_cast<std::size_t>(
            std::floor((cell_ly - core.ly) / bin_height));
        const std::size_t end_x = static_cast<std::size_t>(
            std::floor((cell_ux - core.lx) / bin_width));
        const std::size_t end_y = static_cast<std::size_t>(
            std::floor((cell_uy - core.ly) / bin_height));

        const std::size_t clamped_begin_x = std::min(begin_x, overlap_bins_x_ - 1);
        const std::size_t clamped_begin_y = std::min(begin_y, overlap_bins_y_ - 1);
        const std::size_t clamped_end_x = std::min(end_x, overlap_bins_x_ - 1);
        const std::size_t clamped_end_y = std::min(end_y, overlap_bins_y_ - 1);

        for (std::size_t by = clamped_begin_y; by <= clamped_end_y; ++by) {
            const double bin_ly = core.ly + static_cast<double>(by) * bin_height;
            const double bin_uy = bin_ly + bin_height;
            const double overlap_y = std::min(cell_uy, bin_uy) - std::max(cell_ly, bin_ly);
            if (overlap_y <= 0.0) {
                continue;
            }

            for (std::size_t bx = clamped_begin_x; bx <= clamped_end_x; ++bx) {
                const double bin_lx = core.lx + static_cast<double>(bx) * bin_width;
                const double bin_ux = bin_lx + bin_width;
                const double overlap_x = std::min(cell_ux, bin_ux) - std::max(cell_lx, bin_lx);
                if (overlap_x <= 0.0) {
                    continue;
                }

                bin_demand[by * overlap_bins_x_ + bx] += overlap_x * overlap_y;
            }
        }
    }

    double overflow = 0.0;
    for (double demand : bin_demand) {
        overflow += std::max(0.0, demand - bin_capacity);
    }
    return overflow;
}

}  // namespace plc

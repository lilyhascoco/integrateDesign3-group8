#pragma once

#include <cstddef>
#include <random>

#include "algorithms/initial_placement.h"

namespace plc {

struct RandomPlacementStats {
    std::size_t movable_cells = 0;
    std::size_t placed_cells = 0;
    std::size_t skipped_cells = 0;
    double movable_area = 0.0;
    double estimated_overlap_area = 0.0;
    double estimated_overlap_ratio = 0.0;
};

class RandomPlacement : public InitialPlacement {
public:
    explicit RandomPlacement(unsigned seed = std::random_device{}(),
                             std::size_t overlap_bins_x = 64,
                             std::size_t overlap_bins_y = 64);

    void run(const PlacementDB& db, PlacementState* state) override;

    const RandomPlacementStats& stats() const { return stats_; }

    // Bin overflow proxy shared by initializer comparisons; not exact geometric overlap.
    double estimate_overlap_area(const PlacementDB& db, const PlacementState& state) const;

private:
    std::mt19937 rng_;
    std::size_t overlap_bins_x_;
    std::size_t overlap_bins_y_;
    RandomPlacementStats stats_;
};

}  // namespace plc

#pragma once

#include <cstddef>
#include <random>

#include "algorithms/initial_placement.h"

namespace plc {

struct LocalSAConfig {
    std::size_t bins_x = 50;
    std::size_t bins_y = 50;
    std::size_t max_iterations = 100000;
    std::size_t cooling_interval = 500;
    double initial_temperature = 100.0;
    double cooling_factor = 0.95;
    double minimum_temperature = 0.1;
};

struct LocalSAStats {
    std::size_t attempted_moves = 0;
    std::size_t accepted_moves = 0;
    std::size_t rejected_moves = 0;
    std::size_t skipped_moves = 0;
    double initial_hpwl = 0.0;
    double final_hpwl = 0.0;
    double final_temperature = 0.0;
};

// Local simulated annealing using same-bin cell-center swaps.
class LocalSAPlacer : public InitialPlacement {
public:
    explicit LocalSAPlacer(LocalSAConfig config = {},
                           unsigned seed = std::random_device{}());

    void run(const PlacementDB& db, PlacementState* state) override;

    const LocalSAStats& stats() const { return stats_; }

private:
    double compute_net_hpwl(const PlacementDB& db,
                            const PlacementState& state,
                            NetId net_id) const;

    LocalSAConfig config_;
    std::mt19937 rng_;
    LocalSAStats stats_;
};

}  // namespace plc

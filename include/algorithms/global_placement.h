#pragma once

#include <cstddef>

#include "database/placement_db.h"

namespace plc {

// Parameters of the electrostatics-based global placer (ePlace style):
// minimize HPWL + lambda * density using Nesterov's accelerated gradient.
struct GlobalPlacementConfig {
    double target_density = 0.9;      // target utilization rho_t
    double target_overflow = 0.10;    // stopping criterion tau_min
    double lambda_reference = 3.5e5;  // Delta HPWL_ref in the lambda update
    std::size_t max_iterations = 200;
    std::size_t max_bins = 512;       // upper bound of the m x m bin grid
    unsigned seed = 2024u;
    bool verbose = true;
};

struct GlobalPlacementStats {
    std::size_t iterations = 0;
    std::size_t bins_x = 0;
    std::size_t bins_y = 0;
    std::size_t fillers = 0;
    double initial_hpwl = 0.0;
    double final_hpwl = 0.0;
    double initial_overflow = 0.0;
    double final_overflow = 0.0;
    double final_lambda = 0.0;
};

class GlobalPlacer {
public:
    explicit GlobalPlacer(GlobalPlacementConfig config = {});

    void run(const PlacementDB& db, PlacementState* state);

    const GlobalPlacementStats& stats() const { return stats_; }

private:
    GlobalPlacementConfig config_;
    GlobalPlacementStats stats_;
};

}  // namespace plc

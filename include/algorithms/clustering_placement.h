#pragma once
#include "algorithms/initial_placement.h"
namespace plc {
struct ClusteringConfig {
    std::size_t target_clusters = 0; // 0: ceil(movable / 32)
    std::size_t max_cluster_cells = 64;
    std::size_t clique_limit = 32;
    std::size_t coarse_iterations = 20;
};
struct ClusteringStats {
    std::size_t movable_cells = 0, placed_cells = 0, skipped_cells = 0;
    std::size_t clusters = 0, merges = 0, approximated_nets = 0;
};
class ClusteringPlacement : public InitialPlacement {
public:
    explicit ClusteringPlacement(ClusteringConfig config = {}) : config_(config) {}
    void run(const PlacementDB& db, PlacementState* state) override;
    const ClusteringStats& stats() const { return stats_; }
private:
    ClusteringConfig config_;
    ClusteringStats stats_;
};
} // namespace plc

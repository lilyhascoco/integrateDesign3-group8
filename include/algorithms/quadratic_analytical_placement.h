#pragma once

#include <cstddef>
#include <vector>

#include "algorithms/initial_placement.h"
#include "common/types.h"

namespace plc {

struct QuadraticConfig {
    // Outer reweighting iterations for w = 1 / distance.
    std::size_t max_outer_iterations = 5;
    // Maximum BiCGSTAB iterations per outer iteration.
    std::size_t solver_max_iterations = 700;
    // Relative infinity-norm residual target: ||b - A*x||inf / max(||b||inf, 1).
    double target_error = 1e-4;
    // Minimum distance protects the reweighting denominator.
    double min_distance = 1e-3;
    // Small legal-placement margin inside the core boundary.
    double eps = 1e-3;
};

struct QuadraticStats {
    std::size_t movable_cells = 0;
    std::size_t outer_iterations = 0;
    std::size_t total_solver_iterations = 0;
    double initial_hpwl = 0.0;
    double final_hpwl = 0.0;
    double final_x_error = 0.0;
    double final_y_error = 0.0;
};

// Row-major compressed sparse matrix for BiCGSTAB multiplication.
struct CSRMatrix {
    std::size_t n = 0;
    std::vector<double> values;
    std::vector<std::size_t> row_ptr;
    std::vector<std::size_t> col_idx;

    CSRMatrix() = default;
    explicit CSRMatrix(std::size_t order) : n(order), row_ptr(order + 1, 0) {}

    // y = A * x
    void multiply(const std::vector<double>& x, std::vector<double>* y) const;
};

// BiCGSTAB solver for the CSR linear system A * x = b.
double BiCGSTAB(const CSRMatrix& A,
                const std::vector<double>& b,
                std::vector<double>* x,
                double tolerance,
                std::size_t max_iterations,
                std::size_t* iterations_used = nullptr);

// Quadratic analytical placement using a CSR matrix and BiCGSTAB.
class QuadraticAnalyticalPlacement : public InitialPlacement {
public:
    explicit QuadraticAnalyticalPlacement(QuadraticConfig config = {});

    void run(const PlacementDB& db, PlacementState* state) override;

    const QuadraticStats& stats() const { return stats_; }

private:
    // Build movable-cell indices and append solver-only net-star indices.
    void build_movable_index_map(const PlacementDB& db);

    // Rebuild x/y systems from the current cell-center coordinates.
    void build_system(const PlacementDB& db,
                      const PlacementState& state,
                      const std::vector<double>& centers_x,
                      const std::vector<double>& centers_y,
                      CSRMatrix* A_x,
                      CSRMatrix* A_y,
                      std::vector<double>* b_x,
                      std::vector<double>* b_y,
                      std::vector<double>* new_centers_x,
                      std::vector<double>* new_centers_y) const;

    // Clamp real movable-cell centers to the core.
    void clip_to_core(const PlacementDB& db,
                      std::vector<double>* centers_x,
                      std::vector<double>* centers_y) const;

    // Write only real movable-cell centers back as lower-left positions.
    void write_centers_to_state(const PlacementDB& db,
                                const std::vector<double>& centers_x,
                                const std::vector<double>& centers_y,
                                PlacementState* state) const;

    // Compute an absolute pin coordinate from a cell-center coordinate.
    double pin_abs_x(const Cell& cell, const Pin& pin, double cell_center_x) const;
    double pin_abs_y(const Cell& cell, const Pin& pin, double cell_center_y) const;

    QuadraticConfig config_;
    QuadraticStats stats_;

    // Real movable-cell indices only; used by writeback and boundary clipping.
    std::vector<CellId> index_to_cell_;
    std::vector<std::size_t> cell_to_index_;

    // Solver-only star index per net, or max<size_t> when unused.
    std::vector<std::size_t> net_to_star_index_;
    std::size_t movable_count_ = 0;
};

}  // namespace plc

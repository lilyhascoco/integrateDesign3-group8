#pragma once

#include <cstddef>
#include <vector>

#include "algorithms/initial_placement.h"
#include "common/types.h"

namespace plc {

struct QuadraticConfig {
    // 外层迭代次数：权重 w = 1/distance 依赖当前位置，需要反复构建和求解
    std::size_t max_outer_iterations = 30;
    // 求解器最大迭代次数（BiCGSTAB 内部）
    std::size_t solver_max_iterations = 500;
    // 外层迭代收敛误差阈值
    double target_error = 1e-4;
    // 最小距离：防止 weight = 1/distance 在两点接近时爆炸
    double min_distance = 1e-3;
    // 一个小的安全边距，避免单元紧贴 core 边界
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

// Compressed Sparse Row 稀疏矩阵，行优先存储。
// 仅提供 BiCGSTAB 所需的矩阵-向量乘法能力。
struct CSRMatrix {
    std::size_t n = 0;                   // 矩阵阶数
    std::vector<double> values;          // 非零元素值（行优先）
    std::vector<std::size_t> row_ptr;    // row_ptr[n+1]，第 i 行非零元范围 [row_ptr[i], row_ptr[i+1])
    std::vector<std::size_t> col_idx;    // 每个 value 对应的列号

    CSRMatrix() = default;
    explicit CSRMatrix(std::size_t order) : n(order), row_ptr(order + 1, 0) {}

    // y = A * x
    void multiply(const std::vector<double>& x, std::vector<double>* y) const;
};

// BiCGSTAB (Biconjugate Gradient Stabilized) 迭代求解器。
// 求解 A*x = b，其中 A 是 CSR 格式的方阵。
// 返回值：最终残差（infinity norm）；若迭代次数用尽返回 last_error。
double BiCGSTAB(const CSRMatrix& A,
                const std::vector<double>& b,
                std::vector<double>* x,
                double tolerance,
                std::size_t max_iterations,
                std::size_t* iterations_used = nullptr);

// 二次方程解析布局（Quadratic Analytical Placement）。
// 通过构建线长的二次能量函数，最小化得到可移动单元中心的解析解。
// 不依赖 Eigen，使用纯 C++17 实现 CSR + BiCGSTAB，保证 Ubuntu 上可直接编译。
class QuadraticAnalyticalPlacement : public InitialPlacement {
public:
    explicit QuadraticAnalyticalPlacement(QuadraticConfig config = {});

    void run(const PlacementDB& db, PlacementState* state) override;

    const QuadraticStats& stats() const { return stats_; }

private:
    // 为可移动单元建立一个“在矩阵中的索引”到 CellId 的映射，
    // 固定单元（terminal）不进入矩阵，而是作为常数写入 b 向量。
    void build_movable_index_map(const PlacementDB& db);

    // 在某次迭代中，基于当前 cell center 重新计算 pin 的绝对坐标，
    // 然后构建 A 和 b。A、b 都要按 x/y 两个方向分别构建。
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

    // 将求得的单元中心裁剪回 core 区域
    void clip_to_core(const PlacementDB& db,
                      std::vector<double>* centers_x,
                      std::vector<double>* centers_y) const;

    // 将求解得到的 cell center 写回 PlacementState（由 center 转回 lower-left）
    void write_centers_to_state(const PlacementDB& db,
                                const std::vector<double>& centers_x,
                                const std::vector<double>& centers_y,
                                PlacementState* state) const;

    // 给定 cell center 坐标，计算指定 pin 的绝对坐标
    double pin_abs_x(const Cell& cell, const Pin& pin, double cell_center_x) const;
    double pin_abs_y(const Cell& cell, const Pin& pin, double cell_center_y) const;

    QuadraticConfig config_;
    QuadraticStats stats_;

    // 从 CellId 到矩阵索引（movable cells）
    std::vector<CellId> index_to_cell_;          // size = num_movable
    std::vector<std::size_t> cell_to_index_;     // size = num_cells, 固定单元存 num_movable
};

}  // namespace plc

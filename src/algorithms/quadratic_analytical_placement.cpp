#include <stdexcept>
#include "algorithms/quadratic_analytical_placement.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "metrics/hpwl.h"

namespace plc {

namespace {

// Jacobi 预条件：对每一行按对角线元素做归一化，让矩阵对角线统一为 1。
// 行缩放不改变方程的解 x，但能大幅改善 BiCGSTAB 的收敛速度。
void ApplyJacobiPreconditioner(CSRMatrix* A, std::vector<double>* b) {
    for (std::size_t i = 0; i < A->n; ++i) {
        double diag = 0.0;
        const std::size_t begin = A->row_ptr[i];
        const std::size_t end = A->row_ptr[i + 1];
        for (std::size_t k = begin; k < end; ++k) {
            if (A->col_idx[k] == i) {
                diag = A->values[k];
                break;
            }
        }
        if (std::fabs(diag) < 1e-300) diag = 1.0;
        const double inv = 1.0 / diag;
        for (std::size_t k = begin; k < end; ++k) {
            A->values[k] *= inv;
        }
        (*b)[i] *= inv;
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// CSR matrix – 行优先存储，只提供 A*x 矩阵向量乘法
// ---------------------------------------------------------------------------
void CSRMatrix::multiply(const std::vector<double>& x,
                         std::vector<double>* y) const {
    const std::size_t m = n;
    y->assign(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        double acc = 0.0;
        const std::size_t begin = row_ptr[i];
        const std::size_t end = row_ptr[i + 1];
        for (std::size_t k = begin; k < end; ++k) {
            acc += values[k] * x[col_idx[k]];
        }
        (*y)[i] = acc;
    }
}

// ---------------------------------------------------------------------------
// BiCGSTAB – 稳定双共轭梯度法，用于求解 Ax = b
// ---------------------------------------------------------------------------
double BiCGSTAB(const CSRMatrix& A,
                const std::vector<double>& b,
                std::vector<double>* x,
                double tolerance,
                std::size_t max_iterations,
                std::size_t* iterations_used) {
    const std::size_t n = b.size();
    if (n == 0 || A.n != n) {
        return 0.0;
    }

    std::vector<double> r(n), r_hat(n), p(n), v(n), s(n), t(n);

    // r = b - A*x0
    std::vector<double> ax(n);
    A.multiply(*x, &ax);
    for (std::size_t i = 0; i < n; ++i) {
        r[i] = b[i] - ax[i];
    }
    r_hat = r;
    double rhs_norm = 0.0;
    for (double value : b) {
        rhs_norm = std::max(rhs_norm, std::fabs(value));
    }
    rhs_norm = std::max(rhs_norm, 1.0);

    auto relative_infinity_norm = [rhs_norm](const std::vector<double>& residual) {
        double norm = 0.0;
        for (double value : residual) {
            norm = std::max(norm, std::fabs(value));
        }
        return norm / rhs_norm;
    };

    double rho_old = 1.0;
    double alpha = 1.0;
    double omega = 1.0;

    double r_dot_sq = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        r_dot_sq += r_hat[i] * r[i];
    }

    double last_error = relative_infinity_norm(r);
    tolerance = std::max(tolerance, 1e-16);
    max_iterations = std::max<std::size_t>(1, max_iterations);
    if (last_error < tolerance) {
        if (iterations_used) *iterations_used = 0;
        return last_error;
    }

    for (std::size_t iter = 0; iter < max_iterations; ++iter) {
        const double rho_new = r_dot_sq;
        if (std::fabs(rho_new) < 1e-300) {
            // 接近 breakdown，直接退出
            if (iterations_used) *iterations_used = iter;
            return last_error;
        }

        if (iter == 0) {
            for (std::size_t i = 0; i < n; ++i) {
                p[i] = r[i];
            }
        } else {
            const double beta = (rho_new / rho_old) * (alpha / omega);
            for (std::size_t i = 0; i < n; ++i) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
        }

        A.multiply(p, &v);

        double v_dot_rhat = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            v_dot_rhat += v[i] * r_hat[i];
        }
        if (std::fabs(v_dot_rhat) < 1e-300) {
            if (iterations_used) *iterations_used = iter;
            return last_error;
        }

        alpha = rho_new / v_dot_rhat;
        for (std::size_t i = 0; i < n; ++i) {
            s[i] = r[i] - alpha * v[i];
        }

        // 检查 s 的范数，若已经很小，直接用 x += alpha * p 并退出
        if (relative_infinity_norm(s) < tolerance) {
            for (std::size_t i = 0; i < n; ++i) {
                (*x)[i] += alpha * p[i];
            }
            if (iterations_used) *iterations_used = iter + 1;
            return relative_infinity_norm(s);
        }

        A.multiply(s, &t);

        double t_dot_sq = 0.0, t_norm_sq = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            t_dot_sq += t[i] * s[i];
            t_norm_sq += t[i] * t[i];
        }
        omega = (t_norm_sq < 1e-300) ? 0.0 : t_dot_sq / t_norm_sq;

        for (std::size_t i = 0; i < n; ++i) {
            (*x)[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }

        // Use a scale-independent residual so the same tolerance is meaningful
        // for both the sample and million-pin benchmarks.
        last_error = relative_infinity_norm(r);

        // 更新 r_dot_sq 供下一次迭代用
        r_dot_sq = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            r_dot_sq += r_hat[i] * r[i];
        }
        rho_old = rho_new;

        if (last_error < tolerance) {
            if (iterations_used) *iterations_used = iter + 1;
            return last_error;
        }

        if (std::fabs(omega) < 1e-300) {
            // omega ≈ 0，BiCGSTAB breakdown
            if (iterations_used) *iterations_used = iter + 1;
            return last_error;
        }
    }

    if (iterations_used) *iterations_used = max_iterations;
    return last_error;
}

// ---------------------------------------------------------------------------
// QuadraticAnalyticalPlacement
// ---------------------------------------------------------------------------
QuadraticAnalyticalPlacement::QuadraticAnalyticalPlacement(QuadraticConfig config)
    : config_(config) {
    config_.max_outer_iterations = std::max<std::size_t>(1, config_.max_outer_iterations);
    config_.solver_max_iterations = std::max<std::size_t>(1, config_.solver_max_iterations);
    config_.min_distance = std::max(0.0, config_.min_distance);
    config_.eps = std::max(0.0, config_.eps);
}

void QuadraticAnalyticalPlacement::build_movable_index_map(const PlacementDB& db) {
    index_to_cell_.clear();
    cell_to_index_.assign(db.num_cells(), std::numeric_limits<std::size_t>::max());
    index_to_cell_.reserve(db.num_cells());

    for (const auto& cell : db.design.cells) {
        if (cell.is_movable()) {
            const std::size_t idx = index_to_cell_.size();
            cell_to_index_[cell.id] = idx;
            index_to_cell_.push_back(cell.id);
        }
    }

    movable_count_ = index_to_cell_.size();
    net_to_star_index_.assign(db.num_nets(), std::numeric_limits<std::size_t>::max());
    std::size_t next_index = movable_count_;
    for (const auto& net : db.design.nets) {
        if (net.degree() < 3) {
            continue;
        }
        const bool touches_movable = std::any_of(
            net.pin_ids.begin(), net.pin_ids.end(), [&](PinId pin_id) {
                return db.cell(db.pin(pin_id).cell_id).is_movable();
            });
        if (touches_movable) {
            net_to_star_index_[net.id] = next_index++;
        }
    }
}

double QuadraticAnalyticalPlacement::pin_abs_x(const Cell& cell,
                                               const Pin& pin,
                                               double cell_center_x) const {
    // cell_center_x 已经是中心点，pin.offset_x 是相对于 cell center 的偏移
    return cell_center_x + pin.offset_x;
}

double QuadraticAnalyticalPlacement::pin_abs_y(const Cell& cell,
                                               const Pin& pin,
                                               double cell_center_y) const {
    return cell_center_y + pin.offset_y;
}

void QuadraticAnalyticalPlacement::build_system(const PlacementDB& db,
                                                  const PlacementState& state,
                                                  const std::vector<double>& centers_x,
                                                  const std::vector<double>& centers_y,
                                                  CSRMatrix* A_x,
                                                  CSRMatrix* A_y,
                                                  std::vector<double>* b_x,
                                                  std::vector<double>* b_y,
                                                  std::vector<double>* new_centers_x,
                                                  std::vector<double>* new_centers_y) const {
    const std::size_t n = centers_x.size();
    new_centers_x->assign(centers_x.begin(), centers_x.end());
    new_centers_y->assign(centers_y.begin(), centers_y.end());

    // Matrix indices [0, movable_count_) are real cells. Net-star indices are
    // appended afterwards and exist only while solving this linear system.
    std::vector<std::vector<std::pair<std::size_t, double>>> row_entries_x(n);
    std::vector<std::vector<std::pair<std::size_t, double>>> row_entries_y(n);
    b_x->assign(n, 0.0);
    b_y->assign(n, 0.0);

    auto add_movable_star_edge =
        [&](std::size_t cell_index, std::size_t star_index, double edge_weight,
            double pin_offset, std::vector<std::vector<std::pair<std::size_t, double>>>* rows,
            std::vector<double>* rhs) {
            (*rows)[cell_index].push_back({cell_index, edge_weight});
            (*rows)[cell_index].push_back({star_index, -edge_weight});
            (*rows)[star_index].push_back({cell_index, -edge_weight});
            (*rows)[star_index].push_back({star_index, edge_weight});
            (*rhs)[cell_index] += 2.0 * edge_weight * pin_offset;
            (*rhs)[star_index] -= 2.0 * edge_weight * pin_offset;
        };

    // A two-pin net is already a single O(1) edge. Keeping it direct avoids
    // adding a numerically unnecessary auxiliary variable and preserves the
    // prior two-pin system exactly.
    auto add_two_pin_edge = [&](const Pin& pin_a, const Cell& cell_a,
                                const Pin& pin_b, const Cell& cell_b,
                                double net_weight) {
        if (cell_a.id == cell_b.id) {
            return;
        }
        const bool a_movable = cell_a.is_movable();
        const bool b_movable = cell_b.is_movable();
        const double center_ax = a_movable ? centers_x[cell_to_index_[cell_a.id]]
                                           : state.x(cell_a.id) + cell_a.width * 0.5;
        const double center_ay = a_movable ? centers_y[cell_to_index_[cell_a.id]]
                                           : state.y(cell_a.id) + cell_a.height * 0.5;
        const double center_bx = b_movable ? centers_x[cell_to_index_[cell_b.id]]
                                           : state.x(cell_b.id) + cell_b.width * 0.5;
        const double center_by = b_movable ? centers_y[cell_to_index_[cell_b.id]]
                                           : state.y(cell_b.id) + cell_b.height * 0.5;
        const double pa_x = pin_abs_x(cell_a, pin_a, center_ax);
        const double pa_y = pin_abs_y(cell_a, pin_a, center_ay);
        const double pb_x = pin_abs_x(cell_b, pin_b, center_bx);
        const double pb_y = pin_abs_y(cell_b, pin_b, center_by);
        const double w_x = net_weight / std::max(std::fabs(pa_x - pb_x), config_.min_distance);
        const double w_y = net_weight / std::max(std::fabs(pa_y - pb_y), config_.min_distance);

        if (a_movable && b_movable) {
            const std::size_t i = cell_to_index_[cell_a.id];
            const std::size_t j = cell_to_index_[cell_b.id];
            for (auto* rows : {&row_entries_x, &row_entries_y}) {
                const double w = (rows == &row_entries_x) ? w_x : w_y;
                (*rows)[i].push_back({i, w});
                (*rows)[i].push_back({j, -w});
                (*rows)[j].push_back({i, -w});
                (*rows)[j].push_back({j, w});
            }
            (*b_x)[i] += 2.0 * w_x * (pin_a.offset_x - pin_b.offset_x);
            (*b_x)[j] += 2.0 * w_x * (pin_b.offset_x - pin_a.offset_x);
            (*b_y)[i] += 2.0 * w_y * (pin_a.offset_y - pin_b.offset_y);
            (*b_y)[j] += 2.0 * w_y * (pin_b.offset_y - pin_a.offset_y);
        } else if (a_movable) {
            const std::size_t i = cell_to_index_[cell_a.id];
            row_entries_x[i].push_back({i, w_x});
            row_entries_y[i].push_back({i, w_y});
            (*b_x)[i] += 2.0 * w_x * (pin_a.offset_x - pb_x);
            (*b_y)[i] += 2.0 * w_y * (pin_a.offset_y - pb_y);
        } else if (b_movable) {
            const std::size_t j = cell_to_index_[cell_b.id];
            row_entries_x[j].push_back({j, w_x});
            row_entries_y[j].push_back({j, w_y});
            (*b_x)[j] += 2.0 * w_x * (pin_b.offset_x - pa_x);
            (*b_y)[j] += 2.0 * w_y * (pin_b.offset_y - pa_y);
        }
    };

    for (const auto& net : db.design.nets) {
        const std::size_t star_index = net_to_star_index_[net.id];
        if (net.degree() == 2) {
            const Pin& pin_a = db.pin(net.pin_ids[0]);
            const Pin& pin_b = db.pin(net.pin_ids[1]);
            add_two_pin_edge(pin_a, db.cell(pin_a.cell_id), pin_b, db.cell(pin_b.cell_id),
                             net.weight);
            continue;
        }
        if (star_index == std::numeric_limits<std::size_t>::max()) {
            continue;
        }

        double min_x = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double min_y = std::numeric_limits<double>::max();
        double max_y = std::numeric_limits<double>::lowest();
        for (PinId pin_id : net.pin_ids) {
            const Pin& pin = db.pin(pin_id);
            const Cell& cell = db.cell(pin.cell_id);
            const bool movable = cell.is_movable();
            const double center_x = movable ? centers_x[cell_to_index_[cell.id]]
                                            : state.x(cell.id) + cell.width * 0.5;
            const double center_y = movable ? centers_y[cell_to_index_[cell.id]]
                                            : state.y(cell.id) + cell.height * 0.5;
            const double pin_x = pin_abs_x(cell, pin, center_x);
            const double pin_y = pin_abs_y(cell, pin, center_y);
            min_x = std::min(min_x, pin_x);
            max_x = std::max(max_x, pin_x);
            min_y = std::min(min_y, pin_y);
            max_y = std::max(max_y, pin_y);
        }

        // Eliminating a star center with d equal-weight spokes creates a
        // uniform clique whose pair weight is spoke_weight / d. For d == 2
        // this exactly reproduces the prior pairwise reweighted objective.
        const double degree = static_cast<double>(net.degree());
        const double edge_weight_x =
            degree * net.weight / std::max(max_x - min_x, config_.min_distance);
        const double edge_weight_y =
            degree * net.weight / std::max(max_y - min_y, config_.min_distance);

        for (PinId pin_id : net.pin_ids) {
            const Pin& pin = db.pin(pin_id);
            const Cell& cell = db.cell(pin.cell_id);
            if (cell.is_movable()) {
                const std::size_t cell_index = cell_to_index_[cell.id];
                add_movable_star_edge(cell_index, star_index, edge_weight_x, pin.offset_x,
                                       &row_entries_x, b_x);
                add_movable_star_edge(cell_index, star_index, edge_weight_y, pin.offset_y,
                                       &row_entries_y, b_y);
            } else {
                const double fixed_x = state.x(cell.id) + cell.width * 0.5 + pin.offset_x;
                const double fixed_y = state.y(cell.id) + cell.height * 0.5 + pin.offset_y;
                row_entries_x[star_index].push_back({star_index, edge_weight_x});
                row_entries_y[star_index].push_back({star_index, edge_weight_y});
                (*b_x)[star_index] -= 2.0 * edge_weight_x * fixed_x;
                (*b_y)[star_index] -= 2.0 * edge_weight_y * fixed_y;
            }
        }
    }

    // 将 row_entries_x/y 合并（同一 (i,j) 多次累加）并转成 CSR
    auto finalize_csr = [&](std::vector<std::vector<std::pair<std::size_t, double>>>& row_entries,
                            CSRMatrix* A) {
        A->n = n;
        A->row_ptr.assign(n + 1, 0);

        // 1) 对每行按列排序并累加重复项
        std::vector<std::vector<std::pair<std::size_t, double>>> merged(n);
        for (std::size_t i = 0; i < n; ++i) {
            auto& entries = row_entries[i];
            if (entries.empty()) {
                merged[i].push_back({i, 0.0});  // 对角线保证矩阵非奇异
                continue;
            }
            std::sort(entries.begin(), entries.end(),
                      [](const auto& l, const auto& r) { return l.first < r.first; });
            // 累加
            std::vector<std::pair<std::size_t, double>> acc;
            acc.reserve(entries.size());
            std::size_t prev_col = entries[0].first;
            double prev_val = entries[0].second;
            for (std::size_t k = 1; k < entries.size(); ++k) {
                if (entries[k].first == prev_col) {
                    prev_val += entries[k].second;
                } else {
                    acc.push_back({prev_col, prev_val});
                    prev_col = entries[k].first;
                    prev_val = entries[k].second;
                }
            }
            acc.push_back({prev_col, prev_val});
            // 保证对角元存在（即使值为 0 也放一个）
            bool has_diag = false;
            for (auto& p : acc) {
                if (p.first == i) {
                    has_diag = true;
                    if (std::fabs(p.second) < 1e-300) p.second = 1e-10;  // 避免对角线为 0
                    break;
                }
            }
            if (!has_diag) {
                // 插入对角线元素，保持列有序
                for (auto it = acc.begin(); it != acc.end(); ++it) {
                    if (it->first > i) {
                        acc.insert(it, {i, 1e-10});
                        has_diag = true;
                        break;
                    }
                }
                if (!has_diag) acc.push_back({i, 1e-10});
            }
            merged[i] = std::move(acc);
        }

        // 2) 填充 row_ptr / values / col_idx
        std::size_t total_nnz = 0;
        for (std::size_t i = 0; i < n; ++i) total_nnz += merged[i].size();

        A->row_ptr[0] = 0;
        A->values.reserve(total_nnz);
        A->col_idx.reserve(total_nnz);
        for (std::size_t i = 0; i < n; ++i) {
            for (auto& [col, val] : merged[i]) {
                A->col_idx.push_back(col);
                A->values.push_back(val);
            }
            A->row_ptr[i + 1] = A->col_idx.size();
        }
    };

    finalize_csr(row_entries_x, A_x);
    finalize_csr(row_entries_y, A_y);
}

void QuadraticAnalyticalPlacement::clip_to_core(const PlacementDB& db,
                                                 std::vector<double>* centers_x,
                                                 std::vector<double>* centers_y) const {
    const auto& core = db.design.core;
    const double eps = config_.eps;

    for (std::size_t i = 0; i < index_to_cell_.size(); ++i) {
        const Cell& cell = db.cell(index_to_cell_[i]);
        double& cx = (*centers_x)[i];
        double& cy = (*centers_y)[i];

        // An epsilon margin cannot exceed the available slack. For a cell
        // exactly as tall/wide as the core the only legal center has no margin.
        const double slack_x = core.ux - core.lx - cell.width;
        const double slack_y = core.uy - core.ly - cell.height;
        if (slack_x < 0 || slack_y < 0 || !std::isfinite(cx) || !std::isfinite(cy))
            throw std::runtime_error("Quadratic placement: invalid cell size or non-finite center");
        const double margin_x = std::min(std::max(0.0, eps), slack_x * 0.5);
        const double margin_y = std::min(std::max(0.0, eps), slack_y * 0.5);
        cx = std::clamp(cx, core.lx + cell.width * 0.5 + margin_x,
                            core.ux - cell.width * 0.5 - margin_x);
        cy = std::clamp(cy, core.ly + cell.height * 0.5 + margin_y,
                            core.uy - cell.height * 0.5 - margin_y);
    }
}

void QuadraticAnalyticalPlacement::write_centers_to_state(const PlacementDB& db,
                                                          const std::vector<double>& centers_x,
                                                          const std::vector<double>& centers_y,
                                                          PlacementState* state) const {
    for (std::size_t i = 0; i < index_to_cell_.size(); ++i) {
        const Cell& cell = db.cell(index_to_cell_[i]);
        const double cx = centers_x[i];
        const double cy = centers_y[i];
        // 从 center 转回 lower-left
        state->set_pos(index_to_cell_[i], cx - cell.width * 0.5, cy - cell.height * 0.5);
    }
}

void QuadraticAnalyticalPlacement::run(const PlacementDB& db, PlacementState* state) {
    stats_ = QuadraticStats{};
    if (state == nullptr || db.num_cells() == 0) {
        return;
    }

    const auto& core = db.design.core;
    const double core_width = core.ux - core.lx;
    const double core_height = core.uy - core.ly;
    if (core_width <= 0.0 || core_height <= 0.0) {
        return;
    }

    state->resize(db.num_cells());

    // resize clears positions and orientations; restore the input before solving.
    // Fixed terminals must retain their original coordinates and orientations.
    for (const auto& cell : db.design.cells) {
        state->set_pos(cell.id, cell.init_x, cell.init_y);
        state->set_orient(cell.id, cell.orient);
    }

    // 建立 movable -> index 映射
    build_movable_index_map(db);
    const std::size_t movable_count = index_to_cell_.size();
    stats_.movable_cells = movable_count;
    if (movable_count == 0) return;

    // 初始化：把所有 movable cell 的中心均匀铺在 core 上，
    // 避免全部挤在 core center 导致互距离为 0、权重 w=1/distance 爆炸。
    std::vector<double> centers_x(movable_count), centers_y(movable_count);
    if (movable_count == 1) {
        centers_x[0] = (core.lx + core.ux) * 0.5;
        centers_y[0] = (core.ly + core.uy) * 0.5;
    } else {
        // 近似方形网格，sqrt(n) 行
        const std::size_t cols =
            static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(movable_count))));
        const std::size_t rows = (movable_count + cols - 1) / cols;
        const double margin_x = core_width / static_cast<double>(cols + 1);
        const double margin_y = core_height / static_cast<double>(rows + 1);
        for (std::size_t i = 0; i < movable_count; ++i) {
            const std::size_t r = i / cols;
            const std::size_t c = i % cols;
            centers_x[i] = core.lx + margin_x * static_cast<double>(c + 1);
            centers_y[i] = core.ly + margin_y * static_cast<double>(r + 1);
        }
    }

    // Append solver-only star centers. Each starts at the average current pin
    // location; these entries are deliberately outside index_to_cell_.
    const std::size_t total_unknowns = movable_count +
        static_cast<std::size_t>(std::count_if(
            net_to_star_index_.begin(), net_to_star_index_.end(),
            [](std::size_t index) { return index != std::numeric_limits<std::size_t>::max(); }));
    centers_x.resize(total_unknowns, 0.0);
    centers_y.resize(total_unknowns, 0.0);
    for (const auto& net : db.design.nets) {
        const std::size_t star_index = net_to_star_index_[net.id];
        if (star_index == std::numeric_limits<std::size_t>::max()) {
            continue;
        }
        double sum_x = 0.0;
        double sum_y = 0.0;
        for (PinId pin_id : net.pin_ids) {
            const Pin& pin = db.pin(pin_id);
            const Cell& cell = db.cell(pin.cell_id);
            const bool movable = cell.is_movable();
            const double center_x = movable ? centers_x[cell_to_index_[cell.id]]
                                            : state->x(cell.id) + cell.width * 0.5;
            const double center_y = movable ? centers_y[cell_to_index_[cell.id]]
                                            : state->y(cell.id) + cell.height * 0.5;
            sum_x += pin_abs_x(cell, pin, center_x);
            sum_y += pin_abs_y(cell, pin, center_y);
        }
        centers_x[star_index] = sum_x / static_cast<double>(net.degree());
        centers_y[star_index] = sum_y / static_cast<double>(net.degree());
    }

    // 初始 HPWL
    write_centers_to_state(db, centers_x, centers_y, state);
    stats_.initial_hpwl = HpwlEvaluator{}.compute_total_hpwl(db, *state);

    CSRMatrix A_x, A_y;
    std::vector<double> b_x, b_y;
    std::vector<double> guess_x, guess_y;

    double last_x_error = std::numeric_limits<double>::max();
    double last_y_error = std::numeric_limits<double>::max();

    for (std::size_t iter = 0; iter < config_.max_outer_iterations; ++iter) {
        // 构建 A, b, 初始 guess 来自上一轮结果
        build_system(db, *state, centers_x, centers_y,
                     &A_x, &A_y, &b_x, &b_y,
                     &guess_x, &guess_y);

        // 求解 A_x * x = -b_x / 2
        // 因为能量函数是 E = x^T*A*x + b^T*x，求导得到 2*A*x + b = 0
        // 所以我们真正解的是 A*x = -b/2
        for (auto& val : b_x) val = -val * 0.5;
        for (auto& val : b_y) val = -val * 0.5;

        // Jacobi 预条件（行缩放）：改善 BiCGSTAB 收敛速度
        ApplyJacobiPreconditioner(&A_x, &b_x);
        ApplyJacobiPreconditioner(&A_y, &b_y);

        std::size_t iters_x = 0, iters_y = 0;
        last_x_error = BiCGSTAB(A_x, b_x, &guess_x,
                                 config_.target_error,
                                 config_.solver_max_iterations,
                                 &iters_x);
        last_y_error = BiCGSTAB(A_y, b_y, &guess_y,
                                 config_.target_error,
                                 config_.solver_max_iterations,
                                 &iters_y);

        stats_.total_solver_iterations += iters_x + iters_y;
        stats_.outer_iterations = iter + 1;

        // 更新 centers，并裁剪到 core
        centers_x = guess_x;
        centers_y = guess_y;
        clip_to_core(db, &centers_x, &centers_y);

        // 写回 state，供下一次迭代的 fixed-cell pin 绝对坐标参考
        write_centers_to_state(db, centers_x, centers_y, state);

        // 收敛检查：残差足够小且至少迭代过几次
        if (last_x_error < config_.target_error &&
            last_y_error < config_.target_error &&
            iter > 3) {
            break;
        }
    }

    stats_.final_x_error = last_x_error;
    stats_.final_y_error = last_y_error;
    stats_.final_hpwl = HpwlEvaluator{}.compute_total_hpwl(db, *state);
}

}  // namespace plc

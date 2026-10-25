#include "algorithms/global_placement.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace plc {
namespace {

constexpr double kPi = 3.14159265358979323846;
using Complex = std::complex<double>;

// Iterative radix-2 Cooley-Tukey FFT. The input size must be a power of two.
void fft(std::vector<Complex>& a, bool invert) {
    const std::size_t n = a.size();
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; (j & bit) != 0; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double angle = 2.0 * kPi / static_cast<double>(len) * (invert ? 1.0 : -1.0);
        const Complex wlen(std::cos(angle), std::sin(angle));
        for (std::size_t i = 0; i < n; i += len) {
            Complex w(1.0, 0.0);
            for (std::size_t j = 0; j < len / 2; ++j) {
                const Complex u = a[i + j];
                const Complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        const double inv_n = 1.0 / static_cast<double>(n);
        for (auto& value : a) {
            value *= inv_n;
        }
    }
}

// Raw DCT-II: out[k] = sum_n in[n] * cos(pi*k*(2n+1)/(2N)).
void dct2_1d(const std::vector<double>& in, std::vector<double>& out,
             std::vector<Complex>& buf) {
    const std::size_t n = in.size();
    const std::size_t m = 2 * n;
    buf.assign(m, Complex(0.0, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        buf[i] = in[i];
        buf[m - 1 - i] = in[i];
    }
    fft(buf, false);
    for (std::size_t k = 0; k < n; ++k) {
        const double angle = -kPi * static_cast<double>(k) / static_cast<double>(m);
        out[k] = 0.5 * std::real(buf[k] * Complex(std::cos(angle), std::sin(angle)));
    }
}

// Raw DCT-III: out[n] = sum_k in[k] * cos(pi*k*(2n+1)/(2N)).
void dct3_1d(const std::vector<double>& in, std::vector<double>& out,
             std::vector<Complex>& buf) {
    const std::size_t n = in.size();
    const std::size_t m = 2 * n;
    buf.assign(m, Complex(0.0, 0.0));
    const double first = in[0];
    buf[0] = Complex(first, 0.0);
    for (std::size_t k = 1; k < n; ++k) {
        const double angle = kPi * static_cast<double>(k) / static_cast<double>(m);
        const Complex value = in[k] * Complex(std::cos(angle), std::sin(angle));
        buf[k] = value;
        buf[m - k] = std::conj(value);
    }
    fft(buf, true);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = static_cast<double>(n) * std::real(buf[i]) + 0.5 * first;
    }
}

std::vector<double> dct2_2d(const std::vector<double>& src, std::size_t m) {
    std::vector<double> tmp(m * m, 0.0);
    std::vector<double> dst(m * m, 0.0);
    std::vector<double> in(m, 0.0);
    std::vector<double> out(m, 0.0);
    std::vector<Complex> buf;

    for (std::size_t y = 0; y < m; ++y) {
        for (std::size_t x = 0; x < m; ++x) {
            in[x] = src[y * m + x];
        }
        dct2_1d(in, out, buf);
        for (std::size_t x = 0; x < m; ++x) {
            tmp[y * m + x] = out[x];
        }
    }
    for (std::size_t x = 0; x < m; ++x) {
        for (std::size_t y = 0; y < m; ++y) {
            in[y] = tmp[y * m + x];
        }
        dct2_1d(in, out, buf);
        for (std::size_t y = 0; y < m; ++y) {
            dst[y * m + x] = out[y];
        }
    }
    return dst;
}

std::vector<double> dct3_2d(const std::vector<double>& src, std::size_t m) {
    std::vector<double> tmp(m * m, 0.0);
    std::vector<double> dst(m * m, 0.0);
    std::vector<double> in(m, 0.0);
    std::vector<double> out(m, 0.0);
    std::vector<Complex> buf;

    for (std::size_t y = 0; y < m; ++y) {
        for (std::size_t x = 0; x < m; ++x) {
            in[x] = src[y * m + x];
        }
        dct3_1d(in, out, buf);
        for (std::size_t x = 0; x < m; ++x) {
            tmp[y * m + x] = out[x];
        }
    }
    for (std::size_t x = 0; x < m; ++x) {
        for (std::size_t y = 0; y < m; ++y) {
            in[y] = tmp[y * m + x];
        }
        dct3_1d(in, out, buf);
        for (std::size_t y = 0; y < m; ++y) {
            dst[y * m + x] = out[y];
        }
    }
    return dst;
}

struct Node {
    bool is_cell = false;
    bool filler = false;
    bool fixed = false;
    bool macro = false;
    CellId cell_id = 0;
    double width = 0.0;
    double height = 0.0;
    double area = 0.0;
    double charge = 0.0;
    double x = 0.0;
    double y = 0.0;
};

struct BinGrid {
    std::size_t m = 0;
    double lx = 0.0;
    double ly = 0.0;
    double ux = 0.0;
    double uy = 0.0;
    double bin_w = 0.0;
    double bin_h = 0.0;
    double bin_area = 0.0;
};

std::size_t clamp_bin(double value, double origin, double bin_size, std::size_t m) {
    const double raw = std::floor((value - origin) / bin_size);
    if (raw <= 0.0) {
        return 0;
    }
    if (raw >= static_cast<double>(m - 1)) {
        return m - 1;
    }
    return static_cast<std::size_t>(raw);
}

double clamp_coord(double value, double low, double high) {
    if (high <= low) {
        return low;
    }
    return std::min(std::max(value, low), high);
}

void clamp_node(Node* node, const CoreRegion& core) {
    node->x = clamp_coord(node->x, core.lx, core.ux - node->width);
    node->y = clamp_coord(node->y, core.ly, core.uy - node->height);
}

double clipped_area(double x, double y, double w, double h, const CoreRegion& core) {
    const double x0 = std::max(x, core.lx);
    const double y0 = std::max(y, core.ly);
    const double x1 = std::min(x + w, core.ux);
    const double y1 = std::min(y + h, core.uy);
    if (x1 <= x0 || y1 <= y0) {
        return 0.0;
    }
    return (x1 - x0) * (y1 - y0);
}

// The global placer stores lower-left locations. Pin offsets are center-relative,
// matching HpwlEvaluator's absolute-pin coordinate convention.
double absolute_pin_x(const PlacementDB& db, const std::vector<double>& positions,
                      const Pin& pin) {
    const Cell& cell = db.cell(pin.cell_id);
    return positions[pin.cell_id] + cell.width * 0.5 + pin.offset_x;
}

double absolute_pin_y(const PlacementDB& db, const std::vector<double>& positions,
                      const Pin& pin) {
    const Cell& cell = db.cell(pin.cell_id);
    return positions[pin.cell_id] + cell.height * 0.5 + pin.offset_y;
}

// Accumulate the density contribution of a node into the bin grid. A node that
// is smaller than a bin is smeared uniformly over a bin-sized window centred on
// the node (the LengthScale idea of the task description), so the total mass the
// node injects stays exactly equal to its area * scale. This keeps the density
// map used by the Poisson solve and the density overflow metric consistent.
void add_to_bins(const BinGrid& grid, double x, double y, double w, double h,
                 double scale, std::vector<double>& density) {
    if (w <= 0.0 || h <= 0.0) {
        return;
    }
    const double eff_w = std::max(w, grid.bin_w);
    const double eff_h = std::max(h, grid.bin_h);
    const double eff_x = x + 0.5 * w - 0.5 * eff_w;
    const double eff_y = y + 0.5 * h - 0.5 * eff_h;

    const double x0 = std::max(eff_x, grid.lx);
    const double y0 = std::max(eff_y, grid.ly);
    const double x1 = std::min(eff_x + eff_w, grid.ux);
    const double y1 = std::min(eff_y + eff_h, grid.uy);
    if (x1 <= x0 || y1 <= y0) {
        return;
    }

    // Uniform value over the (clipped) effective rectangle: mass = area.
    const double mass = w * h * scale;
    const double density_per_area = mass / ((x1 - x0) * (y1 - y0));

    const std::size_t bx0 = clamp_bin(x0, grid.lx, grid.bin_w, grid.m);
    const std::size_t bx1 = clamp_bin(x1, grid.lx, grid.bin_w, grid.m);
    const std::size_t by0 = clamp_bin(y0, grid.ly, grid.bin_h, grid.m);
    const std::size_t by1 = clamp_bin(y1, grid.ly, grid.bin_h, grid.m);

    for (std::size_t by = by0; by <= by1; ++by) {
        const double bin_y0 = grid.ly + static_cast<double>(by) * grid.bin_h;
        const double bin_y1 = bin_y0 + grid.bin_h;
        const double overlap_y = std::min(y1, bin_y1) - std::max(y0, bin_y0);
        if (overlap_y <= 0.0) {
            continue;
        }
        for (std::size_t bx = bx0; bx <= bx1; ++bx) {
            const double bin_x0 = grid.lx + static_cast<double>(bx) * grid.bin_w;
            const double bin_x1 = bin_x0 + grid.bin_w;
            const double overlap_x = std::min(x1, bin_x1) - std::max(x0, bin_x0);
            if (overlap_x <= 0.0) {
                continue;
            }
            const std::size_t index = by * grid.m + bx;
            density[index] += overlap_x * overlap_y * density_per_area;
        }
    }
}

double compute_hpwl(const PlacementDB& db, const std::vector<double>& px,
                    const std::vector<double>& py) {
    double total = 0.0;
    for (NetId net_id = 0; net_id < db.num_nets(); ++net_id) {
        const Net& net = db.net(net_id);
        if (net.pin_ids.size() <= 1) {
            continue;
        }
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double max_y = std::numeric_limits<double>::lowest();
        for (PinId pin_id : net.pin_ids) {
            const Pin& pin = db.pin(pin_id);
            const double x = absolute_pin_x(db, px, pin);
            const double y = absolute_pin_y(db, py, pin);
            min_x = std::min(min_x, x);
            min_y = std::min(min_y, y);
            max_x = std::max(max_x, x);
            max_y = std::max(max_y, y);
        }
        total += ((max_x - min_x) + (max_y - min_y)) * net.weight;
    }
    return total;
}

}  // namespace

GlobalPlacer::GlobalPlacer(GlobalPlacementConfig config) : config_(config) {
    config_.target_density = std::min(std::max(config_.target_density, 0.1), 1.0);
    config_.target_overflow = std::min(std::max(config_.target_overflow, 0.0), 1.0);
    config_.max_bins = std::max<std::size_t>(8, config_.max_bins);
    config_.lambda_reference = std::max(config_.lambda_reference, 1.0);
}

void GlobalPlacer::run(const PlacementDB& db, PlacementState* state) {
    stats_ = GlobalPlacementStats{};
    if (state == nullptr) {
        return;
    }
    // Keep the incoming placement (e.g. the random initial placement); only
    // allocate the state when it has not been initialized yet.
    if (state->size() < db.num_cells()) {
        state->resize(db.num_cells());
    }

    const CoreRegion& core = db.design.core;
    const double core_w = core.ux - core.lx;
    const double core_h = core.uy - core.ly;
    if (core_w <= 0.0 || core_h <= 0.0) {
        return;
    }

    // 1. NodesAndFillers: collect movable cells, fixed terminals and fillers.
    std::vector<Node> nodes;
    nodes.reserve(db.num_cells());
    std::vector<double> cell_areas;
    double movable_area = 0.0;
    double fixed_area_in_core = 0.0;

    for (const Cell& cell : db.design.cells) {
        Node node;
        node.is_cell = true;
        node.cell_id = cell.id;
        node.width = cell.width;
        node.height = cell.height;
        node.area = cell.area();
        node.fixed = !cell.is_movable();
        node.macro = (cell.type == CellType::kMacro);
        node.x = state->x(cell.id);
        node.y = state->y(cell.id);
        if (node.fixed) {
            fixed_area_in_core += clipped_area(node.x, node.y, node.width, node.height, core);
        } else {
            movable_area += node.area;
            cell_areas.push_back(node.area);
            clamp_node(&node, core);
        }
        nodes.push_back(node);
    }

    if (movable_area <= 0.0) {
        return;
    }

    // FillerInit: total filler area closes the gap to the target density. The
    // single filler size is the average area of the middle 80% of the cells.
    const double target_total = config_.target_density * core_w * core_h;
    const double needed_filler_area =
        std::max(0.0, target_total - (movable_area + fixed_area_in_core));
    double filler_side = 0.0;
    std::size_t filler_count = 0;
    double filler_total_area = 0.0;

    if (needed_filler_area > 0.0) {
        std::sort(cell_areas.begin(), cell_areas.end());
        const std::size_t count = cell_areas.size();
        std::size_t begin = count / 10;
        std::size_t end = count - count / 10;
        if (end <= begin) {
            begin = 0;
            end = count;
        }
        double sum = 0.0;
        for (std::size_t i = begin; i < end; ++i) {
            sum += cell_areas[i];
        }
        const double mean_area = sum / static_cast<double>(std::max<std::size_t>(1, end - begin));
        const double max_side = std::min(core_w, core_h);
        filler_side = std::sqrt(std::max(mean_area, 1e-12));
        filler_side = std::min(filler_side, max_side);
        filler_side = std::max(filler_side, max_side * 1e-4);
        filler_count = static_cast<std::size_t>(
            std::ceil(needed_filler_area / (filler_side * filler_side)));

        const std::size_t cap = std::max<std::size_t>(2000, nodes.size() * 3);
        if (filler_count > cap) {
            filler_count = cap;
            filler_side = std::sqrt(needed_filler_area / static_cast<double>(filler_count));
        }
        filler_count = std::max<std::size_t>(1, filler_count);
        filler_total_area =
            static_cast<double>(filler_count) * filler_side * filler_side;

        std::mt19937 rng(config_.seed);
        std::uniform_real_distribution<double> x_dist(
            core.lx, std::max(core.lx, core.ux - filler_side));
        std::uniform_real_distribution<double> y_dist(
            core.ly, std::max(core.ly, core.uy - filler_side));
        nodes.reserve(nodes.size() + filler_count);
        for (std::size_t i = 0; i < filler_count; ++i) {
            Node filler;
            filler.filler = true;
            filler.width = filler_side;
            filler.height = filler_side;
            filler.area = filler_side * filler_side;
            filler.x = x_dist(rng);
            filler.y = y_dist(rng);
            nodes.push_back(filler);
        }
    }

    const std::size_t node_count = nodes.size();
    const double total_charge = movable_area + filler_total_area;
    for (Node& node : nodes) {
        node.charge = node.area / std::max(total_charge, 1e-12);
    }

    // BinInit: square-ish m x m grid whose dimension follows sqrt(n').
    std::size_t m = 8;
    while (m * m < node_count && m < config_.max_bins) {
        m <<= 1;
    }
    const std::size_t bin_count = m * m;
    BinGrid grid;
    grid.m = m;
    grid.lx = core.lx;
    grid.ly = core.ly;
    grid.ux = core.ux;
    grid.uy = core.uy;
    grid.bin_w = core_w / static_cast<double>(m);
    grid.bin_h = core_h / static_cast<double>(m);
    grid.bin_area = grid.bin_w * grid.bin_h;

    stats_.bins_x = m;
    stats_.bins_y = m;
    stats_.fillers = filler_count;

    // Fixed terminal density and the non-placeable (dark) area density.
    std::vector<double> terminal_density(bin_count, 0.0);
    for (const Node& node : nodes) {
        if (!node.fixed) {
            continue;
        }
        add_to_bins(grid, node.x, node.y, node.width, node.height,
                    config_.target_density, terminal_density);
    }

    std::vector<double> dark_density(bin_count, 0.0);
    if (!db.design.rows.empty()) {
        std::vector<std::pair<double, double>> intervals;
        intervals.reserve(db.design.rows.size());
        for (const Row& row : db.design.rows) {
            const double low = std::max(row.origin_y, core.ly);
            const double high = std::min(row.origin_y + row.height, core.uy);
            if (high > low) {
                intervals.emplace_back(low, high);
            }
        }
        std::sort(intervals.begin(), intervals.end());
        std::vector<std::pair<double, double>> merged;
        for (const auto& interval : intervals) {
            if (merged.empty() || interval.first > merged.back().second) {
                merged.push_back(interval);
            } else {
                merged.back().second = std::max(merged.back().second, interval.second);
            }
        }
        for (std::size_t by = 0; by < m; ++by) {
            const double bin_y0 = grid.ly + static_cast<double>(by) * grid.bin_h;
            const double bin_y1 = bin_y0 + grid.bin_h;
            double covered = 0.0;
            for (const auto& interval : merged) {
                covered += std::max(0.0, std::min(bin_y1, interval.second) -
                                             std::max(bin_y0, interval.first));
            }
            const double free_height = std::max(0.0, grid.bin_h - covered);
            for (std::size_t bx = 0; bx < m; ++bx) {
                const std::size_t index = by * m + bx;
                dark_density[index] =
                    free_height * grid.bin_w * config_.target_density;
            }
        }
    }

    std::vector<double> pos_x(node_count, 0.0);
    std::vector<double> pos_y(node_count, 0.0);
    for (std::size_t i = 0; i < node_count; ++i) {
        pos_x[i] = nodes[i].x;
        pos_y[i] = nodes[i].y;
    }

    // Reused buffers for the density map, DCT field and gradients.
    std::vector<double> density_map(bin_count, 0.0);
    std::vector<double> density_grad_x(node_count, 0.0);
    std::vector<double> density_grad_y(node_count, 0.0);
    std::vector<double> wl_grad_x(node_count, 0.0);
    std::vector<double> wl_grad_y(node_count, 0.0);
    std::vector<double> grad_x(node_count, 0.0);
    std::vector<double> grad_y(node_count, 0.0);
    std::vector<double> prev_grad_x(node_count, 0.0);
    std::vector<double> prev_grad_y(node_count, 0.0);
    std::vector<double> prev_pos_x(node_count, 0.0);
    std::vector<double> prev_pos_y(node_count, 0.0);
    std::vector<double> major_x(node_count, 0.0);
    std::vector<double> major_y(node_count, 0.0);
    std::vector<double> next_x(node_count, 0.0);
    std::vector<double> next_y(node_count, 0.0);
    std::vector<double> dct_scratch(bin_count, 0.0);
    std::vector<double> potential(bin_count, 0.0);
    std::vector<double> field_x(bin_count, 0.0);
    std::vector<double> field_y(bin_count, 0.0);

    // Build the density map, transform it with the DCT, evaluate the potential
    // and its gradient on the bin centers, then compute the density gradient of
    // every node. Also returns the density overflow ratio tau.
    auto evaluate_density = [&](const std::vector<double>& px,
                                const std::vector<double>& py, double* tau_out) {
        std::fill(density_map.begin(), density_map.end(), 0.0);

        for (std::size_t i = 0; i < node_count; ++i) {
            const Node& node = nodes[i];
            if (node.fixed) {
                continue;
            }
            const double scale =
                (node.filler || node.macro) ? config_.target_density : 1.0;
            add_to_bins(grid, px[i], py[i], node.width, node.height, scale,
                        density_map);
        }
        for (std::size_t b = 0; b < bin_count; ++b) {
            density_map[b] += terminal_density[b] + dark_density[b];
        }

        // Density overflow: total area above the target density, normalised by
        // the total movable (cell + filler) area.
        double overflow_area = 0.0;
        for (std::size_t b = 0; b < bin_count; ++b) {
            overflow_area += std::max(
                0.0, density_map[b] - config_.target_density * grid.bin_area);
        }
        const double normalization =
            std::max(movable_area + filler_total_area, 1e-12);
        *tau_out = overflow_area / normalization;

        for (std::size_t b = 0; b < bin_count; ++b) {
            dct_scratch[b] = density_map[b] / grid.bin_area;
        }

        // a_{j,k}: DCT coefficients of the density distribution.
        const std::vector<double> coefficients = dct2_2d(dct_scratch, m);
        const double inv_m2 = 1.0 / static_cast<double>(m * m);
        std::fill(dct_scratch.begin(), dct_scratch.end(), 0.0);
        for (std::size_t k = 0; k < m; ++k) {
            const double wk = kPi * static_cast<double>(k) / static_cast<double>(m);
            for (std::size_t j = 0; j < m; ++j) {
                if (j == 0 && k == 0) {
                    continue;
                }
                const double wj = kPi * static_cast<double>(j) / static_cast<double>(m);
                dct_scratch[k * m + j] =
                    coefficients[k * m + j] * inv_m2 / (wj * wj + wk * wk);
            }
        }

        // Potential psi at the bin centers (inverse DCT).
        potential = dct3_2d(dct_scratch, m);

        // Electric field = -grad(psi), evaluated with central differences in
        // normalized coordinates. Only the potential gradient is needed for
        // the density force.
        const double step = 1.0 / static_cast<double>(m);
        for (std::size_t by = 0; by < m; ++by) {
            for (std::size_t bx = 0; bx < m; ++bx) {
                const std::size_t index = by * m + bx;
                const double left = (bx > 0) ? potential[index - 1] : potential[index];
                const double right =
                    (bx + 1 < m) ? potential[index + 1] : potential[index];
                const double denom = (bx > 0 && bx + 1 < m) ? 2.0 * step : step;
                field_x[index] = (right - left) / denom;

                const double down = (by > 0) ? potential[index - m] : potential[index];
                const double up =
                    (by + 1 < m) ? potential[index + m] : potential[index];
                const double denom_y = (by > 0 && by + 1 < m) ? 2.0 * step : step;
                field_y[index] = (up - down) / denom_y;
            }
        }

        const double inv_core_w = 1.0 / core_w;
        const double inv_core_h = 1.0 / core_h;
        for (std::size_t i = 0; i < node_count; ++i) {
            const Node& node = nodes[i];
            // Average the bin electric field over the bins covered by the node,
            // weighted by the overlap area. Sampling a single bin would give a
            // zero force whenever a node sits on a local maximum of the field.
            const double x0 = std::max(px[i], grid.lx);
            const double y0 = std::max(py[i], grid.ly);
            const double x1 = std::min(px[i] + node.width, grid.ux);
            const double y1 = std::min(py[i] + node.height, grid.uy);
            double gx = 0.0;
            double gy = 0.0;
            double weight = 0.0;
            if (x1 > x0 && y1 > y0) {
                const std::size_t bx0 = clamp_bin(x0, grid.lx, grid.bin_w, m);
                const std::size_t bx1 = clamp_bin(x1, grid.lx, grid.bin_w, m);
                const std::size_t by0 = clamp_bin(y0, grid.ly, grid.bin_h, m);
                const std::size_t by1 = clamp_bin(y1, grid.ly, grid.bin_h, m);
                for (std::size_t by = by0; by <= by1; ++by) {
                    const double bin_y0 = grid.ly + static_cast<double>(by) * grid.bin_h;
                    const double overlap_y = std::min(y1, bin_y0 + grid.bin_h) -
                                             std::max(y0, bin_y0);
                    if (overlap_y <= 0.0) {
                        continue;
                    }
                    for (std::size_t bx = bx0; bx <= bx1; ++bx) {
                        const double bin_x0 = grid.lx + static_cast<double>(bx) * grid.bin_w;
                        const double overlap_x = std::min(x1, bin_x0 + grid.bin_w) -
                                                 std::max(x0, bin_x0);
                        if (overlap_x <= 0.0) {
                            continue;
                        }
                        const double overlap = overlap_x * overlap_y;
                        const std::size_t index = by * m + bx;
                        gx += field_x[index] * overlap;
                        gy += field_y[index] * overlap;
                        weight += overlap;
                    }
                }
            }
            if (weight > 0.0) {
                gx /= weight;
                gy /= weight;
            }
            density_grad_x[i] = node.charge * gx * inv_core_w;
            density_grad_y[i] = node.charge * gy * inv_core_h;
        }
    };

    double tau = 0.0;
    double hpwl = compute_hpwl(db, pos_x, pos_y);
    double previous_hpwl = hpwl;
    double lambda = 0.0;
    double previous_alpha = 0.0;
    double nesterov_parameter = 1.0;
    std::size_t performed_iterations = 0;

    std::vector<double> t_plus;
    std::vector<double> t_minus;

    stats_.initial_hpwl = hpwl;

    for (std::size_t iteration = 0; iteration < config_.max_iterations; ++iteration) {
        evaluate_density(pos_x, pos_y, &tau);
        hpwl = compute_hpwl(db, pos_x, pos_y);
        if (iteration == 0) {
            stats_.initial_overflow = tau;
            stats_.initial_hpwl = hpwl;
            previous_hpwl = hpwl;
        }
        if (iteration > 0 && tau < config_.target_overflow) {
            performed_iterations = iteration;
            break;
        }

        // Wirelength gradient with the WA smoothing model. gamma follows
        // gamma(tau) = 8 * w_b * 10^((20/9)(tau - 0.1) - 1).
        std::fill(wl_grad_x.begin(), wl_grad_x.end(), 0.0);
        std::fill(wl_grad_y.begin(), wl_grad_y.end(), 0.0);
        const double tau_clamped = std::min(std::max(tau, 0.1), 1.0);
        double gamma = 8.0 * grid.bin_w *
                       std::pow(10.0, (20.0 / 9.0) * (tau_clamped - 0.1) - 1.0);
        gamma = std::max(gamma, 1e-9);

        for (NetId net_id = 0; net_id < db.num_nets(); ++net_id) {
            const Net& net = db.net(net_id);
            const std::size_t degree = net.pin_ids.size();
            if (degree <= 1) {
                continue;
            }
            t_plus.resize(degree);
            t_minus.resize(degree);

            // x direction
            double min_value = std::numeric_limits<double>::max();
            double max_value = std::numeric_limits<double>::lowest();
            for (PinId pin_id : net.pin_ids) {
                const Pin& pin = db.pin(pin_id);
                const double value = absolute_pin_x(db, pos_x, pin);
                min_value = std::min(min_value, value);
                max_value = std::max(max_value, value);
            }
            if (max_value > min_value) {
                double b_plus = 0.0;
                double b_minus = 0.0;
                double c_plus = 0.0;
                double c_minus = 0.0;
                std::size_t index = 0;
                for (PinId pin_id : net.pin_ids) {
                    const Pin& pin = db.pin(pin_id);
                    const double value = absolute_pin_x(db, pos_x, pin);
                    t_plus[index] = std::exp((value - max_value) / gamma);
                    t_minus[index] = std::exp((min_value - value) / gamma);
                    b_plus += t_plus[index];
                    b_minus += t_minus[index];
                    c_plus += value * t_plus[index];
                    c_minus += value * t_minus[index];
                    ++index;
                }
                const double center_plus = c_plus / std::max(b_plus, 1e-30);
                const double center_minus = c_minus / std::max(b_minus, 1e-30);
                index = 0;
                for (PinId pin_id : net.pin_ids) {
                    const Pin& pin = db.pin(pin_id);
                    const double value = absolute_pin_x(db, pos_x, pin);
                    const double positive =
                        t_plus[index] / std::max(b_plus, 1e-30) *
                        (1.0 + (value - center_plus) / gamma);
                    const double negative =
                        t_minus[index] / std::max(b_minus, 1e-30) *
                        (1.0 - (value - center_minus) / gamma);
                    wl_grad_x[pin.cell_id] += net.weight * (positive - negative);
                    ++index;
                }
            }

            // y direction
            min_value = std::numeric_limits<double>::max();
            max_value = std::numeric_limits<double>::lowest();
            for (PinId pin_id : net.pin_ids) {
                const Pin& pin = db.pin(pin_id);
                const double value = absolute_pin_y(db, pos_y, pin);
                min_value = std::min(min_value, value);
                max_value = std::max(max_value, value);
            }
            if (max_value > min_value) {
                double b_plus = 0.0;
                double b_minus = 0.0;
                double c_plus = 0.0;
                double c_minus = 0.0;
                std::size_t index = 0;
                for (PinId pin_id : net.pin_ids) {
                    const Pin& pin = db.pin(pin_id);
                    const double value = absolute_pin_y(db, pos_y, pin);
                    t_plus[index] = std::exp((value - max_value) / gamma);
                    t_minus[index] = std::exp((min_value - value) / gamma);
                    b_plus += t_plus[index];
                    b_minus += t_minus[index];
                    c_plus += value * t_plus[index];
                    c_minus += value * t_minus[index];
                    ++index;
                }
                const double center_plus = c_plus / std::max(b_plus, 1e-30);
                const double center_minus = c_minus / std::max(b_minus, 1e-30);
                index = 0;
                for (PinId pin_id : net.pin_ids) {
                    const Pin& pin = db.pin(pin_id);
                    const double value = absolute_pin_y(db, pos_y, pin);
                    const double positive =
                        t_plus[index] / std::max(b_plus, 1e-30) *
                        (1.0 + (value - center_plus) / gamma);
                    const double negative =
                        t_minus[index] / std::max(b_minus, 1e-30) *
                        (1.0 - (value - center_minus) / gamma);
                    wl_grad_y[pin.cell_id] += net.weight * (positive - negative);
                    ++index;
                }
            }
        }

        // Penalty factor lambda: initialize by balancing both gradient norms,
        // then adapt it from the HPWL variation, lambda_k = mu * lambda_{k-1}.
        if (iteration == 0) {
            double numerator = 0.0;
            double denominator = 0.0;
            for (std::size_t i = 0; i < node_count; ++i) {
                numerator += std::fabs(wl_grad_x[i]) + std::fabs(wl_grad_y[i]);
                denominator += std::fabs(density_grad_x[i]) +
                               std::fabs(density_grad_y[i]);
            }
            lambda = (denominator > 1e-12) ? numerator / denominator : 1.0;
            lambda = std::max(lambda, 1e-12);
        } else {
            const double delta_hpwl = hpwl - previous_hpwl;
            double mu = std::pow(
                1.1, 1.0 - delta_hpwl / config_.lambda_reference);
            mu = std::min(std::max(mu, 0.75), 1.1);
            lambda = std::max(lambda * mu, 1e-12);
        }

        // Total gradient: grad(HPWL) + lambda * grad(density).
        for (std::size_t i = 0; i < node_count; ++i) {
            grad_x[i] = wl_grad_x[i] + lambda * density_grad_x[i];
            grad_y[i] = wl_grad_y[i] + lambda * density_grad_y[i];
        }

        double max_grad = 0.0;
        for (std::size_t i = 0; i < node_count; ++i) {
            if (nodes[i].fixed) {
                continue;
            }
            max_grad = std::max(max_grad, std::fabs(grad_x[i]));
            max_grad = std::max(max_grad, std::fabs(grad_y[i]));
        }

        double alpha = 0.0;
        if (iteration == 0) {
            const double step_limit =
                0.05 * std::max(grid.bin_w, grid.bin_h);
            alpha = (max_grad > 0.0) ? step_limit / max_grad : 0.0;
        } else {
            double grad_diff_sq = 0.0;
            double pos_diff_sq = 0.0;
            for (std::size_t i = 0; i < node_count; ++i) {
                if (nodes[i].fixed) {
                    continue;
                }
                const double dx = grad_x[i] - prev_grad_x[i];
                const double dy = grad_y[i] - prev_grad_y[i];
                grad_diff_sq += dx * dx + dy * dy;
                const double px = pos_x[i] - prev_pos_x[i];
                const double py = pos_y[i] - prev_pos_y[i];
                pos_diff_sq += px * px + py * py;
            }
            alpha = previous_alpha;
            if (pos_diff_sq > 1e-24) {
                const double lipschitz = std::sqrt(grad_diff_sq / pos_diff_sq);
                if (lipschitz > 1e-12) {
                    alpha = 1.0 / lipschitz;
                }
            }
            if (max_grad > 0.0) {
                const double max_step = 0.5 * std::max(grid.bin_w, grid.bin_h);
                const double displacement = alpha * max_grad;
                if (displacement > max_step) {
                    alpha *= max_step / displacement;
                }
                alpha = std::min(std::max(alpha, previous_alpha * 0.25),
                                 previous_alpha * 4.0);
            }
        }
        alpha = std::max(alpha, 0.0);
        previous_alpha = alpha;
        stats_.final_lambda = lambda;

        // Nesterov update: u_{k+1} = v_k - alpha * grad(v_k) and
        // v_{k+1} = u_{k+1} + (a_k - 1) / a_{k+1} * (u_{k+1} - u_k).
        const double next_parameter =
            0.5 * (1.0 + std::sqrt(4.0 * nesterov_parameter * nesterov_parameter + 1.0));
        const double momentum =
            (nesterov_parameter - 1.0) / std::max(next_parameter, 1e-12);
        for (std::size_t i = 0; i < node_count; ++i) {
            if (nodes[i].fixed) {
                next_x[i] = pos_x[i];
                next_y[i] = pos_y[i];
                continue;
            }
            const double ux = clamp_coord(pos_x[i] - alpha * grad_x[i],
                                          core.lx, core.ux - nodes[i].width);
            const double uy = clamp_coord(pos_y[i] - alpha * grad_y[i],
                                          core.ly, core.uy - nodes[i].height);
            next_x[i] = clamp_coord(ux + momentum * (ux - major_x[i]),
                                    core.lx, core.ux - nodes[i].width);
            next_y[i] = clamp_coord(uy + momentum * (uy - major_y[i]),
                                    core.ly, core.uy - nodes[i].height);
            major_x[i] = ux;
            major_y[i] = uy;
        }
        nesterov_parameter = next_parameter;

        prev_grad_x.swap(grad_x);
        prev_grad_y.swap(grad_y);
        prev_pos_x.swap(pos_x);
        prev_pos_y.swap(pos_y);
        pos_x.swap(next_x);
        pos_y.swap(next_y);
        previous_hpwl = hpwl;
        performed_iterations = iteration + 1;

        if (config_.verbose && ((iteration + 1) % 10 == 0)) {
            std::cout << "  [GP] iter " << (iteration + 1) << "  HPWL=" << hpwl
                      << "  overflow=" << tau * 100.0 << "%  lambda=" << lambda
                      << '\n';
        }
    }

    stats_.iterations = performed_iterations;

    // Remove fillers from the solution: sync the cells back to the state.
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const Node& node = nodes[i];
        if (node.is_cell && !node.fixed) {
            state->set_pos(node.cell_id, pos_x[i], pos_y[i]);
        }
    }

    evaluate_density(pos_x, pos_y, &tau);
    stats_.final_overflow = tau;
    stats_.final_hpwl = compute_hpwl(db, pos_x, pos_y);
    stats_.final_lambda = lambda;
}

}  // namespace plc

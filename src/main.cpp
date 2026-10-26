#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "algorithms/clustering_placement.h"
#include "algorithms/global_placement.h"
#include "algorithms/local_sa_placer.h"
#include "algorithms/quadratic_analytical_placement.h"
#include "algorithms/random_placement.h"
#include "io/summary_writer.h"
#include "metrics/hpwl.h"
#include "parser/bookshelf_loader.h"
#include "visualization/cimg_visualizer.h"

namespace {

struct Options {
    std::string method = "random", output;
    unsigned seed = 42;
    bool check_only = false, sa = false;
    bool run_gp = false;          // task 5: run global placement after initializer
    std::string plot;
    plc::PlotOptions plot_options;
    std::size_t compare_timeout_ms = 0;  // 0 disables compare preflight timeout.
    plc::ClusteringConfig clustering;
    plc::GlobalPlacementConfig gp_config;
};

struct MethodSummary {
    std::string method;
    double hpwl = std::numeric_limits<double>::quiet_NaN();
    double overflow_percent = std::numeric_limits<double>::quiet_NaN();
    double runtime_ms = std::numeric_limits<double>::quiet_NaN();
    std::string status = "completed";
};

void usage(const char* program) {
    std::cerr << "Usage: " << program << " <benchmark.aux> [options]\n"
              << "  --method random|cluster|quadratic|compare   (default random)\n"
              << "  --seed N  --clusters N  --max-cluster-cells N\n"
              << "  --sa | --no-sa                               (default: no SA)\n"
              << "  --global                                     (run task-5 global placement)\n"
              << "  --iterations N  --target-density v  --target-overflow v\n"
              << "  --bins N  --quiet                            (GP tuning)\n"
              << "  --compare-timeout-ms N                       (per-method preflight budget)\n"
              << "  --check-only  --output result.pl\n"
              << "  --plot image.png|image.bmp  --plot-width N  --plot-margin N\n"
              << "  --plot-rows|--no-plot-rows  --plot-nets\n"
              << "\n"
              << "Task 4: --method random|cluster|quadratic each produce an initial placement.\n"
              << "        --method compare runs all three independently.\n"
              << "        --sa runs local simulated annealing after the chosen initializer.\n"
              << "\n"
              << "Task 5: --global runs the ePlace-style electrostatics-based global placement\n"
              << "        after the selected initializer (random if none specified).\n"
              << "        Tunables: --iterations, --target-density, --target-overflow, --bins, --quiet.\n"
              << "\n"
              << "  --output result.pl exports the placement (after GP when --global is on).\n"
              << "  compare exports become result.<method>.pl; SA is applied per initializer.\n";
}

Options parse(int argc, char** argv) {
    Options o;
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument("Missing value for " + arg);
            return argv[++i];
        };
        auto number = [&]() -> std::size_t {
            auto v = value(); std::size_t end = 0;
            if (v.empty() || v[0] == '-' || v[0] == '+') throw std::invalid_argument("Invalid number for " + arg);
            auto n = std::stoull(v, &end);
            if (end != v.size() || n > std::numeric_limits<std::size_t>::max())
                throw std::invalid_argument("Invalid number for " + arg);
            return static_cast<std::size_t>(n);
        };
        auto real = [&]() -> double {
            auto v = value();
            if (v.empty()) throw std::invalid_argument("Invalid real for " + arg);
            return std::stod(v);
        };
        if (arg == "--method") o.method = value();
        else if (arg == "--output") o.output = value();
        else if (arg == "--plot") o.plot = value();
        else if (arg == "--plot-width") o.plot_options.width = number();
        else if (arg == "--plot-margin") o.plot_options.margin = number();
        else if (arg == "--plot-rows") o.plot_options.draw_rows = true;
        else if (arg == "--no-plot-rows") o.plot_options.draw_rows = false;
        else if (arg == "--plot-nets") o.plot_options.draw_nets = true;
        else if (arg == "--seed") {
            auto n = number();
            if (n > std::numeric_limits<unsigned>::max()) throw std::invalid_argument("Seed out of range");
            o.seed = static_cast<unsigned>(n);
            o.gp_config.seed = o.seed;
        } else if (arg == "--clusters") o.clustering.target_clusters = number();
        else if (arg == "--max-cluster-cells") {
            o.clustering.max_cluster_cells = number();
            if (!o.clustering.max_cluster_cells) throw std::invalid_argument("Cluster cap must be positive");
        } else if (arg == "--sa") o.sa = true;
        else if (arg == "--no-sa") o.sa = false;
        else if (arg == "--check-only") o.check_only = true;
        // ---- task 5 flags ----
        else if (arg == "--global") o.run_gp = true;
        else if (arg == "--local-sa") { o.run_gp = true; o.sa = true; }
        else if (arg == "--quiet") o.gp_config.verbose = false;
        else if (arg == "--iterations") o.gp_config.max_iterations = number();
        else if (arg == "--target-density") o.gp_config.target_density = real();
        else if (arg == "--target-overflow") o.gp_config.target_overflow = real();
        else if (arg == "--bins") o.gp_config.max_bins = number();
        else if (arg == "--compare-timeout-ms") o.compare_timeout_ms = number();
        else throw std::invalid_argument("Unknown option: " + arg);
    }
    if (o.method != "random" && o.method != "cluster" && o.method != "quadratic" && o.method != "compare")
        throw std::invalid_argument("Unknown method: " + o.method);
    return o;
}

std::string metric(double value, int precision) {
    if (!std::isfinite(value)) return "n/a";
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

// Placement algorithms do not expose cancellation points. Starting an
// over-budget quadratic solve in a detached thread would keep consuming CPU
// while later methods run, so compare uses a conservative preflight budget.
bool exceeds_compare_budget(const plc::PlacementDB& db, const Options& options,
                            const std::string& method) {
    if (options.compare_timeout_ms == 0 || method != "quadratic") return false;

    const plc::QuadraticConfig config;
    const double work = static_cast<double>(db.num_pins()) *
                        static_cast<double>(config.max_outer_iterations) *
                        static_cast<double>(config.solver_max_iterations);
    // Calibrated conservatively for the CSR + BiCGSTAB implementation. This
    // is a scheduling estimate, not a replacement for algorithm convergence.
    const double estimated_ms = work / 20000.0;
    return estimated_ms > static_cast<double>(options.compare_timeout_ms);
}

void print_compare_summary(const std::vector<MethodSummary>& summaries) {
    std::cout << "\n---------------COMPARE SUMMARY---------------\n"
              << std::left << std::setw(12) << "method"
              << std::right << std::setw(18) << "HPWL"
              << std::setw(14) << "overflow%"
              << std::setw(15) << "runtime"
              << std::setw(14) << "status" << '\n';
    for (const auto& summary : summaries) {
        std::cout << std::left << std::setw(12) << summary.method
                  << std::right << std::setw(18) << metric(summary.hpwl, 2)
                  << std::setw(14) << metric(summary.overflow_percent, 2)
                  << std::setw(15) << (std::isfinite(summary.runtime_ms)
                      ? metric(summary.runtime_ms, 3) + " ms" : "n/a")
                  << std::setw(14) << summary.status << '\n';
    }
}

void validate(const plc::PlacementDB& db, const plc::PlacementState& state) {
    std::size_t violations = 0;
    const auto& core = db.design.core;
    for (const auto& cell : db.design.cells) {
        const double x = state.x(cell.id), y = state.y(cell.id);
        if (!std::isfinite(x) || !std::isfinite(y)) { ++violations; continue; }
        if (cell.is_fixed() && (x != cell.init_x || y != cell.init_y || state.orient(cell.id) != cell.orient)) ++violations;
        if (cell.is_movable() && (x < core.lx - 1e-8 || y < core.ly - 1e-8 ||
            x + cell.width > core.ux + 1e-8 || y + cell.height > core.uy + 1e-8)) ++violations;
    }
    std::cout << "    Boundary/Fixed Violations: " << violations << '\n';
    if (violations) throw std::runtime_error("Placement validation failed");
}

void export_initial(const plc::PlacementDB& db, const plc::PlacementState& state,
                    const Options& options, const std::string& method) {
    if (options.output.empty()) return;
    std::filesystem::path path(options.output);
    if (options.method == "compare")
        path = path.parent_path() / (path.stem().string() + "." + method + path.extension().string());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open output: " + path.string());
    out << "UCLA pl 1.0\n" << std::setprecision(17);
    const char* orientations[] = {"N", "S", "E", "W", "FN", "FS", "FE", "FW", "N"};
    for (const auto& cell : db.design.cells) {
        out << cell.name << ' ' << state.x(cell.id) << ' ' << state.y(cell.id)
            << " : " << orientations[static_cast<int>(state.orient(cell.id))];
        if (cell.is_fixed()) out << (cell.place_status == plc::PlaceStatus::kFixedNI ? " /FIXED_NI" : " /FIXED");
        out << '\n';
    }
    out.close();
    if (!out) throw std::runtime_error("Failed to write output: " + path.string());
    std::cout << "    Placement Written To: " << path.string() << '\n';
}

std::string plot_path(const Options& options, const std::string& method,
                      const std::string& stage, bool is_final) {
    if (options.plot.empty()) return {};
    std::filesystem::path path(options.plot);
    if (!path.has_extension()) {
        path += plc::CImgVisualizer::default_extension();
    }

    std::string stem = path.stem().string();
    if (options.method == "compare") {
        stem += "_" + method;
    }
    if (!is_final) {
        stem += "." + stage;
    }
    return (path.parent_path() / (stem + path.extension().string())).string();
}

void plot_stage(const plc::PlacementDB& db, const plc::PlacementState& state,
                const Options& options, const std::string& method,
                const std::string& stage, bool is_final) {
    const std::string path = plot_path(options, method, stage, is_final);
    if (path.empty()) return;
    plc::CImgVisualizer{}.write(db, state, path, options.plot_options);
}

MethodSummary run_method(const plc::PlacementDB& db, const Options& options,
                         const std::string& method) {
    // Each method receives its own copy. No initializer sees another method's result.
    plc::PlacementState state = db.init_state;
    if (state.size() != db.init_state.size()) {
        throw std::runtime_error("Compare state copy has an unexpected size");
    }
    MethodSummary summary;
    summary.method = method;
    const auto begin = std::chrono::steady_clock::now();
    std::cout << "\n---------------INITIAL PLACEMENT: " << method << "---------------\n";
    if (method == "random") {
        plc::RandomPlacement placer(options.seed); placer.run(db, &state);
        std::cout << "    Placed Cells: " << placer.stats().placed_cells
                  << "\n    Skipped Cells: " << placer.stats().skipped_cells << '\n';
    } else if (method == "cluster") {
        plc::ClusteringPlacement placer(options.clustering); placer.run(db, &state);
        const auto& stats = placer.stats();
        std::cout << "    Placed Cells: " << stats.placed_cells << "\n    Skipped Cells: " << stats.skipped_cells
            << "\n    Clusters: " << stats.clusters << "\n    Merges: " << stats.merges
            << "\n    Approximated Nets: " << stats.approximated_nets << '\n';
    } else {
        plc::QuadraticAnalyticalPlacement placer; placer.run(db, &state);
        const auto& stats = placer.stats();
        std::cout << "    Solver Initial HPWL: " << std::fixed << std::setprecision(2) << stats.initial_hpwl
            << "\n    Outer Iterations: " << stats.outer_iterations
            << "\n    Solver Iterations: " << stats.total_solver_iterations
            << "\n    Final X Error: " << std::scientific << stats.final_x_error
            << "\n    Final Y Error: " << stats.final_y_error << '\n';
        if (!std::isfinite(stats.final_x_error) || !std::isfinite(stats.final_y_error) ||
            stats.final_x_error >= plc::QuadraticConfig{}.target_error ||
            stats.final_y_error >= plc::QuadraticConfig{}.target_error) {
            std::cout << "    Solver Warning: residual tolerance not reached; inspect this result.\n";
            summary.status = "incomplete";
        }
    }
    auto ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - begin).count();
    plc::HpwlEvaluator evaluator;
    double hpwl = evaluator.compute_total_hpwl(db, state);
    if (!std::isfinite(hpwl)) throw std::runtime_error("Non-finite HPWL");
    const double overflow = plc::RandomPlacement(options.seed).estimate_overlap_area(db, state);
    summary.hpwl = hpwl;
    summary.overflow_percent = db.movable_cell_area() > 0
        ? 100.0 * overflow / db.movable_cell_area() : 0.0;
    summary.runtime_ms = ms;
    std::cout << "    Movable Cells: " << db.num_movable_cells()
        << "\n    Total HPWL: " << std::fixed << std::setprecision(2) << hpwl
        << "\n    Runtime: " << std::setprecision(3) << ms << " ms"
        << "\n    Bin Overflow Area: " << std::setprecision(2) << overflow
        << "\n    Bin Overflow Ratio: " << summary.overflow_percent << "%\n";
    validate(db, state);
    plot_stage(db, state, options, method, "initial", !options.run_gp && !options.sa);

    // ---- task 5: global placement after initializer ----
    if (options.run_gp) {
        const auto gp_begin = std::chrono::steady_clock::now();
        plc::GlobalPlacementConfig gp_cfg = options.gp_config;
        // Reuse the initializer's seed if not explicitly overridden
        gp_cfg.seed = (gp_cfg.seed == plc::GlobalPlacementConfig{}.seed)
                          ? options.seed : gp_cfg.seed;
        plc::GlobalPlacer gp(gp_cfg);
        gp.run(db, &state);
        const auto gp_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - gp_begin).count();
        const auto& gs = gp.stats();
        std::cout << "\n---------------GLOBAL PLACEMENT---------------\n"
            << "    Iterations: " << gs.iterations
            << "\n    Bin Grid: " << gs.bins_x << " x " << gs.bins_y
            << "\n    Fillers: " << gs.fillers
            << "\n    Initial HPWL: " << std::fixed << std::setprecision(2) << gs.initial_hpwl
            << "\n    Final HPWL: " << std::fixed << std::setprecision(2) << gs.final_hpwl
            << "\n    Initial Density Overflow: " << std::fixed << std::setprecision(2) << gs.initial_overflow * 100.0 << "%"
            << "\n    Final Density Overflow: " << std::fixed << std::setprecision(2) << gs.final_overflow * 100.0 << "%"
            << "\n    Final Lambda: " << std::fixed << std::setprecision(4) << gs.final_lambda
            << "\n    Runtime: " << std::setprecision(3) << gp_ms << " ms\n";
        hpwl = gs.final_hpwl;  // update for SA below
        validate(db, state);
        plot_stage(db, state, options, method, "global", !options.sa);
    }

    export_initial(db, state, options, method);

    if (options.sa) {
        plc::LocalSAPlacer placer({}, options.seed);
        const auto sa_begin = std::chrono::steady_clock::now();
        placer.run(db, &state);
        const auto sa_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sa_begin).count();
        std::cout << "\n---------------LOCAL SA AFTER: " << method
                  << (options.run_gp ? " (post-GP)" : "") << "---------------\n"
            << "    Initial HPWL: " << placer.stats().initial_hpwl
            << "\n    Final HPWL: " << placer.stats().final_hpwl
            << "\n    Accepted Moves: " << placer.stats().accepted_moves << " / " << placer.stats().attempted_moves
            << "\n    Rejected Moves: " << placer.stats().rejected_moves
            << "\n    Skipped Moves: " << placer.stats().skipped_moves
            << "\n    Final Temperature: " << placer.stats().final_temperature
            << "\n    Runtime: " << sa_ms << " ms\n";
        validate(db, state);
        plot_stage(db, state, options, method, "sa", true);
    }

    return summary;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    try {
        const auto options = parse(argc, argv);
        auto result = plc::BookshelfLoader{}.load_from_aux(argv[1]);
        plc::SummaryWriter writer;
        if (!result.report.ok()) { writer.write_validation_report(result.report, std::cout); return 2; }
        writer.write_bookshelf_report(result.db, std::cout);
        if (options.check_only) return 0;
        const std::vector<std::string> methods = options.method == "compare" ?
            std::vector<std::string>{"random", "cluster", "quadratic"} : std::vector<std::string>{options.method};
        std::vector<MethodSummary> summaries;
        summaries.reserve(methods.size());
        for (const auto& method : methods) {
            if (options.method == "compare" && exceeds_compare_budget(result.db, options, method)) {
                std::cout << "\n---------------INITIAL PLACEMENT: " << method << "---------------\n"
                          << "    Compare budget exceeded before start; method skipped.\n";
                MethodSummary summary;
                summary.method = method;
                summary.status = "timeout";
                summaries.push_back(std::move(summary));
                continue;
            }
            try {
                summaries.push_back(run_method(result.db, options, method));
            } catch (const std::exception& error) {
                if (options.method != "compare") throw;
                std::cerr << "    " << method << " incomplete: " << error.what() << '\n';
                MethodSummary summary;
                summary.method = method;
                summary.status = "incomplete";
                summaries.push_back(std::move(summary));
            }
        }
        if (options.method == "compare") print_compare_summary(summaries);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n'; return 3;
    }
    return 0;
}

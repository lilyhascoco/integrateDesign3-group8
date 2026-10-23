#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "algorithms/clustering_placement.h"
#include "algorithms/local_sa_placer.h"
#include "algorithms/quadratic_analytical_placement.h"
#include "algorithms/random_placement.h"
#include "io/summary_writer.h"
#include "metrics/hpwl.h"
#include "parser/bookshelf_loader.h"

namespace {
struct Options {
    std::string method = "random", output;
    unsigned seed = 42;
    bool check_only = false, sa = false;
    plc::ClusteringConfig clustering;
};
void usage(const char* program) {
    std::cerr << "Usage: " << program << " <benchmark.aux> [options]\n"
        << "  --method random|cluster|quadratic|compare (default random)\n"
        << "  --seed N --clusters N --max-cluster-cells N\n"
        << "  --sa | --no-sa (default: no SA)\n"
        << "  --check-only --output result.pl\n"
        << "compare runs all three independently; output becomes result.<method>.pl.\n"
        << "Exports always contain INITIAL placement, before optional SA.\n";
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
        if (arg == "--method") o.method = value();
        else if (arg == "--output") o.output = value();
        else if (arg == "--seed") {
            auto n = number();
            if (n > std::numeric_limits<unsigned>::max()) throw std::invalid_argument("Seed out of range");
            o.seed = static_cast<unsigned>(n);
        } else if (arg == "--clusters") o.clustering.target_clusters = number();
        else if (arg == "--max-cluster-cells") {
            o.clustering.max_cluster_cells = number();
            if (!o.clustering.max_cluster_cells) throw std::invalid_argument("Cluster cap must be positive");
        } else if (arg == "--sa") o.sa = true;
        else if (arg == "--no-sa") o.sa = false;
        else if (arg == "--check-only") o.check_only = true;
        else throw std::invalid_argument("Unknown option: " + arg);
    }
    if (o.method != "random" && o.method != "cluster" && o.method != "quadratic" && o.method != "compare")
        throw std::invalid_argument("Unknown method: " + o.method);
    return o;
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
    std::cout << "    Initial Placement Written To: " << path.string() << '\n';
}
void run_method(const plc::PlacementDB& db, const Options& options, const std::string& method) {
    // Each method receives its own copy. No initializer sees another method's result.
    plc::PlacementState state = db.init_state;
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
            stats.final_y_error >= plc::QuadraticConfig{}.target_error)
            std::cout << "    Solver Warning: residual tolerance not reached; inspect this result.\n";
    }
    const auto ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - begin).count();
    plc::HpwlEvaluator evaluator;
    const double hpwl = evaluator.compute_total_hpwl(db, state);
    if (!std::isfinite(hpwl)) throw std::runtime_error("Non-finite HPWL");
    const double overflow = plc::RandomPlacement(options.seed).estimate_overlap_area(db, state);
    std::cout << "    Movable Cells: " << db.num_movable_cells()
        << "\n    Total HPWL: " << std::fixed << std::setprecision(2) << hpwl
        << "\n    Runtime: " << std::setprecision(3) << ms << " ms"
        << "\n    Bin Overflow Area: " << std::setprecision(2) << overflow
        << "\n    Bin Overflow Ratio: " << (db.movable_cell_area() > 0 ? 100 * overflow / db.movable_cell_area() : 0) << "%\n";
    validate(db, state);
    export_initial(db, state, options, method);
    if (options.sa) {
        plc::LocalSAPlacer placer({}, options.seed);
        const auto sa_begin = std::chrono::steady_clock::now();
        placer.run(db, &state);
        const auto sa_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sa_begin).count();
        std::cout << "\n---------------LOCAL SA AFTER: " << method << "---------------\n"
            << "    Initial HPWL: " << placer.stats().initial_hpwl
            << "\n    Final HPWL: " << placer.stats().final_hpwl
            << "\n    Accepted Moves: " << placer.stats().accepted_moves << " / " << placer.stats().attempted_moves
            << "\n    Rejected Moves: " << placer.stats().rejected_moves
            << "\n    Skipped Moves: " << placer.stats().skipped_moves
            << "\n    Final Temperature: " << placer.stats().final_temperature
            << "\n    Runtime: " << sa_ms << " ms\n";
        validate(db, state);
    }
}
} // namespace
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
        for (const auto& method : methods) run_method(result.db, options, method);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n'; return 3;
    }
    return 0;
}

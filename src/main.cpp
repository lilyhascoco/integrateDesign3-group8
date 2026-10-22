#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "algorithms/global_placement.h"
#include "algorithms/local_sa_placer.h"
#include "algorithms/random_placement.h"
#include "io/summary_writer.h"
#include "metrics/hpwl.h"
#include "parser/bookshelf_loader.h"

namespace {

void PrintUsage(const char* program) {
    std::cerr
        << "Usage: " << program << " <benchmark.aux> [options]\n"
        << "Options:\n"
        << "  --check-only            Parse and validate the benchmark only\n"
        << "  --iterations <N>        Maximum global placement iterations (default 200)\n"
        << "  --target-density <v>    Target density in (0, 1] (default 0.9)\n"
        << "  --target-overflow <v>   Stop when density overflow is below v (default 0.1)\n"
        << "  --bins <N>              Upper bound of the m x m bin grid (default 512)\n"
        << "  --seed <N>              Random seed (default 2024)\n"
        << "  --output <file.pl>      Write the global placement result\n"
        << "  --local-sa              Run local simulated annealing after global placement\n"
        << "  --quiet                 Disable the per-iteration global placement log\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        PrintUsage(argv[0]);
        return 1;
    }

    const std::string aux_path = argv[1];
    bool check_only = false;
    bool run_local_sa = false;
    std::string output_path;

    plc::GlobalPlacementConfig gp_config;
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--check-only") {
            check_only = true;
        } else if (arg == "--local-sa") {
            run_local_sa = true;
        } else if (arg == "--quiet") {
            gp_config.verbose = false;
        } else if (arg == "--iterations" && i + 1 < argc) {
            gp_config.max_iterations = std::stoul(argv[++i]);
        } else if (arg == "--target-density" && i + 1 < argc) {
            gp_config.target_density = std::stod(argv[++i]);
        } else if (arg == "--target-overflow" && i + 1 < argc) {
            gp_config.target_overflow = std::stod(argv[++i]);
        } else if (arg == "--bins" && i + 1 < argc) {
            gp_config.max_bins = std::stoul(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc) {
            gp_config.seed = static_cast<unsigned>(std::stoul(argv[++i]));
        } else if (arg == "--output" && i + 1 < argc) {
            output_path = argv[++i];
        } else {
            PrintUsage(argv[0]);
            return 1;
        }
    }

    plc::BookshelfLoader loader;
    plc::LoadResult result = loader.load_from_aux(aux_path);

    plc::SummaryWriter writer;
    if (!result.report.ok()) {
        writer.write_validation_report(result.report, std::cout);
        return 2;
    }

    writer.write_bookshelf_report(result.db, std::cout);
    if (check_only) {
        return 0;
    }

    plc::PlacementState state = result.db.init_state;
    plc::RandomPlacement random_placer(gp_config.seed);

    const auto random_start_time = std::chrono::steady_clock::now();
    random_placer.run(result.db, &state);
    const auto random_end_time = std::chrono::steady_clock::now();

    plc::HpwlEvaluator hpwl_evaluator;
    const double random_hpwl = hpwl_evaluator.compute_total_hpwl(result.db, state);
    const auto random_elapsed_ms =
        std::chrono::duration<double, std::milli>(random_end_time - random_start_time).count();
    const auto& random_stats = random_placer.stats();

    std::cout << "\n---------------RANDOM PLACEMENT---------------\n";
    std::cout << "    Movable Cells: " << random_stats.movable_cells << '\n';
    std::cout << "    Randomly Placed Cells: " << random_stats.placed_cells << '\n';
    std::cout << "    Skipped Cells: " << random_stats.skipped_cells << '\n';
    std::cout << "    Total HPWL: " << std::fixed << std::setprecision(2) << random_hpwl << '\n';
    std::cout << "    Runtime: " << std::fixed << std::setprecision(3) << random_elapsed_ms
              << " ms\n";
    std::cout << "    Estimated Overlap Area: " << std::fixed << std::setprecision(2)
              << random_stats.estimated_overlap_area << '\n';
    std::cout << "    Estimated Overlap Ratio: " << std::fixed << std::setprecision(2)
              << random_stats.estimated_overlap_ratio * 100.0 << "%\n";

    plc::GlobalPlacer global_placer(gp_config);
    const auto gp_start_time = std::chrono::steady_clock::now();
    global_placer.run(result.db, &state);
    const auto gp_end_time = std::chrono::steady_clock::now();
    const auto& gp_stats = global_placer.stats();
    const auto gp_elapsed_ms =
        std::chrono::duration<double, std::milli>(gp_end_time - gp_start_time).count();

    std::cout << "\n---------------GLOBAL PLACEMENT---------------\n";
    std::cout << "    Iterations: " << gp_stats.iterations << '\n';
    std::cout << "    Bin Grid: " << gp_stats.bins_x << " x " << gp_stats.bins_y << '\n';
    std::cout << "    Fillers: " << gp_stats.fillers << '\n';
    std::cout << "    Initial HPWL: " << std::fixed << std::setprecision(2)
              << gp_stats.initial_hpwl << '\n';
    std::cout << "    Final HPWL: " << std::fixed << std::setprecision(2)
              << gp_stats.final_hpwl << '\n';
    std::cout << "    Initial Density Overflow: " << std::fixed << std::setprecision(2)
              << gp_stats.initial_overflow * 100.0 << "%\n";
    std::cout << "    Final Density Overflow: " << std::fixed << std::setprecision(2)
              << gp_stats.final_overflow * 100.0 << "%\n";
    std::cout << "    Final Lambda: " << std::fixed << std::setprecision(4)
              << gp_stats.final_lambda << '\n';
    std::cout << "    Runtime: " << std::fixed << std::setprecision(3) << gp_elapsed_ms
              << " ms\n";

    if (!output_path.empty()) {
        std::ofstream out(output_path);
        if (out) {
            writer.write_placement(result.db, state, out);
            std::cout << "    Placement Written To: " << output_path << '\n';
        } else {
            std::cerr << "Failed to open output file: " << output_path << '\n';
        }
    }

    if (!run_local_sa) {
        return 0;
    }

    plc::LocalSAConfig sa_config;
    plc::LocalSAPlacer local_sa_placer(sa_config, gp_config.seed);
    const auto sa_start_time = std::chrono::steady_clock::now();
    local_sa_placer.run(result.db, &state);
    const auto sa_end_time = std::chrono::steady_clock::now();

    const auto& sa_stats = local_sa_placer.stats();
    const auto sa_elapsed_ms =
        std::chrono::duration<double, std::milli>(sa_end_time - sa_start_time).count();
    std::cout << "\n---------------LOCAL SIMULATED ANNEALING---------------\n";
    std::cout << "    Initial HPWL: " << std::fixed << std::setprecision(2)
              << sa_stats.initial_hpwl << '\n';
    std::cout << "    Final HPWL: " << std::fixed << std::setprecision(2)
              << sa_stats.final_hpwl << '\n';
    std::cout << "    Accepted Moves: " << sa_stats.accepted_moves << " / "
              << sa_stats.attempted_moves << '\n';
    std::cout << "    Rejected Moves: " << sa_stats.rejected_moves << '\n';
    std::cout << "    Skipped Moves: " << sa_stats.skipped_moves << '\n';
    std::cout << "    Final Temperature: " << std::fixed << std::setprecision(4)
              << sa_stats.final_temperature << '\n';
    std::cout << "    Runtime: " << std::fixed << std::setprecision(3) << sa_elapsed_ms
              << " ms\n";
    return 0;
}

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

#include "algorithms/local_sa_placer.h"
#include "algorithms/quadratic_analytical_placement.h"
#include "algorithms/random_placement.h"
#include "io/summary_writer.h"
#include "metrics/hpwl.h"
#include "parser/bookshelf_loader.h"

namespace {

void PrintUsage(const char* program) {
    std::cerr << "Usage: " << program << " <benchmark.aux> [--check-only]\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        PrintUsage(argv[0]);
        return 1;
    }

    const std::string aux_path = argv[1];
    plc::BookshelfLoader loader;
    plc::LoadResult result = loader.load_from_aux(aux_path);

    plc::SummaryWriter writer;
    if (!result.report.ok()) {
        writer.write_validation_report(result.report, std::cout);
        return 2;
    }

    writer.write_bookshelf_report(result.db, std::cout);

    plc::PlacementState random_state = result.db.init_state;
    plc::RandomPlacement random_placer;

    const auto start_time = std::chrono::steady_clock::now();
    random_placer.run(result.db, &random_state);
    const auto end_time = std::chrono::steady_clock::now();

    plc::HpwlEvaluator hpwl_evaluator;
    const double total_hpwl = hpwl_evaluator.compute_total_hpwl(result.db, random_state);
    const auto elapsed_ms =
        std::chrono::duration<double, std::milli>(end_time - start_time).count();
    const auto& stats = random_placer.stats();

    std::cout << "\n---------------RANDOM PLACEMENT---------------\n";
    std::cout << "    Movable Cells: " << stats.movable_cells << '\n';
    std::cout << "    Randomly Placed Cells: " << stats.placed_cells << '\n';
    std::cout << "    Skipped Cells: " << stats.skipped_cells << '\n';
    std::cout << "    Total HPWL: " << std::fixed << std::setprecision(2) << total_hpwl << '\n';
    std::cout << "    Runtime: " << std::fixed << std::setprecision(3) << elapsed_ms << " ms\n";
    std::cout << "    Estimated Overlap Area: " << std::fixed << std::setprecision(2)
              << stats.estimated_overlap_area << '\n';
    std::cout << "    Estimated Overlap Ratio: " << std::fixed << std::setprecision(2)
              << stats.estimated_overlap_ratio * 100.0 << "%\n";

    // Local SA starts from the random placement above. Replace random_state
    // with any existing placement state to use it as a different initializer.
    plc::LocalSAConfig sa_config;
    plc::LocalSAPlacer local_sa_placer(sa_config);
    const auto sa_start_time = std::chrono::steady_clock::now();
    local_sa_placer.run(result.db, &random_state);
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

    // Quadratic Analytical Placement —— 独立运行，不影响 random_state / SA 的结果
    plc::PlacementState quadratic_state = result.db.init_state;
    plc::QuadraticAnalyticalPlacement quadratic_placer;
    const auto qa_start_time = std::chrono::steady_clock::now();
    quadratic_placer.run(result.db, &quadratic_state);
    const auto qa_end_time = std::chrono::steady_clock::now();
    const auto& qa_stats = quadratic_placer.stats();
    const auto qa_elapsed_ms =
        std::chrono::duration<double, std::milli>(qa_end_time - qa_start_time).count();

    std::cout << "\n---------------QUADRATIC ANALYTICAL PLACEMENT---------------\n";
    std::cout << "    Movable Cells: " << qa_stats.movable_cells << '\n';
    std::cout << "    Outer Iterations: " << qa_stats.outer_iterations << '\n';
    std::cout << "    Total Solver Iterations: " << qa_stats.total_solver_iterations << '\n';
    std::cout << "    Initial HPWL: " << std::fixed << std::setprecision(2)
              << qa_stats.initial_hpwl << '\n';
    std::cout << "    Final HPWL: " << std::fixed << std::setprecision(2)
              << qa_stats.final_hpwl << '\n';
    std::cout << "    Final X Error: " << std::fixed << std::setprecision(6)
              << qa_stats.final_x_error << '\n';
    std::cout << "    Final Y Error: " << std::fixed << std::setprecision(6)
              << qa_stats.final_y_error << '\n';
    std::cout << "    Runtime: " << std::fixed << std::setprecision(3) << qa_elapsed_ms
              << " ms\n";

    return 0;
}
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include "algorithms/clustering_placement.h"

#include "algorithms/local_sa_placer.h"
#include "algorithms/random_placement.h"
#include "io/summary_writer.h"
#include "metrics/hpwl.h"
#include "parser/bookshelf_loader.h"

namespace {

void PrintUsage(const char* program) {
    std::cerr << "Usage: " << program << " <benchmark.aux> [--method random|cluster|compare] [--seed N] [--clusters N] [--max-cluster-cells N] [--no-sa] [--output result.pl] [--check-only]\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        PrintUsage(argv[0]);
        return 1;
    }

    std::string method="random", output;
    unsigned seed=42;
    bool check_only=false, run_sa=true;
    plc::ClusteringConfig cluster_config;
    try {
        for(int i=2;i<argc;++i) {
            const std::string arg=argv[i];
            auto value=[&]() -> std::string {
                if(i+1>=argc) throw std::invalid_argument("Missing value for "+arg);
                return argv[++i];
            };
            auto number=[&]() -> std::size_t {
                const auto v=value(); std::size_t end=0;
                if(v.empty() || v[0]=='-') throw std::invalid_argument("Invalid number");
                auto n=std::stoull(v,&end);
                if(end!=v.size()) throw std::invalid_argument("Invalid number");
                return n;
            };
            if(arg=="--method") method=value();
            else if(arg=="--seed") seed=static_cast<unsigned>(number());
            else if(arg=="--clusters") cluster_config.target_clusters=number();
            else if(arg=="--max-cluster-cells") {
                cluster_config.max_cluster_cells=number();
                if(!cluster_config.max_cluster_cells) throw std::invalid_argument("Cluster cap must be positive");
            }
            else if(arg=="--output") output=value();
            else if(arg=="--check-only") check_only=true;
            else if(arg=="--no-sa") run_sa=false;
            else throw std::invalid_argument("Unknown option: "+arg);
        }
        if(method!="random" && method!="cluster" && method!="compare")
            throw std::invalid_argument("Method must be random, cluster or compare");
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
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
    if(check_only) return 0;
    plc::RandomPlacement random_placer(seed);

    const auto start_time = std::chrono::steady_clock::now();
    if(method!="cluster") random_placer.run(result.db, &random_state);
    const auto end_time = std::chrono::steady_clock::now();

    plc::HpwlEvaluator hpwl_evaluator;
    const double total_hpwl = hpwl_evaluator.compute_total_hpwl(result.db, random_state);
    const auto elapsed_ms =
        std::chrono::duration<double, std::milli>(end_time - start_time).count();
    const auto& stats = random_placer.stats();

    if(method!="cluster") {
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

    }
    if(method!="random") {
        plc::ClusteringPlacement placer(cluster_config);
        const auto begin=std::chrono::steady_clock::now();
        try {placer.run(result.db,&random_state);}
        catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 3;}
        const double overflow=random_placer.estimate_overlap_area(result.db,random_state);
        const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count();
        const auto& cs=placer.stats();
        std::cout<<"\n---------------CLUSTERING PLACEMENT---------------\n"
            <<"    Movable Cells: "<<cs.movable_cells<<"\n    Placed Cells: "<<cs.placed_cells
            <<"\n    Skipped Cells: "<<cs.skipped_cells<<"\n    Clusters: "<<cs.clusters
            <<"\n    Merges: "<<cs.merges<<"\n    Approximated Nets: "<<cs.approximated_nets
            <<"\n    Total HPWL: "<<std::fixed<<std::setprecision(2)<<hpwl_evaluator.compute_total_hpwl(result.db,random_state)
            <<"\n    Runtime: "<<std::setprecision(3)<<ms<<" ms\n"
            <<"    Estimated Overlap Area: "<<std::setprecision(2)<<overflow<<'\n'
            <<"    Estimated Overlap Ratio: "<<(result.db.movable_cell_area()>0 ? 100*overflow/result.db.movable_cell_area() : 0)<<"%\n";
    }
    std::size_t violations=0;
    const auto& core=result.db.design.core;
    for(const auto& cell:result.db.design.cells) {
        double x=random_state.x(cell.id),y=random_state.y(cell.id);
        if(!std::isfinite(x)||!std::isfinite(y)) ++violations;
        if(cell.is_fixed() && (x!=cell.init_x || y!=cell.init_y || random_state.orient(cell.id)!=cell.orient)) ++violations;
        if(cell.is_movable() && (x<core.lx-1e-8 || y<core.ly-1e-8 || x+cell.width>core.ux+1e-8 || y+cell.height>core.uy+1e-8)) ++violations;
    }
    std::cout<<"    Boundary/Fixed Violations: "<<violations<<'\n';
    if(violations) return 3;
    if(!output.empty()) {
        std::ofstream out(output);
        if(!out) {std::cerr<<"Cannot open output: "<<output<<'\n';return 4;}
        out<<"UCLA pl 1.0\n"<<std::setprecision(17);
        const char* orientations[]={"N","S","E","W","FN","FS","FE","FW","N"};
        for(const auto& cell:result.db.design.cells) {
            out<<cell.name<<' '<<random_state.x(cell.id)<<' '<<random_state.y(cell.id)
               <<" : "<<orientations[static_cast<int>(random_state.orient(cell.id))];
            if(cell.is_fixed()) out<<(cell.place_status==plc::PlaceStatus::kFixedNI ? " /FIXED_NI" : " /FIXED");
            out<<'\n';
        }
        out.close();
        if(!out) {std::cerr<<"Failed to write output\n";return 4;}
    }
    if(!run_sa) return 0;
    // Continue from the selected initializer (cluster in compare mode).
    plc::LocalSAConfig sa_config;
    plc::LocalSAPlacer local_sa_placer(sa_config, seed);
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
    return 0;
}

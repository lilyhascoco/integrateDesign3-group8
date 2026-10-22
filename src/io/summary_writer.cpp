#include "io/summary_writer.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <string>

#include "metrics/area_stats.h"

namespace plc {

namespace {

std::string file_name(const std::string& path) {
    return std::filesystem::path(path).filename().string();
}

std::string format_number(double value) {
    if (std::fabs(value - std::round(value)) < 1e-6) {
        return std::to_string(static_cast<long long>(std::llround(value)));
    }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << value;
    return oss.str();
}

std::string format_percent(double numerator, double denominator) {
    if (denominator <= 0.0) {
        return "0.00";
    }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << (numerator * 100.0 / denominator);
    return oss.str();
}

std::string orient_token(Orient orient) {
    switch (orient) {
        case Orient::kN: return "N";
        case Orient::kS: return "S";
        case Orient::kE: return "E";
        case Orient::kW: return "W";
        case Orient::kFN: return "FN";
        case Orient::kFS: return "FS";
        case Orient::kFE: return "FE";
        case Orient::kFW: return "FW";
        default: return "N";
    }
}

double clipped_fixed_area_in_core(const PlacementDB& db) {
    double area = 0.0;
    const auto& core = db.design.core;
    for (const auto& cell : db.design.cells) {
        if (!cell.is_fixed()) {
            continue;
        }
        const double lx = std::max(cell.init_x, core.lx);
        const double ly = std::max(cell.init_y, core.ly);
        const double ux = std::min(cell.init_x + cell.width, core.ux);
        const double uy = std::min(cell.init_y + cell.height, core.uy);
        if (ux > lx && uy > ly) {
            area += (ux - lx) * (uy - ly);
        }
    }
    return area;
}

}  // namespace

void SummaryWriter::write_bookshelf_report(const PlacementDB& db, std::ostream& os) const {
    const double core_area = (db.design.core.ux - db.design.core.lx) *
                             (db.design.core.uy - db.design.core.ly);
    const double total_area = db.total_cell_area();
    const double movable_area = db.movable_cell_area();
    const double cell_area = movable_area;
    const double fixed_area = total_area - movable_area;
    const double fixed_area_in_core = clipped_fixed_area_in_core(db);

    std::size_t max_net_degree = 0;
    std::size_t degree_1 = 0;
    std::size_t degree_2 = 0;
    std::size_t degree_3_to_10 = 0;
    std::size_t degree_11_to_100 = 0;
    std::size_t degree_gt_100 = 0;
    for (const auto& net : db.design.nets) {
        const std::size_t degree = net.degree();
        max_net_degree = std::max(max_net_degree, degree);
        if (degree == 1) {
            ++degree_1;
        } else if (degree == 2) {
            ++degree_2;
        } else if (degree >= 3 && degree <= 10) {
            ++degree_3_to_10;
        } else if (degree >= 11 && degree <= 100) {
            ++degree_11_to_100;
        } else if (degree > 100) {
            ++degree_gt_100;
        }
    }

    os << "Use BOOKSHELF placement format\n";
    os << "ReadAUXFile: " << db.design.info.aux_path << '\n';
    os << "ReadSCLFile: " << file_name(db.design.info.scl_path) << '\n';
    os << "    CoreRegion: lower left: (" << format_number(db.design.core.lx) << ","
       << format_number(db.design.core.ly) << ") to upper right: ("
       << format_number(db.design.core.ux) << "," << format_number(db.design.core.uy) << ")\n";
    os << "    NumRows: " << db.num_rows() << '\n';
    os << "ReadNodesFile: " << file_name(db.design.info.nodes_path) << '\n';
    os << "    NumModules: " << db.num_cells() << '\n';
    os << "    NumNodes: " << db.num_movable_cells() << '\n';
    os << "    Terminals: " << db.num_fixed_cells() << '\n';
    os << "ReadNetsFile: " << file_name(db.design.info.nets_path) << '\n';
    os << "    Nets: " << db.num_nets() << '\n';
    os << "    Pins: " << db.num_pins() << '\n';
    os << "    Max net degree= " << max_net_degree << '\n';
    os << "    total pin number= " << db.num_pins() << '\n';
    os << "ReadPLFile: " << file_name(db.design.info.pl_path) << "\n\n";

    os << "---------------SUMMARIES---------------\n";
    os << "    Area:\n";
    os << "    Core Area: " << format_number(core_area) << '\n';
    os << "    Cell Area: " << format_number(cell_area)
       << " (cellArea / coreArea = " << format_percent(cell_area, core_area) << "%)\n";
    os << "    Movable Area: " << format_number(movable_area)
       << " (movableArea / coreArea = " << format_percent(movable_area, core_area) << "%)\n";
    os << "    Fixed Area: " << format_number(fixed_area)
       << " (fixedArea / coreArea = " << format_percent(fixed_area, core_area) << "%)\n";
    os << "    Fixed Area in Core: " << format_number(fixed_area_in_core)
       << " (fixedAreaInCore / coreArea = "
       << format_percent(fixed_area_in_core, core_area) << "%)\n";
    os << "    There are " << degree_1 << " nets (has 1 pins)\n";
    os << "    There are " << degree_2 << " nets (has 2 pins)\n";
    os << "    There are " << degree_3_to_10 << " nets (has 3-10 pins)\n";
    os << "    There are " << degree_11_to_100 << " nets (has 11-100 pins)\n";
    os << "    There are " << degree_gt_100 << " nets (has >100 pins)\n";
}

void SummaryWriter::write_design_summary(const PlacementDB& db, std::ostream& os) const {
    AreaStatsCalculator area_calc;
    const AreaStats area = area_calc.compute(db);

    os << "Design: " << db.design.info.design_name << '\n';
    os << "Cells: " << db.num_cells() << " (movable=" << db.num_movable_cells()
       << ", fixed=" << db.num_fixed_cells() << ")\n";
    os << "Nets: " << db.num_nets() << ", Pins: " << db.num_pins() << '\n';
    os << "Rows: " << db.num_rows() << '\n';
    os << "Core: [" << db.design.core.lx << ", " << db.design.core.ly << "] - ["
       << db.design.core.ux << ", " << db.design.core.uy << "]\n";
    os << std::fixed << std::setprecision(4);
    os << "Area: total=" << area.total_cell_area << ", movable=" << area.movable_area
       << ", fixed=" << area.fixed_area << ", util=" << area.utilization << '\n';
}

void SummaryWriter::write_validation_report(const ValidationReport& report,
                                            std::ostream& os) const {
    os << "Validation: " << (report.ok() ? "OK" : "FAILED") << '\n';
    if (!report.errors.empty()) {
        os << "Errors:\n";
        for (const auto& msg : report.errors) {
            os << "  - " << msg << '\n';
        }
    }
    if (!report.warnings.empty()) {
        os << "Warnings:\n";
        for (const auto& msg : report.warnings) {
            os << "  - " << msg << '\n';
        }
    }
}

void SummaryWriter::write_placement(const PlacementDB& db, const PlacementState& state,
                                    std::ostream& os) const {
    os << "UCLA pl 1.0\n";
    os << "# Generated by the electrostatics-based global placer\n";
    os << std::fixed << std::setprecision(3);
    for (const auto& cell : db.design.cells) {
        os << cell.name << ' ' << state.x(cell.id) << ' ' << state.y(cell.id)
           << " : " << orient_token(state.orient(cell.id));
        if (!cell.is_movable()) {
            os << " /FIXED";
        }
        os << '\n';
    }
}

}  // namespace plc

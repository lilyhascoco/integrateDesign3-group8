#include "parser/bookshelf_loader.h"

#include <filesystem>

#include "parser/aux_parser.h"
#include "parser/nets_parser.h"
#include "parser/nodes_parser.h"
#include "parser/pl_parser.h"
#include "parser/scl_parser.h"
#include "parser/wts_parser.h"

namespace plc {

LoadResult BookshelfLoader::load_from_aux(const std::string& aux_path) const {
    LoadResult result;

    BookshelfInputSpec spec;
    AuxParser aux_parser;
    if (!aux_parser.parse(aux_path, &spec, &result.report)) {
        return result;
    }

    result.db.design.info.design_name = spec.benchmark_name;
    result.db.design.info.benchmark_dir = spec.benchmark_dir;
    result.db.design.info.aux_path = spec.aux_path;
    result.db.design.info.scl_path = spec.scl_path;
    result.db.design.info.nodes_path = spec.nodes_path;
    result.db.design.info.nets_path = spec.nets_path;
    result.db.design.info.pl_path = spec.pl_path;
    result.db.design.info.wts_path = spec.wts_path;
    if (!run_parsers(spec, &result.db, &result.report)) {
        return result;
    }
    build_initial_state(&result.db);
    run_post_checks(spec, &result.db, &result.report);
    return result;
}

bool BookshelfLoader::run_parsers(const BookshelfInputSpec& spec, PlacementDB* db, ValidationReport* report) const {
    SclParser scl_parser;
    NodesParser nodes_parser;
    PlParser pl_parser;
    NetsParser nets_parser;
    WtsParser wts_parser;

    bool ok = true;
    if (!spec.scl_path.empty()) {
        ok = scl_parser.parse(spec.scl_path, db, report) && ok;
    }
    ok = nodes_parser.parse(spec.nodes_path, db, report) && ok;
    if (!spec.pl_path.empty()) {
        ok = pl_parser.parse(spec.pl_path, db, report) && ok;
    }
    ok = nets_parser.parse(spec.nets_path, db, report) && ok;
    if (!spec.wts_path.empty()) {
        ok = wts_parser.parse(spec.wts_path, db, report) && ok;
    }
    return ok;
}

void BookshelfLoader::build_initial_state(PlacementDB* db) const {
    db->init_state.resize(db->num_cells());
    for (const auto& cell : db->design.cells) {
        db->init_state.set_pos(cell.id, cell.init_x, cell.init_y);
        db->init_state.set_orient(cell.id, cell.orient);
    }
}

void BookshelfLoader::run_post_checks(const BookshelfInputSpec& spec, PlacementDB* db, ValidationReport* report) const {
    for (const auto& cell : db->design.cells) {
        if (cell.place_status == PlaceStatus::kUnplaced) {
            report->add_warning("Cell remains unplaced after parsing: " + cell.name);
        }
        if (!db->design.rows.empty() && cell.is_movable()) {
            const bool x_ok = cell.init_x >= db->design.core.lx && cell.init_x <= db->design.core.ux;
            const bool y_ok = cell.init_y >= db->design.core.ly && cell.init_y <= db->design.core.uy;
            if (!x_ok || !y_ok) {
                report->add_warning("Movable cell initial position outside core: " + cell.name);
            }
        }
    }

    if (!std::filesystem::exists(spec.nodes_path)) {
        report->add_error("nodes file does not exist on disk");
    }
    if (!std::filesystem::exists(spec.nets_path)) {
        report->add_error("nets file does not exist on disk");
    }
}

}  // namespace plc

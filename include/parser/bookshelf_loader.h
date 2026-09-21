#pragma once

#include <string>

#include "database/placement_db.h"
#include "database/validation_report.h"
#include "parser/bookshelf_input_spec.h"

namespace plc {

struct LoadResult {
    PlacementDB db;
    ValidationReport report;
};

class BookshelfLoader {
public:
    LoadResult load_from_aux(const std::string& aux_path) const;

private:
    bool run_parsers(const BookshelfInputSpec& spec, PlacementDB* db, ValidationReport* report) const;
    void build_initial_state(PlacementDB* db) const;
    void run_post_checks(const BookshelfInputSpec& spec, PlacementDB* db, ValidationReport* report) const;
};

}  // namespace plc

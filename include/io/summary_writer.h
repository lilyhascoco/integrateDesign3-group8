#pragma once

#include <ostream>

#include "database/placement_db.h"
#include "database/validation_report.h"

namespace plc {

class SummaryWriter {
public:
    void write_bookshelf_report(const PlacementDB& db, std::ostream& os) const;
    void write_design_summary(const PlacementDB& db, std::ostream& os) const;
    void write_validation_report(const ValidationReport& report, std::ostream& os) const;
};

}  // namespace plc

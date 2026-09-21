#include "database/validation_report.h"

namespace plc {

void ValidationReport::add_error(const std::string& msg) {
    errors.push_back(msg);
}

void ValidationReport::add_warning(const std::string& msg) {
    warnings.push_back(msg);
}

}  // namespace plc

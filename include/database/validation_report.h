#pragma once

#include <string>
#include <vector>

namespace plc {

struct ValidationReport {
    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    bool ok() const { return errors.empty(); }

    void add_error(const std::string& msg);
    void add_warning(const std::string& msg);
};

}  // namespace plc

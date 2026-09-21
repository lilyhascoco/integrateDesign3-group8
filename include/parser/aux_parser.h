#pragma once

#include <string>

#include "database/validation_report.h"
#include "parser/bookshelf_input_spec.h"

namespace plc {

class AuxParser {
public:
    bool parse(const std::string& aux_path, BookshelfInputSpec* spec, ValidationReport* report) const;
};

}  // namespace plc

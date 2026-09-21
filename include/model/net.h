#pragma once

#include <string>
#include <vector>

#include "common/types.h"

namespace plc {

struct Net {
    NetId id = 0;
    std::string name;
    double weight = 1.0;
    std::vector<PinId> pin_ids;

    std::size_t degree() const { return pin_ids.size(); }
};

}  // namespace plc

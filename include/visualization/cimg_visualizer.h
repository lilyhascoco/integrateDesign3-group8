#pragma once

#include <cstddef>
#include <string>

#include "database/placement_state.h"
#include "database/placement_db.h"

namespace plc {

struct PlotOptions {
    std::size_t width = 1600;
    std::size_t margin = 24;
    bool draw_rows = true;
    bool draw_nets = false;
};

// Renders a PlacementState without modifying placement coordinates or orientation.
class CImgVisualizer {
public:
    static bool is_available();
    static bool png_supported();
    static std::string default_extension();

    void write(const PlacementDB& db,
               const PlacementState& state,
               const std::string& output_path,
               const PlotOptions& options) const;
};

}  // namespace plc

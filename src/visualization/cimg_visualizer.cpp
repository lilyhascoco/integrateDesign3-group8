#include "visualization/cimg_visualizer.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#if PLC_HAS_CIMG
#define cimg_display 0
#if PLC_HAS_PNG
#define cimg_use_png 1
#endif
#include <CImg.h>
#endif

namespace plc {
namespace {

constexpr std::size_t kMinimumImageDimension = 64;
constexpr std::size_t kMaximumImageHeight = 4096;
constexpr std::size_t kMaximumNetsToDraw = 4000;
constexpr std::size_t kMaximumNetDegreeToDraw = 32;

struct Region {
    double lx = 0.0;
    double ly = 0.0;
    double ux = 0.0;
    double uy = 0.0;
};

bool finite(double value) {
    return std::isfinite(value);
}

void include_rect(Region* region, double lx, double ly, double ux, double uy) {
    if (!finite(lx) || !finite(ly) || !finite(ux) || !finite(uy) || ux < lx || uy < ly) {
        return;
    }
    region->lx = std::min(region->lx, lx);
    region->ly = std::min(region->ly, ly);
    region->ux = std::max(region->ux, ux);
    region->uy = std::max(region->uy, uy);
}

Region chip_region(const PlacementDB& db, const PlacementState& state) {
    const auto& core = db.design.core;
    if (!finite(core.lx) || !finite(core.ly) || !finite(core.ux) || !finite(core.uy) ||
        core.ux <= core.lx || core.uy <= core.ly) {
        throw std::runtime_error("Cannot plot: invalid core region");
    }

    Region region{core.lx, core.ly, core.ux, core.uy};
    for (const auto& cell : db.design.cells) {
        const double x = state.x(cell.id);
        const double y = state.y(cell.id);
        if (cell.width <= 0.0 || cell.height <= 0.0) {
            continue;
        }
        include_rect(&region, x, y, x + cell.width, y + cell.height);
    }
    return region;
}

#if PLC_HAS_CIMG
int to_pixel(double value) {
    const double lower = static_cast<double>(std::numeric_limits<int>::min() / 2);
    const double upper = static_cast<double>(std::numeric_limits<int>::max() / 2);
    return static_cast<int>(std::lround(std::max(lower, std::min(upper, value))));
}

struct Transform {
    Region region;
    int image_height = 0;
    int margin = 0;
    double scale = 1.0;

    int x(double x) const {
        return to_pixel(static_cast<double>(margin) + (x - region.lx) * scale);
    }

    // Placement coordinates are bottom-up; CImg coordinates are top-down.
    int y(double y) const {
        return to_pixel(static_cast<double>(image_height - margin) - (y - region.ly) * scale);
    }
};

std::filesystem::path normalized_path(const std::string& output_path) {
    std::filesystem::path path(output_path);
    if (path.empty()) {
        throw std::runtime_error("Cannot plot: empty output path");
    }
    if (!path.has_extension()) {
        path += CImgVisualizer::default_extension();
    }
    const std::string extension = path.extension().string();
    if (extension != ".png" && extension != ".bmp") {
        throw std::runtime_error("Cannot plot: output extension must be .png or .bmp");
    }
#if !PLC_HAS_PNG
    if (extension == ".png") {
        throw std::runtime_error(
            "Cannot write PNG: libpng was not found at build time; use a .bmp path or install libpng-dev");
    }
#endif
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    return path;
}

void draw_cell(cimg_library::CImg<unsigned char>* image,
               const Cell& cell,
               const PlacementState& state,
               const Transform& transform) {
    const double x0 = state.x(cell.id);
    const double y0 = state.y(cell.id);
    const double x1 = x0 + cell.width;
    const double y1 = y0 + cell.height;
    if (!finite(x0) || !finite(y0) || !finite(x1) || !finite(y1) ||
        cell.width <= 0.0 || cell.height <= 0.0) {
        return;
    }

    const int px0 = transform.x(x0);
    const int px1 = transform.x(x1);
    const int py0 = transform.y(y0);
    const int py1 = transform.y(y1);
    const int left = std::min(px0, px1);
    const int right = std::max(px0, px1);
    const int top = std::min(py0, py1);
    const int bottom = std::max(py0, py1);
    const unsigned char movable[] = {49, 130, 189};
    const unsigned char fixed[] = {214, 39, 40};
    const unsigned char macro[] = {128, 76, 167};
    const unsigned char* color = cell.type == CellType::kMacro ? macro :
                                 cell.is_fixed() ? fixed : movable;
    image->draw_rectangle(left, top, right, bottom, color, 0.92F);
}

void draw_nets(cimg_library::CImg<unsigned char>* image,
               const PlacementDB& db,
               const PlacementState& state,
               const Transform& transform) {
    const unsigned char net_color[] = {150, 150, 150};
    std::size_t drawn = 0;
    for (const auto& net : db.design.nets) {
        if (drawn >= kMaximumNetsToDraw) {
            break;
        }
        if (net.degree() < 2 || net.degree() > kMaximumNetDegreeToDraw) {
            continue;
        }

        const Pin& root_pin = db.pin(net.pin_ids.front());
        const Cell& root_cell = db.cell(root_pin.cell_id);
        const double root_x = state.x(root_cell.id) + root_cell.width * 0.5 + root_pin.offset_x;
        const double root_y = state.y(root_cell.id) + root_cell.height * 0.5 + root_pin.offset_y;
        if (!finite(root_x) || !finite(root_y)) {
            continue;
        }

        for (std::size_t i = 1; i < net.pin_ids.size(); ++i) {
            const Pin& pin = db.pin(net.pin_ids[i]);
            const Cell& cell = db.cell(pin.cell_id);
            const double x = state.x(cell.id) + cell.width * 0.5 + pin.offset_x;
            const double y = state.y(cell.id) + cell.height * 0.5 + pin.offset_y;
            if (finite(x) && finite(y)) {
                image->draw_line(transform.x(root_x), transform.y(root_y),
                                 transform.x(x), transform.y(y), net_color, 0.25F);
            }
        }
        ++drawn;
    }
}
#endif

}  // namespace

bool CImgVisualizer::is_available() {
#if PLC_HAS_CIMG
    return true;
#else
    return false;
#endif
}

bool CImgVisualizer::png_supported() {
#if PLC_HAS_CIMG && PLC_HAS_PNG
    return true;
#else
    return false;
#endif
}

std::string CImgVisualizer::default_extension() {
    return png_supported() ? ".png" : ".bmp";
}

void CImgVisualizer::write(const PlacementDB& db,
                            const PlacementState& state,
                            const std::string& output_path,
                            const PlotOptions& options) const {
#if !PLC_HAS_CIMG
    (void)db;
    (void)state;
    (void)output_path;
    (void)options;
    throw std::runtime_error(
        "CImg visualization is unavailable: install cimg-dev (Ubuntu) or provide CImg.h, then reconfigure CMake");
#else
    if (state.size() < db.design.cells.size()) {
        throw std::runtime_error("Cannot plot: PlacementState does not contain all cells");
    }
    const Region region = chip_region(db, state);
    const double region_width = region.ux - region.lx;
    const double region_height = region.uy - region.ly;
    if (!finite(region_width) || !finite(region_height) || region_width <= 0.0 || region_height <= 0.0) {
        throw std::runtime_error("Cannot plot: invalid chip region dimensions");
    }

    const std::size_t width = std::max(kMinimumImageDimension, options.width);
    const std::size_t margin = std::min(options.margin, width / 4);
    const std::size_t drawable_width = std::max<std::size_t>(1, width - 2 * margin);
    const double desired_height = static_cast<double>(drawable_width) * region_height / region_width +
                                  2.0 * static_cast<double>(margin);
    const std::size_t height = std::max(kMinimumImageDimension,
        std::min(kMaximumImageHeight, static_cast<std::size_t>(std::ceil(desired_height))));
    const std::size_t drawable_height = std::max<std::size_t>(1, height - 2 * margin);
    const double scale = std::min(static_cast<double>(drawable_width) / region_width,
                                  static_cast<double>(drawable_height) / region_height);
    if (!finite(scale) || scale <= 0.0) {
        throw std::runtime_error("Cannot plot: invalid coordinate scale");
    }

    cimg_library::CImg<unsigned char> image(static_cast<int>(width), static_cast<int>(height), 1, 3, 255);
    const Transform transform{region, static_cast<int>(height), static_cast<int>(margin), scale};
    const unsigned char row_color[] = {232, 232, 232};
    const unsigned char core_color[] = {30, 30, 30};

    if (options.draw_rows) {
        for (const auto& row : db.design.rows) {
            if (!finite(row.origin_y) || !finite(row.x_min()) || !finite(row.x_max()) ||
                row.height <= 0.0 || row.x_max() <= row.x_min()) {
                continue;
            }
            const int top = std::min(transform.y(row.origin_y), transform.y(row.origin_y + row.height));
            const int bottom = std::max(transform.y(row.origin_y), transform.y(row.origin_y + row.height));
            image.draw_rectangle(transform.x(row.x_min()), top, transform.x(row.x_max()), bottom,
                                 row_color, 0.75F);
        }
    }

    if (options.draw_nets) {
        draw_nets(&image, db, state, transform);
    }
    for (const auto& cell : db.design.cells) {
        draw_cell(&image, cell, state, transform);
    }

    const auto& core = db.design.core;
    const int left = std::min(transform.x(core.lx), transform.x(core.ux));
    const int right = std::max(transform.x(core.lx), transform.x(core.ux));
    const int top = std::min(transform.y(core.ly), transform.y(core.uy));
    const int bottom = std::max(transform.y(core.ly), transform.y(core.uy));
    image.draw_rectangle(left, top, right, bottom, core_color, 1.0F, ~0U);

    const std::filesystem::path path = normalized_path(output_path);
    try {
        image.save(path.string().c_str());
    } catch (const cimg_library::CImgException& error) {
        throw std::runtime_error("Cannot write plot '" + path.string() + "': " + error.what());
    }
    std::cout << "    Plot Written To: " << path.string() << '\n';
#endif
}

}  // namespace plc

cmake_minimum_required(VERSION 3.16)

project(placement_course_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/utf-8)
endif()

add_executable(placer
    src/algorithms/clustering_placement.cpp
    src/main.cpp
    src/algorithms/global_placement.cpp
    src/algorithms/local_sa_placer.cpp
    src/algorithms/quadratic_analytical_placement.cpp
    src/algorithms/random_placement.cpp
    src/database/placement_db.cpp
    src/database/placement_state.cpp
    src/database/validation_report.cpp
    src/io/summary_writer.cpp
    src/metrics/area_stats.cpp
    src/metrics/hpwl.cpp
    src/parser/aux_parser.cpp
    src/parser/bookshelf_loader.cpp
    src/parser/nets_parser.cpp
    src/parser/nodes_parser.cpp
    src/parser/parser_common.cpp
    src/parser/pl_parser.cpp
    src/parser/scl_parser.cpp
    src/parser/wts_parser.cpp
    src/visualization/cimg_visualizer.cpp
)

target_include_directories(placer PRIVATE include)
target_compile_definitions(placer PRIVATE cimg_display=0)

find_path(CIMG_INCLUDE_DIR NAMES CImg.h PATH_SUFFIXES CImg)
if(CIMG_INCLUDE_DIR)
    message(STATUS "CImg found: ${CIMG_INCLUDE_DIR}")
    target_include_directories(placer PRIVATE "${CIMG_INCLUDE_DIR}")
    target_compile_definitions(placer PRIVATE PLC_HAS_CIMG=1)

    find_package(PNG QUIET)
    if(PNG_FOUND)
        message(STATUS "PNG found: PNG output is enabled")
        target_link_libraries(placer PRIVATE PNG::PNG)
        target_compile_definitions(placer PRIVATE PLC_HAS_PNG=1 cimg_use_png=1)
    else()
        message(WARNING "libpng not found: CImg visualization supports .bmp only")
        target_compile_definitions(placer PRIVATE PLC_HAS_PNG=0)
    endif()
else()
    message(WARNING
        "CImg.h was not found: --plot will report an actionable error. "
        "Install cimg-dev (Ubuntu) or provide CImg.h and reconfigure CMake.")
    target_compile_definitions(placer PRIVATE PLC_HAS_CIMG=0 PLC_HAS_PNG=0)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(clustering_test tests/clustering_test.cpp
        src/algorithms/quadratic_analytical_placement.cpp src/algorithms/local_sa_placer.cpp
        src/algorithms/clustering_placement.cpp src/database/placement_db.cpp
        src/database/placement_state.cpp src/metrics/hpwl.cpp)
    target_include_directories(clustering_test PRIVATE include)
    add_test(NAME clustering_and_hpwl COMMAND clustering_test)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME cli_integration COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_integration.py" "$<TARGET_FILE:placer>")
    endif()
endif()

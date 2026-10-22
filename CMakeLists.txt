cmake_minimum_required(VERSION 3.16)

project(placement_course_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(placer
    src/main.cpp
    src/algorithms/global_placement.cpp
    src/algorithms/local_sa_placer.cpp
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
)

target_include_directories(placer PRIVATE include)

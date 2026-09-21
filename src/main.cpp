#include <iostream>
#include <string>

#include "io/summary_writer.h"
#include "parser/bookshelf_loader.h"

namespace {

void PrintUsage(const char* program) {
    std::cerr << "Usage: " << program << " <benchmark.aux> [--check-only]\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        PrintUsage(argv[0]);
        return 1;
    }

    const std::string aux_path = argv[1];
    plc::BookshelfLoader loader;
    plc::LoadResult result = loader.load_from_aux(aux_path);

    plc::SummaryWriter writer;
    if (!result.report.ok()) {
        writer.write_validation_report(result.report, std::cout);
        return 2;
    }

    writer.write_bookshelf_report(result.db, std::cout);
    return 0;
}

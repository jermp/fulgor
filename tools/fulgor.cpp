#include <iostream>

#include "../include/index.hpp"
#include "../src/index.cpp"
#include "../include/index_types.hpp"

#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"

#include "run_ggcat.cpp"
#include "invert.cpp"
#include "sort_unique.cpp"
#include "permute_unitigs.cpp"
#include "build.cpp"
#include "pseudoalign.cpp"

using namespace fulgor;

int stats(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    index_type index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    index.print_stats();
    return 0;
}

int print_filenames(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    index_type index;
    essentials::logger("loading cc_index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    auto nd = index.num_docs();
    for (size_t i = 0; i != nd; ++i) std::cout << i << '\t' << index.filename(i) << '\n';
    return 0;
}

int help(char* arg0) {
    std::cout << "== Fulgor: a colored compacted de Bruijn graph index =========================="
              << std::endl
              << std::endl;
    std::cout << "Usage: " << arg0 << " <tool> ...\n\n"
              << "Available tools:\n"
              << "  invert          \t invert the reference-to-unitig Cuttlefish's file \n"
              << "  sort-unique     \t deduplicate the color sets \n"
              << "  permute-unitigs \t permute unitigs according to the color they map to \n"
              << "  build           \t build a fulgor index \n"
              << "  pseudoalign     \t pseudoalign reads to references using a fulgor index \n"
              << "  stats           \t print index statistics \n"
              << "  print-filenames \t print all reference filenames " << std::endl;
    return 1;
}

int main(int argc, char** argv) {
    if (argc < 2) return help(argv[0]);
    auto tool = std::string(argv[1]);
    if (tool == "run_ggcat") {
        return run_ggcat(argc - 1, argv + 1);
    } else if (tool == "invert") {
        return invert(argc - 1, argv + 1);
    } else if (tool == "sort-unique") {
        return sort_unique(argc - 1, argv + 1);
    } else if (tool == "permute-unitigs") {
        return permute_unitigs(argc - 1, argv + 1);
    } else if (tool == "build") {
        return build(argc - 1, argv + 1);
    } else if (tool == "pseudoalign") {
        return pseudoalign(argc - 1, argv + 1);
    } else if (tool == "stats") {
        return stats(argc - 1, argv + 1);
    } else if (tool == "print-filenames") {
        return print_filenames(argc - 1, argv + 1);
    }
    std::cout << "Unsupported tool '" << tool << "'.\n" << std::endl;
    return help(argv[0]);
}

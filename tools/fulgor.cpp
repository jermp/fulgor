#include <iostream>
#include <filesystem>

#include "external/sshash/external/gz/zip_stream.hpp"
#include "external/sshash/external/gz/zip_stream.cpp"
#include "external/sshash/src/builder/build.cpp"
#include "external/sshash/src/dictionary.cpp"
#include "external/sshash/src/info.cpp"
#include "external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "external/FQFeeder/include/FastxParser.hpp"
#include "external/FQFeeder/src/FastxParser.cpp"

#include "include/index_types.hpp"
#include "src/index.cpp"
#include "src/color_sets.cpp"

#include "util.cpp"
#include "build.cpp"
#include "permute.cpp"
#include "pseudoalign.cpp"
#include "kmer_conservation.cpp"
#include "kmer_matches.cpp"

int help(char* arg0) {
    std::cout << "== Fulgor: a colored de Bruijn graph index"
              << " (v"
              << essentials::version_number(constants::current_version_number::major,
                                            constants::current_version_number::minor,
                                            constants::current_version_number::patch)
                     .to_string()
              << ')' << " ======================================="  //
              << std::endl
              << std::endl;

    std::cout << "Usage: " << arg0 << " <tool> ...\n\n";

    std::cout << "Construction:\n"
              << "  build              build an index\n"
              << "  color              build a meta- or a diff- or a meta-diff- index\n"
              << "  permute            permute the reference names of an index\n"
              << std::endl;

    std::cout << "Queries:\n"
              << "  pseudoalign        perform pseudoalignment to an index\n"
              << "  kmer-conservation  print color set info for each positive kmer in query\n"
              << "  kmer-matches       print positive kmers per query and number of kmer matches "
                 "per color\n"
              << std::endl;

    std::cout
        << "Debug:\n"
        << "  check              perform an in-depth check to verify that an index was built "
           "correctly\n"
        << "  verify             verify that index works correctly with current library version\n"
        << "  stats              print index statistics\n"
        << "  print-filenames    print all reference filenames\n"
        << "  dump               write unitigs and color sets of an index in text format\n"
        << "  load               build an index from dump output\n"
        << std::endl;

    std::cout << "Other:\n"
              << "  help               print this helper and exit gracefully\n"
              << std::endl;

    return 1;
}

int main(int argc, char** argv) {
    if (argc < 2) return help(argv[0]);

    const auto tool = std::string(argv[1]);

    using ToolFunction = int (*)(int, char**);
    const std::unordered_map<std::string, ToolFunction> tool_map{{
        {"build", build},
        {"pseudoalign", pseudoalign},
        {"kmer-conservation", kmer_conservation},
        {"kmer-matches", kmer_matches},
        {"check", check},
        {"verify", verify},
        {"stats", stats},
        {"print-filenames", print_filenames},
        {"permute", permute},
        {"dump", dump},
        {"load", load},
        {"color", color},
    }};

    if (tool == "help") {
        help(argv[0]);
        return 0;
    }
    if (tool == "load") {
        std::cerr << "Operation temporarily disabled" << std::endl;
        return 1;
    }

    // 3. Look up the tool in the map
    const auto it = tool_map.find(tool);
    if (it != tool_map.end()) {
        // Execute the function dynamically
        return it->second(argc - 1, argv + 1);
    }

    std::cout << "Unsupported tool '" << tool << "'." << std::endl;

    return help(argv[0]);
}

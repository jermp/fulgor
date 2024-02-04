#include <iostream>

#include "../include/index_types.hpp"
#include "../src/index.cpp"
#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"

#include "build.cpp"
#include "pseudoalign.cpp"

using namespace fulgor;

template <typename FulgorIndex>
void print_stats(std::string const& index_filename) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    index.print_stats();
}

template <typename FulgorIndex>
void print_filenames(std::string const& index_filename) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    for (uint64_t i = 0; i != index.num_docs(); ++i) {
        std::cout << i << '\t' << index.filename(i) << '\n';
    }
}

template <typename FulgorIndex>
void dump_colors(std::string const& index_filename, std::string const& output_filename) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    std::ofstream os(output_filename.c_str());
    if (!os.is_open()) throw std::runtime_error("cannot open output file");
    index.dump_colors(os);
    os.close();
}

int stats(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    if (sshash::util::ends_with(index_filename,
                                constants::meta_colored_fulgor_filename_extension)) {
        print_stats<meta_index_type>(index_filename);
    } else if (sshash::util::ends_with(index_filename, constants::fulgor_filename_extension)) {
        print_stats<index_type>(index_filename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int print_filenames(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    if (sshash::util::ends_with(index_filename,
                                constants::meta_colored_fulgor_filename_extension)) {
        print_filenames<meta_index_type>(index_filename);
    } else if (sshash::util::ends_with(index_filename, constants::fulgor_filename_extension)) {
        print_filenames<index_type>(index_filename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int dump_colors(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("output_filename", "The output filename where to write colors.", "-o", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    auto output_filename = parser.get<std::string>("output_filename");
    if (sshash::util::ends_with(index_filename,
                                constants::meta_colored_fulgor_filename_extension)) {
        dump_colors<meta_index_type>(index_filename, output_filename);
    } else if (sshash::util::ends_with(index_filename, constants::fulgor_filename_extension)) {
        dump_colors<index_type>(index_filename, output_filename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int help(char* arg0) {
    std::cout << "== Fulgor: a (meta-) colored compacted de Bruijn graph index "
                 "============================="
              << std::endl
              << std::endl;

    std::cout << "Usage: " << arg0 << " <tool> ...\n\n";

    std::cout << "Tools:\n"
              << "  build              build a Fulgor index\n"
              << "  pseudoalign        pseudoalign reads to references\n"
              << "  stats              print index statistics\n"
              << "  print-filenames    print all reference filenames\n"
              << std::endl;

    std::cout
        << "Advanced tools:\n"
        << "  partition          partition a Fulgor index and build a meta-colored Fulgor index\n"
        << "  dump-colors        write colors to an output file in text format" << std::endl;

    return 1;
}

int main(int argc, char** argv) {
    if (argc < 2) return help(argv[0]);

    auto tool = std::string(argv[1]);

    /* basic tools */
    if (tool == "build") {
        return build(argc - 1, argv + 1);
    } else if (tool == "pseudoalign") {
        return pseudoalign(argc - 1, argv + 1);
    } else if (tool == "stats") {
        return stats(argc - 1, argv + 1);
    } else if (tool == "print-filenames") {
        return print_filenames(argc - 1, argv + 1);
    }

    /* advanced tools */
    else if (tool == "partition") {
        return partition(argc - 1, argv + 1);
    } else if (tool == "dump-colors") {
        return dump_colors(argc - 1, argv + 1);
    }

    std::cout << "Unsupported tool '" << tool << "'.\n" << std::endl;

    return help(argv[0]);
}

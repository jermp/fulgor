#include <iostream>
#include <numeric>
#include <vector>

#include "../include/index.hpp"
#include "../src/index.cpp"
#include "../include/index_types.hpp"

#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"

using namespace fulgor;

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("cc_index_filename", "The fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    auto cc_index_filename = parser.get<std::string>("cc_index_filename");

    index_type index;
    essentials::logger("loading index from disk...");
    essentials::load(index, cc_index_filename.c_str());
    essentials::logger("DONE");
    index.print_stats();

    return 0;
}
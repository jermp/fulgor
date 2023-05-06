#include <iostream>
#include <numeric>
#include <vector>

#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/index.hpp"
#include "../include/index_types.hpp"
#include "../src/index.cpp"

using namespace fulgor;

int main(int argc, char** argv) {
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
    for (size_t i = 0; i < nd; ++i) { std::cout << i << '\t' << index.filename(i) << '\n'; }

    return 0;
}

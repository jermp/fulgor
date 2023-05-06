#include <iostream>

#include "../include/inverted_index_builder.hpp"
#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"

using namespace fulgor;

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("file_base_name", "Cuttlefish input file_base_name.", "-i", true);
    parser.add("RAM",
               "RAM limit in GiB. Default value is " +
                   std::to_string(constants::default_ram_limit_in_GiB) + ".",
               "-g", false);
    parser.add("num_refs", "The number of references to index.", "-n", false);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("verbose", "Verbose output during construction.", "--verbose", false, true);

    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    build_configuration build_config;
    build_config.file_base_name = parser.get<std::string>("file_base_name");
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    if (parser.parsed("RAM")) build_config.ram_limit_in_GiB = parser.get<double>("RAM");
    if (parser.parsed("num_refs")) build_config.num_refs = parser.get<uint32_t>("num_refs");
    if (parser.parsed("verbose")) {
        assert(parser.get<bool>("verbose"));
        build_config.verbose = true;
    }

    inverted_index_builder ii_builder(build_config);
    ii_builder.build();

    return 0;
}
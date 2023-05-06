#include "../include/list_sorter.hpp"

using namespace fulgor;

int sort_unique(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("file_base_name", "Cuttlefish input file_base_name.", "-i", true);
    parser.add("RAM",
               "RAM limit in GiB. Default value is " +
                   std::to_string(constants::default_ram_limit_in_GiB) + ".",
               "-g", false);
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
    build_config.verbose = parser.get<bool>("verbose");

    mm::file_source<uint64_t> mm_index_file(build_config.file_base_name + ".inv_idx",
                                            mm::advice::sequential);
    inverted_index::iterator<> it(mm_index_file.data(), mm_index_file.size());
    uint32_t num_docs = it.num_docs();
    std::cout << "num_docs: " << num_docs << std::endl;

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();

    list_sorter sorter(build_config, num_docs);
    uint64_t i = 0;
    while (it.has_next()) {
        auto list = it.list();
        sorter.add(list.seg_id(), list.begin(), list.end(), list.size());
        it.advance_to_next_list();
        ++i;
    }
    mm_index_file.close();
    sorter.finalize();

    timer.stop();
    std::cout << "** sorting took " << timer.elapsed() << " seconds / " << timer.elapsed() / 60
              << " minutes" << std::endl;

    return 0;
}
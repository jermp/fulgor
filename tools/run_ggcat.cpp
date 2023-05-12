#include "../include/GGCAT.hpp"

using namespace fulgor;

int run_ggcat(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("filenames_list", "Filenames list.", "-l", true);
    parser.add("output_basename", "Output basename.", "-o", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(sshash::constants::max_k) + ").",
               "-k", true);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("RAM",
               "RAM limit in GiB. Default value is " +
                   std::to_string(constants::default_ram_limit_in_GiB) + ".",
               "-g", false);
    parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    std::string tmp_dirname = constants::default_tmp_dirname;
    if (parser.parsed("tmp_dirname")) {
        tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(tmp_dirname);
    }

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();

    std::vector<std::string> filenames;
    {
        std::ifstream in(parser.get<std::string>("filenames_list"));
        if (!in.is_open()) throw std::runtime_error("error in opening file");
        std::string filename;
        while (in >> filename) filenames.push_back(filename);
        std::cout << "about to process " << filenames.size() << " files..." << std::endl;
    }

    uint64_t mem_gigas = parser.get<uint64_t>("RAM");
    uint64_t k = parser.get<uint64_t>("k");
    uint64_t num_threads = parser.get<uint64_t>("num_threads");
    std::string output_basename = parser.get<std::string>("output_basename");
    GGCAT builder(filenames, mem_gigas, k, num_threads, tmp_dirname, output_basename);

    timer.stop();
    std::cout << "** running GGCAT took " << timer.elapsed() << " seconds / "
              << timer.elapsed() / 60 << " minutes" << std::endl;

    essentials::logger("DONE");

    return 0;
}
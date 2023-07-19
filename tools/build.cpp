using namespace fulgor;

int build(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("filenames_list", "Filenames list.", "-l", true);
    parser.add("file_base_name", "File basename.", "-o", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(sshash::constants::max_k) + ").",
               "-k", true);
    parser.add("m", "Minimizer length (must be < k).", "-m", true);
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
    parser.add("verbose", "Verbose output during construction.", "--verbose", false, true);
    parser.add("check", "Check correctness after index construction (it might take some time).",
               "--check", false, true);

    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    build_configuration build_config;
    build_config.file_base_name = parser.get<std::string>("file_base_name");
    auto k = parser.get<uint64_t>("k");
    auto m = parser.get<uint64_t>("m");
    build_config.k = k;
    build_config.m = m;
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    build_config.verbose = parser.get<bool>("verbose");
    build_config.check = parser.get<bool>("check");
    build_config.filenames_list = parser.get<std::string>("filenames_list");
    if (parser.get<uint64_t>("RAM")) {
        build_config.ram_limit_in_GiB = parser.get<uint64_t>("RAM");
    }
    if (parser.parsed("num_threads")) {
        build_config.num_threads = parser.get<uint64_t>("num_threads");
    }

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();

    index_type index;
    typename index_type::builder builder(build_config);
    builder.build(index);
    index.print_stats();

    timer.stop();
    essentials::logger("DONE");
    std::cout << "** building the index took " << timer.elapsed() << " seconds / "
              << timer.elapsed() / 60 << " minutes" << std::endl;

    std::string output_filename =
        build_config.file_base_name + "." + index_type::color_classes_type::type() + ".index";
    essentials::logger("saving index to disk...");
    essentials::save(index, output_filename.c_str());
    essentials::logger("DONE");

    return 0;
}
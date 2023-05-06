using namespace fulgor;

int build(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("file_base_name", "Cuttlefish input file_base_name.", "-i", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(sshash::constants::max_k) + ").",
               "-k", true);
    parser.add("m", "Minimizer length (must be < k).", "-m", true);
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
    auto k = parser.get<uint64_t>("k");
    auto m = parser.get<uint64_t>("m");
    build_config.k = k;
    build_config.m = m;
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    build_config.verbose = parser.get<bool>("verbose");

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();
    index_type index;
    index.build(build_config);
    index.print_stats();
    timer.stop();
    std::cout << "** building the index took " << timer.elapsed() << " seconds / "
              << timer.elapsed() / 60 << " minutes" << std::endl;

    std::string output_filename =
        build_config.file_base_name + "." + index_type::color_classes_type::type() + ".index";
    essentials::logger("saving index to disk...");
    essentials::save(index, output_filename.c_str());
    essentials::logger("DONE");

    return 0;
}
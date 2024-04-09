using namespace fulgor;

int permute(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename",
               "The Fulgor index filename from which we permute the reference names.", "-i", true);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("output_filename", "Output file where to save the permuted filenames.", "-o", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    build_configuration build_config;
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }

    auto index_filename = parser.get<std::string>("index_filename");

    if (!sshash::util::ends_with(index_filename, "." + constants::fulgor_filename_extension)) {
        std::cerr << "Error: the file to partition must have extension \"."
                  << constants::fulgor_filename_extension
                  << "\". Have you first built a Fulgor index with the tool \"build\"?"
                  << std::endl;
        return 1;
    }

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();

    index_type index;
    essentials::logger("step 1. loading index to be partitioned...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");

    permuter p(build_config);
    p.permute(index);
    auto const& filenames = p.filenames();

    std::ofstream out(parser.get<std::string>("output_filename").c_str());
    if (!out.is_open()) {
        std::cerr << "cannot open output filename" << std::endl;
        return 1;
    }
    for (auto const& fn : filenames) out << fn << '\n';
    out.close();

    timer.stop();
    essentials::logger("DONE");
    std::cout << "** permuting the reference names took " << timer.elapsed() << " seconds / "
              << timer.elapsed() / 60 << " minutes" << std::endl;

    return 0;
}
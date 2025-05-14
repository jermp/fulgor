using namespace fulgor;

void meta_color(build_configuration const& build_config) {
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();

    meta_index_type index;
    typename meta_index_type::meta_builder builder(build_config);
    builder.build(index);
    index.print_stats();

    timer.stop();
    essentials::logger("DONE");
    std::cout << "** building the index took " << timer.elapsed() << " seconds / "
              << timer.elapsed() / 60 << " minutes" << std::endl;

    std::string output_filename = build_config.index_filename_to_partition.substr(
                                      0, build_config.index_filename_to_partition.length() -
                                             constants::fulgor_filename_extension.length() - 1) +
                                  "." + constants::meta_colored_fulgor_filename_extension;
    essentials::logger("saving index to disk...");
    essentials::save(index, output_filename.c_str());
    essentials::logger("DONE");
}

void diff_color(build_configuration const& build_config) {
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();

    differential_index_type index;
    typename differential_index_type::differential_builder builder(build_config);
    builder.build(index);
    index.print_stats();

    timer.stop();
    essentials::logger("DONE");
    std::cout << "** building the index took " << timer.elapsed() << " seconds / "
              << timer.elapsed() / 60 << " minutes" << std::endl;

    std::string output_filename = build_config.index_filename_to_partition.substr(
                                      0, build_config.index_filename_to_partition.length() -
                                             constants::fulgor_filename_extension.length() - 1) +
                                  "." + constants::diff_colored_fulgor_filename_extension;
    essentials::logger("saving index to disk...");
    essentials::save(index, output_filename.c_str());
    essentials::logger("DONE");
}

void meta_diff_color(build_configuration const& build_config)  //
{
    { /* first build a meta-colored Fulgor index */
        meta_color(build_config);
    }

    build_configuration meta_diff_build_config = build_config;
    meta_diff_build_config.index_filename_to_partition =
        build_config.index_filename_to_partition.substr(
            0, build_config.index_filename_to_partition.length() -
                   constants::fulgor_filename_extension.length() - 1) +
        "." + constants::meta_colored_fulgor_filename_extension;

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();

    meta_differential_index_type index;
    typename meta_differential_index_type::meta_differential_builder builder(
        meta_diff_build_config);
    builder.build(index);
    index.print_stats();

    timer.stop();
    essentials::logger("DONE");
    std::cout << "** building the index took " << timer.elapsed() << " seconds / "
              << timer.elapsed() / 60 << " minutes" << std::endl;

    std::string output_filename =
        meta_diff_build_config.index_filename_to_partition.substr(
            0, meta_diff_build_config.index_filename_to_partition.length() -
                   constants::meta_colored_fulgor_filename_extension.length() - 1) +
        "." + constants::meta_diff_colored_fulgor_filename_extension;
    essentials::logger("saving index to disk...");
    essentials::save(index, output_filename.c_str());
    essentials::logger("DONE");
}

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
    parser.add("force", "Re-build the index even when an index with the same name is found.",
               "--force", false, true);
    parser.add("meta", "Build a meta-colored index.", "--meta", false, true);
    parser.add("diff", "Build a differential-colored index.", "--diff", false, true);

    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    build_configuration build_config;
    build_config.file_base_name = parser.get<std::string>("file_base_name");
    std::string output_filename =
        build_config.file_base_name + "." + constants::fulgor_filename_extension;
    build_config.index_filename_to_partition = output_filename;
    bool force = parser.get<bool>("force");
    build_config.meta_colored = parser.get<bool>("meta");
    build_config.diff_colored = parser.get<bool>("diff");

    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    if (parser.parsed("num_threads")) {
        build_config.num_threads = parser.get<uint64_t>("num_threads");
    }

    if (std::filesystem::exists(output_filename) and !force) {
        std::cerr << "An index with the name '" << output_filename << "' alreay exists."
                  << std::endl;
        std::cerr << "Use option '--force' to re-build the index." << std::endl;
        if (build_config.meta_colored and build_config.diff_colored) {
            std::cerr << "Consider using: \"./fulgor color -i " << output_filename << " -d "
                      << build_config.tmp_dirname << " -t "
                      << std::to_string(build_config.num_threads) << " --diff --meta\""
                      << std::endl;
        } else if (build_config.meta_colored) {
            std::cerr << "Consider using: \"./fulgor color -i " << output_filename << " -d "
                      << build_config.tmp_dirname << " -t "
                      << std::to_string(build_config.num_threads) << " --meta\"" << std::endl;
        } else if (build_config.diff_colored) {
            std::cerr << "Consider using: \"./fulgor color -i " << output_filename << " -d "
                      << build_config.tmp_dirname << " -t "
                      << std::to_string(build_config.num_threads) << " --diff\"" << std::endl;
        }
        return 1;
    }

    auto k = parser.get<uint64_t>("k");
    auto m = parser.get<uint64_t>("m");
    build_config.k = k;
    build_config.m = m;
    build_config.verbose = parser.get<bool>("verbose");
    build_config.check = parser.get<bool>("check");
    build_config.filenames_list = parser.get<std::string>("filenames_list");
    if (parser.get<uint64_t>("RAM")) {
        build_config.ram_limit_in_GiB = parser.get<uint64_t>("RAM");
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

    essentials::logger("saving index to disk...");
    essentials::save(index, output_filename.c_str());
    essentials::logger("DONE");

    if (build_config.meta_colored and build_config.diff_colored) {
        meta_diff_color(build_config);
    } else if (build_config.meta_colored) {
        meta_color(build_config);
    } else if (build_config.diff_colored) {
        diff_color(build_config);
    }

    return 0;
}

int color(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename to partition.", "-i", true);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    parser.add("check", "Check correctness after index construction (it might take some time).",
               "--check", false, true);
    parser.add("meta", "Build a meta-colored index.", "--meta", false, true);
    parser.add("diff", "Build a differential-colored index.", "--diff", false, true);

    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    build_configuration build_config;
    build_config.index_filename_to_partition = parser.get<std::string>("index_filename");
    if (!sshash::util::ends_with(build_config.index_filename_to_partition,
                                 "." + constants::fulgor_filename_extension)) {
        std::cerr << "Error: the file to partition must have extension \"."
                  << constants::fulgor_filename_extension
                  << "\". Have you first built a Fulgor index with the tool \"build\"?"
                  << std::endl;
        return 1;
    }

    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    if (parser.parsed("num_threads")) {
        build_config.num_threads = parser.get<uint64_t>("num_threads");
    }
    build_config.check = parser.get<bool>("check");
    build_config.meta_colored = parser.get<bool>("meta");
    build_config.diff_colored = parser.get<bool>("diff");

    if (build_config.meta_colored and build_config.diff_colored) {
        meta_diff_color(build_config);
    } else if (build_config.meta_colored) {
        meta_color(build_config);
    } else if (build_config.diff_colored) {
        diff_color(build_config);
    } else {
        std::cerr << "Either \"--meta\" or \"--diff\" should be specified." << std::endl;
        return 1;
    }

    return 0;
}

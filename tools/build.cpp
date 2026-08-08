using namespace fulgor;

[[deprecated("Fix or remove")]] void meta_color(build_configuration& build_config,
                                                const bool force)  //
{
    build_config.output_filename.replace_extension(constants::mfur_filename_extension);

    if (std::filesystem::exists(build_config.output_filename)) {
        std::cerr << "An index with the name '" << build_config.output_filename
                  << "' already exists." << std::endl;
        if (force) {
            std::cerr << "Option '--force' specified: re-building the index." << std::endl;
        } else {
            std::cerr << "Use option '--force' to re-build the index." << std::endl;
            return;
        }
    }

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();

    mfur_index_t index;
    // mfur_index_t::meta_builder builder(build_config);
    // builder.build(index);

    timer.stop();
    essentials::logger("BUILDING DONE");
    std::cout << "** building the index took " << timer.elapsed() << " seconds / "
              << timer.elapsed() / 60 << " minutes" << std::endl;

    essentials::logger("saving index to disk...");
    // essentials::save(index, output_filename.c_str());
    essentials::logger("DONE");
    essentials::load(index, build_config.output_filename.c_str());

    if (build_config.verbose) index.print_stats();
    // if (build_config.check) builder.check(index);
}

void diff_color(build_configuration& build_config, const bool force)  //
{
    build_config.output_filename.replace_extension(constants::dfur_filename_extension);

    if (std::filesystem::exists(build_config.output_filename)) {
        std::cerr << "An index with the name '" << build_config.output_filename
                  << "' already exists." << std::endl;
        if (force) {
            std::cerr << "Option '--force' specified: re-building the index." << std::endl;
        } else {
            std::cerr << "Use option '--force' to re-build the index." << std::endl;
            return;
        }
    }

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();

    dfur_index_t index;
    dfur_index_t::differential_builder builder(build_config);
    builder.build(index);

    timer.stop();
    essentials::logger("BUILDING DONE");
    std::cout << "** building the index took " << timer.elapsed() << " seconds / "
              << timer.elapsed() / 60 << " minutes" << std::endl;

    essentials::logger("saving index to disk...");
    essentials::save(index, build_config.output_filename.c_str());
    essentials::logger("DONE");

    if (build_config.verbose) {
        index.print_stats();
    }
    if (build_config.check) {
        builder.check(index);
    }
}

void meta_diff_color(build_configuration& build_config, const bool force)  //
{
    build_config.output_filename.replace_extension(constants::mdfur_filename_extension);

    if (std::filesystem::exists(build_config.output_filename)) {
        std::cerr << "An index with the name '" << build_config.output_filename
                  << "' already exists." << std::endl;
        if (force) {
            std::cerr << "Option '--force' specified: re-building the index." << std::endl;
        } else {
            std::cerr << "Use option '--force' to re-build the index." << std::endl;
            return;
        }
    }

    std::filesystem::path meta_filename = build_config.index_filename_to_partition;
    meta_filename.replace_extension(constants::mfur_filename_extension);

    /* first build a meta-colored Fulgor index */
    if (!std::filesystem::exists(meta_filename)) {
        meta_color(build_config, force);
    } else {
        std::cout << ".mfur file found, skipping meta partitioning" << std::endl;
    }

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> build_timer;
    build_timer.start();

    build_configuration meta_diff_build_config = build_config;
    meta_diff_build_config.index_filename_to_partition = meta_filename;

    mdfur_index_t index;
    mdfur_index_t::meta_differential_builder builder(meta_diff_build_config);
    builder.build(index);
    build_timer.stop();

    essentials::logger("BUILD DONE");
    std::cout << "** building the index took " << build_timer.elapsed() << " seconds / "
              << build_timer.elapsed() / 60 << " minutes" << std::endl;

    essentials::logger("saving index to disk...");
    essentials::save(index, build_config.output_filename.c_str());
    essentials::logger("DONE");

    if (build_config.verbose) {
        index.print_stats();
    }
    if (build_config.check) {
        builder.check(index);
    }
}

int build(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("filenames_list", "Filenames list.", "-l", true);
    parser.add("file_base_name", "File basename.", "-o", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(kmer_type::max_k) + ").", "-k",
               true);
    parser.add("m", "Minimizer length (must be < k).", "-m", true);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory (default is directory '" +
            constants::default_tmp_dirname + "').",
        "-d", false);
    parser.add("RAM",
               "RAM limit in GiB (default is " +
                   std::to_string(constants::default_ram_limit_in_GiB) + ").",
               "-g", false);
    parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    parser.add("verbose", "Verbose output during construction.", "--verbose", false, true);
    parser.add("check", "Check correctness after index construction (it might take some time).",
               "--check", false, true);
    parser.add("force", "Re-build the index even when an index with the same name is found.",
               "--force", false, true);
    parser.add("meta", "Build a meta-colored index.", "--meta", false, true);
    parser.add("stats", "Prints index stats after construction.", "--stats", false, true);

    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    build_configuration build_config;
    build_config.output_filename = parser.get<std::string>("file_base_name");
    build_config.output_filename.replace_extension(constants::hfur_filename_extension);
    build_config.index_filename_to_partition = build_config.output_filename;

    const bool force = parser.get<bool>("force");
    build_config.meta_colored = parser.get<bool>("meta");

    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    if (parser.parsed("num_threads")) {
        build_config.num_threads = parser.get<uint64_t>("num_threads");
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

    std::variant<hfur_index_t, mfur_index_t> index;
    if (build_config.meta_colored) {
        build_config.output_filename.replace_extension(constants::mfur_filename_extension);
        index = mfur_index_t();
    } else {
        index = hfur_index_t();
    }

    if (std::filesystem::exists(build_config.output_filename)) {
        std::cerr << "An index with the name '" << build_config.output_filename
                  << "' already exists." << std::endl;
        if (force) {
            std::cerr << "Option '--force' specified: re-building the index." << std::endl;
        } else {
            std::cerr << "Use option '--force' to re-build the index." << std::endl;
            std::string color_flag = "";
            if (build_config.meta_colored) {
                color_flag += "--meta ";
            }

            std::cerr << "Consider using: \"./fulgor color -i " << build_config.output_filename
                      << " -d " << build_config.tmp_dirname << " -t "
                      << std::to_string(build_config.num_threads) << " " << color_flag << "\""
                      << std::endl;
            return 1;
        }
    }

    std::visit(
        [&]<typename Index>(Index index_) {
            {
                util::timed_phase timer("Building the index");

                typename Index::builder b(build_config);
                b.build();

                essentials::logger("Index stored at " + build_config.output_filename.string());
            }

            if (parser.get<bool>("stats")) {
                essentials::mmap(index_, build_config.output_filename.c_str());
                index_.print_stats();
            }

            if (build_config.verbose) {
                const auto rpt_head = "---------- Build Report ----------\n";
                const auto rpt_ln0 =
                    std::format("Index built at: {}\n",
                                std::filesystem::absolute(build_config.output_filename).string());
                const auto rpt_ln1 = std::format(
                    "Check correctness by using the tool \"check -i {} -q 0.01 -p 0.01 "
                    "--verbose\"\n",
                    build_config.output_filename.string());
                const auto rpt_foot = "----------------------------------\n";

                std::cout << rpt_head + rpt_ln0 + rpt_ln1 + rpt_foot << std::endl;

                util::timed_phase::print_breakdown();
            }
        },
        index);

    // if (build_config.check) builder.check(index);

    // if (build_config.meta_colored and build_config.diff_colored) {
    //     meta_diff_color(build_config, force);
    // } else if (build_config.meta_colored) {
    //     meta_color(build_config, force);
    // } else if (build_config.diff_colored) {
    //     diff_color(build_config, force);
    // }

    return 0;
}

int color(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename to partition.", "-i", true);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory (default is directory '" +
            constants::default_tmp_dirname + "').",
        "-d", false);
    parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    parser.add("RAM",
               "RAM limit in GiB (default is " +
                   std::to_string(constants::default_ram_limit_in_GiB) + ").",
               "-g", false);
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
    build_config.index_filename_to_partition = parser.get<std::string>("index_filename");
    build_config.output_filename = build_config.index_filename_to_partition;

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
    build_config.verbose = parser.get<bool>("verbose");
    const bool force = parser.get<bool>("force");
    if (parser.get<uint64_t>("RAM")) {
        build_config.ram_limit_in_GiB = parser.get<uint64_t>("RAM");
    }

    if (build_config.meta_colored && build_config.diff_colored &&
        build_config.index_filename_to_partition.extension() !=
            constants::hfur_filename_extension &&
        build_config.index_filename_to_partition.extension() !=
            constants::mfur_filename_extension) {
        const auto mess = std::format(
            "Error: the file to partition must have extension \"{}\" or \"{}\". Have "
            "you first built a Fulgor index with the tool \"build\"?",
            constants::hfur_filename_extension, constants::mfur_filename_extension);
        std::cerr << mess << std::endl;
        return 1;
    }
    if (build_config.index_filename_to_partition.extension() !=
        constants::hfur_filename_extension) {
        std::cerr << "Error: the file to partition must have extension \"."
                  << constants::hfur_filename_extension
                  << "\". Have you first built a Fulgor index with the tool \"build\"?"
                  << std::endl;
        return 1;
    }

    if (build_config.meta_colored and build_config.diff_colored) {
        meta_diff_color(build_config, force);
    } else if (build_config.meta_colored) {
        meta_color(build_config, force);
    } else if (build_config.diff_colored) {
        diff_color(build_config, force);
    } else {
        std::cerr << "Either \"--meta\" or \"--diff\" should be specified." << std::endl;
        return 1;
    }

    return 0;
}

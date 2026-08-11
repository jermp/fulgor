#include <variant>

using namespace fulgor;

bool is_meta(std::filesystem::path const& index_filename) {
    return index_filename.extension() == constants::mfur_filename_extension;
}

bool is_meta_diff(std::filesystem::path const& index_filename) {
    return index_filename.extension() == constants::mdfur_filename_extension;
}

bool is_diff(std::filesystem::path const& index_filename) {
    return index_filename.extension() == constants::dfur_filename_extension;
}

bool is_hybrid(std::filesystem::path const& index_filename) {
    return index_filename.extension() == constants::hfur_filename_extension;
}

template <typename FulgorIndex>
void verify(std::string const& index_filename) {
    FulgorIndex index;
    essentials::version_number vnum(constants::current_version_number::major,   //
                                    constants::current_version_number::minor,   //
                                    constants::current_version_number::patch);  //
    essentials::loader l(index_filename.c_str());
    l.visit(vnum);
    std::cout << "read version number = " << vnum.to_string() << std::endl;
    util::check_version_number(vnum);
    essentials::logger("OK: Fulgor index is compatible with current library version.");
}

template <typename FulgorIndex>
void print_stats(std::string const& index_filename) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::mmap(index, index_filename.c_str());
    essentials::logger("DONE");
    index.print_stats();
}

template <typename FulgorIndex>
void print_filenames(std::string const& index_filename) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::mmap(index, index_filename.c_str());
    essentials::logger("DONE");
    for (uint64_t i = 0; i != index.num_colors(); ++i) {
        std::cout << i << '\t' << index.filename(i) << '\n';
    }
}

template <typename FulgorIndex>
void dump(std::string const& index_filename, build_configuration const& build_config) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::mmap(index, index_filename.c_str());
    essentials::logger("DONE");
    index.dump(build_config);
}

template <typename BaseIndex, typename TargetIndex>
int check(BaseIndex base, TargetIndex target, uint64_t num_threads, bool verbose) {
    if (target.num_colors() != base.num_colors()) {
        std::cerr << "Number of colors mismatch" << target.num_colors() << std::endl;
        return 1;
    }
    const uint64_t num_colors = base.num_colors();

    if (target.num_color_sets() != base.num_color_sets()) {
        std::cout << "Number of color sets mismatch" << std::endl;
        return 1;
    }
    const uint64_t num_color_sets = base.num_color_sets();

    if (target.get_k2u().num_strings() != base.get_k2u().num_strings()) {
        std::cout << "Number of contigs mismatch" << std::endl;
        return 1;
    }
    const uint64_t num_unitigs = base.get_k2u().num_strings();

    if (target.get_k2u().num_kmers() != base.get_k2u().num_kmers()) {
        std::cout << "Number of kmers mismatch" << std::endl;
        return 1;
    }
    const uint64_t num_kmers = base.get_k2u().num_kmers();

    std::vector<uint32_t> base_to_target(num_colors);
    {
        std::vector<uint32_t> base_permutation(num_colors), target_permutation(num_colors);
        std::iota(base_permutation.begin(), base_permutation.end(), 0);
        std::iota(target_permutation.begin(), target_permutation.end(), 0);
        std::sort(base_permutation.begin(), base_permutation.end(), [&](uint32_t a, uint32_t b) {
            return base.get_filenames()[a] < base.get_filenames()[b];
        });
        std::sort(target_permutation.begin(), target_permutation.end(),
                  [&](uint32_t a, uint32_t b) {
                      return target.get_filenames()[a] < target.get_filenames()[b];
                  });

        for (uint64_t i = 0; i != num_colors; ++i) {
            base_to_target[base_permutation[i]] = target_permutation[i];
        }
    }

    std::vector<bool> checked_color_sets(num_color_sets, false);
    std::atomic<uint64_t> num_checked_unitigs(0), num_checked_color_sets(0), num_checked_kmers(0);

    if (verbose) {
        std::cout << "Unitigs checked:     0/" << num_unitigs << std::endl;
        std::cout << "Color sets checked:  0/" << num_color_sets << std::endl;
        std::cout << "Kmers checked:       0/" << num_kmers << std::endl;
    }

    bool errors_found = false;
    auto exe = [&](uint64_t unitig_begin, uint64_t unitig_end) {
        for (uint64_t unitig_id = unitig_begin; unitig_id != unitig_end; ++unitig_id) {
            if (verbose && ++num_checked_unitigs % 1000 == 0) {
                std::cout << "\033[3A";  // Move up 3 lines

                std::cout << "\033[2K";  // Clear line
                std::cout << "Unitigs checked:     " << num_checked_unitigs << "/" << num_unitigs
                          << std::endl;

                std::cout << "\033[2K";  // Clear line
                std::cout << "Color sets checked:  " << num_checked_color_sets << "/"
                          << num_color_sets << std::endl;

                std::cout << "\033[2K";  // Clear line
                std::cout << "Kmers checked:       " << num_checked_kmers << "/" << num_kmers
                          << std::endl;
            }

            auto it = target.get_k2u().at_string_id(unitig_id);
            auto [_, kmer] = it.next();
            const uint64_t base_string_id = base.get_k2u().lookup(kmer).string_id;
            const uint64_t target_string_id = target.get_k2u().lookup(kmer).string_id;
            ++num_checked_kmers;

            while (it.has_next()) {
                ++num_checked_kmers;
                auto [_, kmer] = it.next();
                const uint64_t curr_target_string_id = target.get_k2u().lookup(kmer).string_id;
                const uint64_t curr_base_string_id = base.get_k2u().lookup(kmer).string_id;
                if (target_string_id != curr_target_string_id) {  // should never happen
                    errors_found = true;
                    std::cerr << "\033[1;31m"
                              << "expected unitig " << target_string_id << " but found "
                              << curr_target_string_id << "\033[0m" << std::endl;
                }
                if (base_string_id != curr_base_string_id) {
                    errors_found = true;
                    std::cerr << "\033[1;31m"
                              << "expected unitig " << base_string_id << " but found "
                              << curr_base_string_id << "\033[0m" << std::endl;
                }
            }

            uint64_t base_color_set_id = base.u2c(base_string_id);
            uint64_t target_color_set_id = target.u2c(target_string_id);

            if (checked_color_sets[target_color_set_id]) continue;
            checked_color_sets[target_color_set_id] = true;
            ++num_checked_color_sets;

            auto base_it = base.color_set(base_color_set_id);
            auto target_it = target.color_set(target_color_set_id);

            if (target_it.size() != base_it.size()) {
                errors_found = true;
                std::cerr << "\033[1;31m"
                          << "Error while checking color set " << target_color_set_id
                          << ", different sizes: expected " << base_it.size() << " but got "
                          << target_it.size() << "\033[0m" << std::endl;
                continue;
            }

            std::vector<uint32_t> permuted_base;
            for (uint64_t j = 0; j < base_it.size(); ++j, ++base_it) {
                permuted_base.push_back(base_to_target[*base_it]);
            }
            std::sort(permuted_base.begin(), permuted_base.end());
            auto pbase_it = permuted_base.begin();

            for (uint64_t j = 0; j < target_it.size(); ++j, ++pbase_it, ++target_it) {
                auto base_val = *pbase_it;
                auto target_val = *target_it;
                if (base_val != target_val) {
                    errors_found = true;
                    std::cerr << "\033[1;31m"
                              << "Error while checking color set " << target_color_set_id
                              << ", mismatch at position " << j << ": expected "
                              << base_to_target[base_val] << " but got " << target_val << "\033[0m"
                              << std::endl;
                }
            }
        }
    };

    kmeans::thread_pool threads(num_threads);
    const uint64_t load_per_thread = num_unitigs / (num_threads << 10);
    uint64_t start = 0, end = load_per_thread;
    while (end < num_unitigs) {
        threads.enqueue([&, start, end] { exe(start, std::min(end, num_unitigs)); });
        start = end;
        end = std::min(end + load_per_thread, num_unitigs);
    }
    threads.enqueue([&, start, num_unitigs] { exe(start, num_unitigs); });  // last one
    threads.wait();

    if (verbose) {
        std::cout << "\033[3A";  // Move up 3 lines

        std::cout << "\033[2K";  // Clear line
        std::cout << "Unitigs checked:     " << num_checked_unitigs << "/" << num_unitigs
                  << std::endl;

        std::cout << "\033[2K";  // Clear line
        std::cout << "Color sets checked:  " << num_checked_color_sets << "/" << num_color_sets
                  << std::endl;

        std::cout << "\033[2K";  // Clear line
        std::cout << "Kmers checked:       " << num_checked_kmers << "/" << num_kmers << std::endl;
    }

    return errors_found ? 1 : 0;
}

int verify(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    if (is_meta(index_filename)) {
        verify<mfur_index_t>(index_filename);
    } else if (is_meta_diff(index_filename)) {
        verify<mdfur_index_t>(index_filename);
    } else if (is_diff(index_filename)) {
        verify<dfur_index_t>(index_filename);
    } else if (is_hybrid(index_filename)) {
        verify<hfur_index_t>(index_filename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int stats(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    if (is_meta(index_filename)) {
        print_stats<mfur_index_t>(index_filename);
    } else if (is_meta_diff(index_filename)) {
        print_stats<mdfur_index_t>(index_filename);
    } else if (is_diff(index_filename)) {
        print_stats<dfur_index_t>(index_filename);
    } else if (is_hybrid(index_filename)) {
        print_stats<hfur_index_t>(index_filename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int print_filenames(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    if (is_meta_diff(index_filename)) {
        print_filenames<mdfur_index_t>(index_filename);
    } else if (is_meta(index_filename)) {
        print_filenames<mfur_index_t>(index_filename);
    } else if (is_diff(index_filename)) {
        print_filenames<dfur_index_t>(index_filename);
    } else if (is_hybrid(index_filename)) {
        print_filenames<hfur_index_t>(index_filename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int dump(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("output_basename",
               "Output basename. If the supplied basename is F, the output will consist in four "
               "files: F.unitigs.fa, F.color_sets.txt, F.filenames.txt, and F.metadata.txt. (If "
               "this is not supplied, the basename of the index is used instead.)",
               "-o", false);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    auto index_filename = parser.get<std::string>("index_filename");
    std::string output_basename("");
    if (parser.parsed("output_basename")) {
        output_basename = parser.get<std::string>("output_basename");
        assert(output_basename.length() != 0);
    }

    build_configuration build_config;
    build_config.output_filename = output_basename;

    if (is_meta_diff(index_filename)) {
        build_config.output_filename.replace_extension(constants::mdfur_filename_extension);
        dump<mdfur_index_t>(index_filename, build_config);
    } else if (is_meta(index_filename)) {
        build_config.output_filename.replace_extension(constants::mfur_filename_extension);
        dump<mfur_index_t>(index_filename, build_config);
    } else if (is_diff(index_filename)) {
        build_config.output_filename.replace_extension(constants::dfur_filename_extension);
        dump<dfur_index_t>(index_filename, build_config);
    } else if (is_hybrid(index_filename)) {
        build_config.output_filename.replace_extension(constants::hfur_filename_extension);
        dump<hfur_index_t>(index_filename, build_config);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int load(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add(
        "input_basename",
        "Input basename. If the supplied basename is F, the program attempts to read four files: "
        "F.unitigs.fa, F.color_sets.txt, F.filenames.txt, and F.metadata.txt. (These files are the "
        "output of the `dump` tool.)",
        "-i", true);
    parser.add("output_basename",
               "Output basename. If not provided, the input basename will be used.", "-o", false);
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
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    build_configuration build_config;

    auto m = parser.get<uint64_t>("m");
    build_config.m = m;
    if (parser.get<uint64_t>("RAM")) {
        build_config.ram_limit_in_GiB = parser.get<uint64_t>("RAM");
    }
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    if (parser.parsed("num_threads")) {
        build_config.num_threads = parser.get<uint64_t>("num_threads");
    }
    build_config.verbose = parser.get<bool>("verbose");

    build_config.base_filename = parser.get<std::string>("input_basename");
    assert(build_config.base_filename.string().length() != 0);
    build_config.output_filename = build_config.base_filename;
    if (parser.parsed("output_basename")) {
        build_config.output_filename = parser.get<std::string>("output_basename");
    }
    build_config.output_filename.replace_extension(constants::hfur_filename_extension);

    essentials::logger("loading the index...");
    hfur_index_t::loader loader(build_config);
    loader.load();

    essentials::logger("DONE");

    if (build_config.verbose) {
        util::timed_phase::print_breakdown();
    }

    return 0;
}

template <typename Index>
uint64_t probabilistic_check(Index& index, double const file_prob, double const kmer_prob,
                             uint32_t const num_threads, const bool verbose) {
    const uint64_t k = index.k();
    const uint32_t num_colors = index.num_colors();

    std::atomic<uint32_t> next_color_id(0);
    std::atomic<uint64_t> total_files_sampled(0);
    std::atomic<uint64_t> total_kmers_checked(0), total_raw_kmers(0);
    std::atomic<uint64_t> num_errors(0);
    std::mutex out_mtx;

    std::cout << "\rChecked 0 files, 0 k-mers" << std::flush;

    auto worker = [&](const int thread_id) {
        std::random_device rd;
        std::mt19937 gen(rd() ^ (static_cast<uint32_t>(thread_id) << 16));

        std::bernoulli_distribution file_dis(file_prob);
        std::geometric_distribution<size_t> kmer_dis(kmer_prob);

        while (true) {
            uint32_t color = next_color_id.fetch_add(1, std::memory_order_relaxed);
            if (color >= num_colors) break;

            if (!file_dis(gen)) continue;

            total_files_sampled.fetch_add(1, std::memory_order_relaxed);
            std::string_view filename = index.filename(color);

            try {
                const std::vector file_vec = {std::string(filename)};

                fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(file_vec, 1, 1);
                parser.start();

                auto rg = parser.getReadGroup();
                while (parser.refill(rg)) {
                    for (const auto& record : rg) {
                        const std::string& seq = record.seq;
                        if (seq.length() < k) continue;

                        const uint64_t num_kmers = seq.length() - k + 1;
                        total_raw_kmers += num_kmers;

                        for (uint64_t i = kmer_dis(gen); i < num_kmers; i += 1 + kmer_dis(gen)) {
                            std::string_view kmer_view(seq.data() + i, k);
                            if (kmer_view.find('N') != std::string_view::npos) continue;
                            ++total_kmers_checked;

                            std::vector<uint32_t> color_set_ids;
                            index.fetch_color_set_ids(std::string(kmer_view), color_set_ids);
                            assert(color_set_ids.size() <= 1);
                            if (color_set_ids.empty()) {
                                std::lock_guard lock(out_mtx);
                                const auto msg = std::format("[File {}] K-mer {} not found",
                                                             filename, kmer_view);
                                std::cout << msg << std::endl;
                                ++num_errors;
                                continue;
                            }

                            auto color_set = index.color_set(color_set_ids.front());
                            while (*color_set < color) {
                                ++color_set;
                            }

                            if (*color_set != color) {
                                std::lock_guard lock(out_mtx);
                                std::cout << std::format(
                                    "[File {}] k-mer {} not found (exp: {}, got:{})\n", filename,
                                    kmer_view, color, *color_set);
                                ++num_errors;
                            }
                        }
                    }
                }
                parser.stop();
                const auto progress =
                    std::format("\rChecked {} files, {} k-mers", total_files_sampled.load(),
                                total_kmers_checked.load());
                std::cout << progress << std::flush;
            } catch (const std::exception& e) {
                if (verbose) {
                    const auto msg = std::format("[Thread {}] Error reading file {}: {}", thread_id,
                                                 filename, e.what());
                    std::cerr << msg << std::endl;
                }
            }
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (uint64_t t = 0; t < num_threads; ++t) {
        threads.emplace_back(worker, t);
    }

    for (auto& t : threads) {
        t.join();
    }

    if (verbose) {
        const auto rpt_head = "\n--- Probabilistic Check Report ---\n";
        const auto rpt_ln0 =
            std::format("Files sampled: {} / {} ({}%)\n", total_files_sampled.load(), num_colors,
                        100. * total_files_sampled / num_colors);
        const auto rpt_ln1 =
            std::format("Total kmers checked: {} / {} ({}%)\n", total_kmers_checked.load(),
                        total_raw_kmers.load(), 100. * total_kmers_checked / total_raw_kmers);
        const auto rpt_foot = "----------------------------------\n";

        std::cout << rpt_head + rpt_ln0 + rpt_ln1 + rpt_foot << std::endl;
    }

    return num_errors;
}

int probabilistic_check(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index to be checked for correctness.", "-i", true);
    parser.add(
        "file_prob",
        "Probability of each file to be checked. Value must be in (0, 1]. (1 means all files)",
        "-q", true);
    parser.add(
        "kmer_prob",
        "Probability of each kmer to be checked. Value must be in (0, 1]. (1 means all kmers)",
        "-p", true);
    parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    parser.add("verbose", "Verbose output during processing (default is false).", "--verbose",
               false, true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    std::filesystem::path index_filename = parser.get<std::string>("index_filename");
    auto file_prob = parser.get<double>("file_prob");
    auto kmer_prob = parser.get<double>("kmer_prob");
    bool verbose = parser.get<bool>("verbose");
    uint64_t num_threads = parser.parsed("num_threads") ? parser.get<uint64_t>("num_threads") : 1;

    if (file_prob < 0.0 || file_prob > 1.0 || kmer_prob <= 0.0 || kmer_prob > 1.0) {
        throw std::invalid_argument(
            "Probabilities must be within valid bounds (0 < kmer_prob <= 1).");
    }

    std::variant<hfur_index_t, mdfur_index_t, mfur_index_t, dfur_index_t> index;
    if (is_meta_diff(index_filename)) {
        index = mdfur_index_t();
    } else if (is_meta(index_filename)) {
        index = mfur_index_t();
    } else if (is_diff(index_filename)) {
        index = dfur_index_t();
    } else if (is_hybrid(index_filename)) {
        index = hfur_index_t();
    } else {
        std::cerr << "Wrong index filename supplied." << std::endl;
        return 1;
    }

    uint64_t with_errors = 0;
    std::visit(
        [&index_filename, &with_errors, file_prob, kmer_prob, num_threads, verbose](auto&& index) {
            if (verbose) essentials::logger("*** START: loading the base index");
            essentials::mmap(index, index_filename.c_str());
            if (verbose) essentials::logger("*** DONE: loading the base index");

            with_errors = probabilistic_check(index, file_prob, kmer_prob, num_threads, verbose);
        },
        index);

    if (with_errors) {
        essentials::logger("*** Completed with errors, try to rebuild the index");
    } else {
        essentials::logger("*** Completed successfully!");
    }

    return 0;
}
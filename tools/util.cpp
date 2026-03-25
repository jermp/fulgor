using namespace fulgor;

bool is_meta(std::string const& index_filename) {
    return sshash::util::ends_with(index_filename,
                                   constants::mfur_filename_extension);
}

bool is_meta_diff(std::string const& index_filename) {
    return sshash::util::ends_with(index_filename,
                                   constants::mdfur_filename_extension);
}

bool is_diff(std::string const& index_filename) {
    return sshash::util::ends_with(index_filename,
                                   constants::dfur_filename_extension);
}

bool is_hybrid(std::string const& index_filename) {
    return sshash::util::ends_with(index_filename, constants::hfur_filename_extension);
}

template <typename FulgorIndex>
void verify(std::string const& index_filename) {
    FulgorIndex index;
    essentials::version_number vnum(constants::current_version_number::x,   //
                                    constants::current_version_number::y,   //
                                    constants::current_version_number::z);  //
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
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    index.print_stats();
}

template <typename FulgorIndex>
void print_filenames(std::string const& index_filename) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    for (uint64_t i = 0; i != index.num_colors(); ++i) {
        std::cout << i << '\t' << index.filename(i) << '\n';
    }
}

template <typename FulgorIndex>
void dump(std::string const& index_filename, std::string const& basename) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    index.dump(basename);
}

template<typename BaseIndex, typename TargetIndex>
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

    if (target.get_k2u().num_contigs() != base.get_k2u().num_contigs()) {
        std::cout << "Number of contigs mismatch" << std::endl;
        return 1;
    }
    const uint64_t num_unitigs = base.get_k2u().num_contigs();

    if (target.get_k2u().size() != base.get_k2u().size()) {
        std::cout << "Number of kmers mismatch" << std::endl;
        return 1;
    }
    const uint64_t num_kmers = base.get_k2u().size();

    std::vector<uint32_t> base_to_target(num_colors);
    {
        std::vector<uint32_t> base_permutation(num_colors), target_permutation(num_colors);
        std::iota(base_permutation.begin(), base_permutation.end(), 0);
        std::iota(target_permutation.begin(), target_permutation.end(), 0);
        std::sort(base_permutation.begin(), base_permutation.end(), [&](uint32_t a, uint32_t b) {
            return base.get_filenames()[a] < base.get_filenames()[b];
        });
        std::sort(target_permutation.begin(), target_permutation.end(), [&](uint32_t a, uint32_t b) {
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

    auto exe = [&](uint64_t unitig_begin, uint64_t unitig_end) {
        for (uint64_t unitig_id = unitig_begin; unitig_id != unitig_end; ++unitig_id) {
            if (verbose && ++num_checked_unitigs % 1000 == 0) {
                std::cout << "\033[3A"; // Move up 3 lines

                std::cout << "\033[2K"; // Clear line
                std::cout << "Unitigs checked:     " << num_checked_unitigs << "/" << num_unitigs << std::endl;

                std::cout << "\033[2K"; // Clear line
                std::cout << "Color sets checked:  " << num_checked_color_sets << "/" << num_color_sets << std::endl;

                std::cout << "\033[2K"; // Clear line
                std::cout << "Kmers checked:       " << num_checked_kmers << "/" << num_kmers << std::endl;
            }

            auto it = target.get_k2u().at_contig_id(unitig_id);
            auto [_, kmer] = it.next();
            const uint64_t base_contig_id = base.get_k2u().lookup_advanced(kmer.c_str()).contig_id;
            const uint64_t target_contig_id = target.get_k2u().lookup_advanced(kmer.c_str()).contig_id;
            ++num_checked_kmers;

            while (it.has_next()) {
                ++num_checked_kmers;
                auto [_, kmer] = it.next();
                const uint64_t curr_target_contig_id = target.get_k2u().lookup_advanced(kmer.c_str()).contig_id;
                const uint64_t curr_base_contig_id = base.get_k2u().lookup_advanced(kmer.c_str()).contig_id;
                if (target_contig_id != curr_target_contig_id) { // should never happen
                    std::cerr << "\033[1;31m" << "expected unitig " << target_contig_id << " but found " << curr_target_contig_id
                              << "\033[0m" << std::endl;
                }
                if (base_contig_id != curr_base_contig_id) {
                    std::cerr << "\033[1;31m" << "expected unitig " << base_contig_id << " but found " << curr_base_contig_id
                              << "\033[0m" << std::endl;
                }
            }

            uint64_t base_color_set_id = base.u2c(base_contig_id);
            uint64_t target_color_set_id = target.u2c(target_contig_id);

            if (checked_color_sets[target_color_set_id]) continue;
            checked_color_sets[target_color_set_id] = true;
            ++num_checked_color_sets;

            auto base_it = base.color_set(base_color_set_id);
            auto target_it = target.color_set(target_color_set_id);

            if (target_it.size() != base_it.size()) {
                std::cerr << "\033[1;31m" << "Error while checking color set " << target_color_set_id
                          << ", different sizes: expected " << base_it.size()
                          << " but got " << target_it.size() << "\033[0m" << std::endl;
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
                    std::cerr << "\033[1;31m" << "Error while checking color set " << target_color_set_id
                              << ", mismatch at position " << j << ": expected " << base_to_target[base_val]
                              << " but got " << target_val << "\033[0m" << std::endl;
                }
            }
        }
    };

    kmeans::thread_pool threads(num_threads);
    const uint64_t load_per_thread = num_unitigs / (num_threads << 10);
    uint64_t start = 0, end = load_per_thread;
    while (end < num_unitigs) {
        threads.enqueue([&, start, end]{ exe(start, std::min(end, num_unitigs)); });
        start = end;
        end = std::min(end + load_per_thread, num_unitigs);
    }
    threads.enqueue([&, start, num_unitigs]{ exe(start, num_unitigs); }); // last one
    threads.wait();

    std::cout << "\033[3A"; // Move up 3 lines

    std::cout << "\033[2K"; // Clear line
    std::cout << "Unitigs checked:     " << num_checked_unitigs << "/" << num_unitigs << std::endl;

    std::cout << "\033[2K"; // Clear line
    std::cout << "Color sets checked:  " << num_checked_color_sets << "/" << num_color_sets << std::endl;

    std::cout << "\033[2K"; // Clear line
    std::cout << "Kmers checked:       " << num_checked_kmers << "/" << num_kmers << std::endl;

    return 0;
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
               "Output basename. If the supplied basename is F, the output will consist in three "
               "files: F.unitigs.fa, F.color_sets.txt, and F.metadata.txt. (If this is not "
               "supplied, the basename of the index is used instead.)",
               "-o", false);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    auto index_filename = parser.get<std::string>("index_filename");
    std::string output_basename("");
    if (parser.parsed("output_basename")) {
        output_basename = parser.get<std::string>("output_basename");
        assert(output_basename.length() != 0);
    }

    if (is_meta_diff(index_filename)) {
        std::string basename{index_filename.data(),
                             index_filename.length() -
                                 constants::mdfur_filename_extension.length() -
                                 1};
        dump<mdfur_index_t>(
            index_filename, output_basename.length() == 0 ? basename : output_basename);
    } else if (is_meta(index_filename)) {
        std::string basename{index_filename.data(),
                             index_filename.length() -
                                 constants::mfur_filename_extension.length() - 1};
        dump<mfur_index_t>(index_filename,
                              output_basename.length() == 0 ? basename : output_basename);
    } else if (is_diff(index_filename)) {
        std::string basename{index_filename.data(),
                             index_filename.length() -
                                 constants::dfur_filename_extension.length() - 1};
        dump<dfur_index_t>(index_filename,
                                      output_basename.length() == 0 ? basename : output_basename);
    } else if (is_hybrid(index_filename)) {
        std::string basename{
            index_filename.data(),
            index_filename.length() - constants::hfur_filename_extension.length() - 1};
        dump<hfur_index_t>(index_filename,
                         output_basename.length() == 0 ? basename : output_basename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int check(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("base_filename", "The *correct* Fulgor index to be checked against.", "--base", true);
    parser.add("target_filename","The Fulgor index to be checked for correctness. Cannot be a .fur index.", "--target", true);
    parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    parser.add("verbose", "Verbose output during query (default is false).", "--verbose", false,
           true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    auto base_filename = parser.get<std::string>("base_filename");
    auto target_filename = parser.get<std::string>("target_filename");
    bool verbose = parser.get<bool>("verbose");
    uint64_t num_threads = parser.parsed("num_threads") ? parser.get<uint64_t>("num_threads") : 1;

    std::variant<hfur_index_t, mdfur_index_t, mfur_index_t, dfur_index_t> base_index;
    if (is_meta_diff(base_filename)) {
        base_index = mdfur_index_t();
    } else if (is_meta(base_filename)) {
        base_index = mfur_index_t();
    } else if (is_diff(base_filename)) {
        base_index = dfur_index_t();
    } else if (is_hybrid(base_filename)) {
        base_index = hfur_index_t();
    } else {
        std::cerr << "Wrong base index filename supplied." << std::endl;
        return 1;
    }

    std::variant<mdfur_index_t, mfur_index_t, dfur_index_t> target_index;
    if (is_meta_diff(target_filename)) {
        target_index = mdfur_index_t();
    } else if (is_meta(target_filename)) {
        target_index = mfur_index_t();
    } else if (is_diff(target_filename)) {
        target_index = dfur_index_t();
    } else {
        std::cerr << "Wrong target index filename supplied." << std::endl;
        return 1;
    }

    std::visit([&base_filename, &target_filename, num_threads, verbose]
                  (auto&& base, auto&& target) {
        if (verbose) essentials::logger("*** START: loading the base index");
        essentials::load(base, base_filename.c_str());
        if (verbose) essentials::logger("*** DONE: loading the base index");

        if (verbose) essentials::logger("*** START: loading the target index");
        essentials::load(target, target_filename.c_str());
        if (verbose) essentials::logger("*** DONE: loading the target index");

        check(base, target, num_threads, verbose);
    }, base_index, target_index);

    return 0;
}


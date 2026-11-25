#include "include/util.hpp"

#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
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

int verify(std::string const& index_filename) {
    const auto vnum = util::read_version_number(index_filename);
    std::cout << "read version number = " << vnum.to_string() << std::endl;
    util::check_version_number(vnum);
    essentials::logger("OK: Fulgor index is compatible with current library version.");
    return 0;
}

template <typename FulgorIndex>
void print_stats(std::string const& index_filename) {
    FulgorIndex index;
    util::load_index(index, index_filename, true, true);
    index.print_stats();
}

template <typename FulgorIndex>
void print_filenames(std::string const& index_filename) {
    FulgorIndex index;
    util::load_index(index, index_filename, true, true);
    for (uint64_t i = 0; i != index.num_colors(); ++i) {
        std::cout << i << '\t' << index.filename(i) << '\n';
    }
}

template <typename FulgorIndex>
void dump(std::string const& index_filename, build_configuration const& build_config) {
    FulgorIndex index;
    util::load_index(index, index_filename, true, true);
    index.dump(build_config);
}

int verify(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    const auto index_filename = parser.get<std::string>("index_filename");

    const bool known_extension = is_meta(index_filename) || is_meta_diff(index_filename) ||
                                 is_diff(index_filename) || is_hybrid(index_filename);
    if (!known_extension) {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }

    return verify(index_filename);
}

int stats(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    const auto index_filename = parser.get<std::string>("index_filename");
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

    const auto index_filename = parser.get<std::string>("index_filename");
    std::filesystem::path output_basename;
    if (parser.parsed("output_basename")) {
        output_basename = parser.get<std::string>("output_basename");
    }
    if (output_basename.empty()) {
        output_basename = std::filesystem::path(index_filename).stem();
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

    const auto m = parser.get<uint64_t>("m");
    build_config.m = m;
    if (parser.parsed("RAM")) {
        const uint64_t ram_limit = parser.get<uint64_t>("RAM");
        if (ram_limit == 0) {
            std::cerr << "-g must be > 0" << std::endl;
            return 1;
        }
        build_config.ram_limit_in_GiB = ram_limit;
    }
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    if (parser.parsed("num_threads")) {
        build_config.num_threads = parser.get<uint64_t>("num_threads");
        if (build_config.num_threads == 0) {
            std::cerr << "Number of threads must be >= 1." << std::endl;
            return 1;
        }
    }
    build_config.verbose = parser.get<bool>("verbose");

    build_config.base_filename = parser.get<std::string>("input_basename");
    assert(!build_config.base_filename.string().empty());
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
uint64_t check(Index& index, double const file_prob, double const kmer_prob,
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
                            ++total_kmers_checked;
                            if (!util::is_dna(kmer_view)) continue;

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
                            color_set.next_geq(color);

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

int check(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index to be checked for correctness.", "-i", true);
    parser.add(
        "file_prob",
        "Probability of each file to be checked. Value must be in (0, 1]. (1 means all files)",
        "-fp", true);
    parser.add(
        "kmer_prob",
        "Probability of each kmer to be checked. Value must be in (0, 1]. (1 means all kmers)",
        "-kp", true);
    parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    parser.add("mmap",
               "Use memory mapping instead of loading the whole index in RAM. Use this option if "
               "the index does not fit in memory. WARNING: significantly slows down query speed. "
               "Not supported for meta-differential indexes.",
               "--mmap", false, true);
    parser.add("verbose", "Verbose output during processing (default is false).", "--verbose",
               false, true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    std::filesystem::path index_filename = parser.get<std::string>("index_filename");
    auto file_prob = parser.get<double>("file_prob");
    auto kmer_prob = parser.get<double>("kmer_prob");
    bool verbose = parser.get<bool>("verbose");
    bool mmap = parser.get<bool>("mmap");

    uint64_t num_threads = 1;
    if (parser.parsed("num_threads")) {
        num_threads = parser.get<uint64_t>("num_threads");
        if (num_threads == 0) {
            std::cerr << "Number of threads must be >= 1." << std::endl;
            return 1;
        }
    }

    if (file_prob <= 0.0 || file_prob > 1.0 || kmer_prob <= 0.0 || kmer_prob > 1.0) {
        throw std::invalid_argument(
            "Probabilities must be within valid bounds (0 < file_prob, kmer_prob <= 1).");
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
        [&index_filename, &with_errors, file_prob, kmer_prob, num_threads, verbose,
         mmap](auto&& index) {
            util::load_index(index, index_filename, mmap, verbose);

            with_errors = check(index, file_prob, kmer_prob, num_threads, verbose);
        },
        index);

    if (with_errors) {
        essentials::logger("*** Completed with errors, try to rebuild the index");
    } else {
        essentials::logger("*** Completed successfully!");
    }

    return 0;
}

#include "src/ps_utils.cpp"

int uncompress_psa_output(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("input_filename", "The pseudoalignment file to uncompress.", "-i", true);
    // parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    // parser.add("verbose", "Verbose output during processing (default is false).", "--verbose",
    //            false, true);
    parser.add("output_filename",
               "File where output will be written. You can specify \"/dev/stdout\" to write "
               "output to stdout. In this case, it is also recommended to use the --verbose flag "
               "to avoid printing status messages to stdout.",
               "-o", true);
    parser.add("format",
               "Format of the output file: ascii or binary"
               " (default is ascii).",
               "--format", false);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto output_filename = parser.get<std::string>("output_filename");
    auto output_format = parser.parsed("format") ? parser.get<std::string>("format") : "ascii";

    try {
        auto filename = parser.get<std::string>("input_filename");

        std::ifstream file(filename, std::ios::binary | std::ios::ate);
        if (!file.is_open()) throw std::runtime_error("cannot open input file: " + filename);
        const std::streamoff end = file.tellg();
        constexpr uint64_t header_bytes = 3 * sizeof(uint64_t);
        if (end < 0) throw std::runtime_error("cannot determine the size of: " + filename);
        if (static_cast<uint64_t>(end) < header_bytes)
            throw std::runtime_error("input file is too small: " + filename);

        file.seekg(0, std::ios::beg);
        uint64_t num_colors = 0, num_bits = 0, num_words = 0;
        file.read(reinterpret_cast<char*>(&num_colors), sizeof(num_colors));
        file.read(reinterpret_cast<char*>(&num_bits), sizeof(num_bits));
        file.read(reinterpret_cast<char*>(&num_words), sizeof(num_words));
        if (!file) throw std::runtime_error("truncated header in: " + filename);
        if (num_colors == 0 || num_colors > std::numeric_limits<uint32_t>::max())
            throw std::runtime_error("invalid num_colors in header: " + std::to_string(num_colors));

        const uint64_t expected_words = num_bits / 64 + (num_bits % 64 != 0);
        if (num_words != expected_words)
            throw std::runtime_error("inconsistent num_bits/num_words in header");
        const uint64_t file_bytes = static_cast<uint64_t>(end);
        if (num_words > (file_bytes - header_bytes) / sizeof(uint64_t))
            throw std::runtime_error("truncated compressed payload in: " + filename);

        cdbg::streaming_bit_vector sbv(std::move(file), 8);

        auto formatter_var = [&]() -> std::variant<psa_ascii_formatter, psa_binary_formatter> {
            if (output_format == "ascii") {
                return std::variant<psa_ascii_formatter, psa_binary_formatter>{
                    std::in_place_type<psa_ascii_formatter>, output_filename};
            }
            if (output_format == "binary") {
                return std::variant<psa_ascii_formatter, psa_binary_formatter>{
                    std::in_place_type<psa_binary_formatter>, output_filename};
            }
            throw std::invalid_argument("Unknown output format. Supported formats: ascii, binary.");
        }();

        std::visit(
            [&]<typename Formatter>(Formatter& formatter) {
                auto buffer = formatter.buffer();

                std::vector<uint32_t> result;
                while (!sbv.eof()) {
                    const uint64_t query_id = sbv.read_delta();
                    const uint64_t cs_size = sbv.read_delta();
                    if (query_id > std::numeric_limits<uint32_t>::max())
                        throw std::runtime_error("invalid query id in compressed input");
                    if (cs_size > num_colors)
                        throw std::runtime_error("invalid color-set size in compressed input");
                    result.clear();

                    if (cs_size == 0) {
                    } else if (cs_size < num_colors / 4) {
                        uint32_t prev_val = sbv.read_delta();
                        result.push_back(prev_val);
                        for (uint64_t i = 1; i < cs_size; ++i) {
                            auto val = static_cast<uint32_t>(sbv.read_delta() + prev_val + 1);
                            result.push_back(val);
                            prev_val = val;
                        }
                    } else if (cs_size < (3 * num_colors) / 4) {
                        const uint32_t init_pos = sbv.position();
                        for (uint32_t i = 0; i < cs_size; ++i) {
                            result.push_back(sbv.next() - init_pos);
                        }
                        sbv.skip_to(init_pos + num_colors);
                    } else {
                        const uint64_t size = num_colors - cs_size;
                        uint32_t curr = 0;
                        uint32_t next = 0;
                        for (uint64_t i = 0; i < size; ++i) {
                            next += sbv.read_delta();
                            while (curr < next) {
                                result.push_back(curr);
                                curr++;
                            }
                            curr = next + 1;
                            ++next;
                        }
                        while (curr < num_colors) {
                            result.push_back(curr);
                            curr++;
                        }
                    }

                    buffer.write(static_cast<uint32_t>(query_id), result);
                }
            },
            formatter_var);
    } catch (std::exception const& e) {
        std::cerr << "extract-output: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
#include <iostream>
#include <fstream>
#include <sstream>
#include <variant>

#include "src/ps_full_intersection.cpp"
#include "src/ps_threshold_union.cpp"

using namespace fulgor;

enum class pseudoalignment_algorithm : uint8_t { FULL_INTERSECTION, THRESHOLD_UNION };

std::string to_string(pseudoalignment_algorithm algo, double threshold) {
    std::string o;
    switch (algo) {
        case pseudoalignment_algorithm::FULL_INTERSECTION:
            o = "full-intersection";
            break;
        case pseudoalignment_algorithm::THRESHOLD_UNION:
            o = "threshold-union (threshold = " + std::to_string(threshold) + ")";
            break;
    }
    return o;
}

template <typename FulgorIndex, typename Formatter, typename QueryReader>
int pseudoalign(FulgorIndex const& index, QueryReader& query_reader,
                std::atomic<uint64_t>& num_reads, std::atomic<uint64_t>& num_mapped_reads,
                pseudoalignment_algorithm algo, const double threshold, Formatter& formatter,
                std::mutex& iomut, const bool verbose)  //
{
    auto output_buffer = formatter.buffer();
    std::vector<uint32_t> tmp, colors;  // result of pseudoalignment
    std::vector<uint32_t> color_set_ids;
    std::stringstream ss;

    auto qg = query_reader.get_query_group();

    while (qg.refill()) {
        while (qg.has_next()){
            typename QueryReader::query_t query;
            qg.value(query);

            switch (algo) {
                case pseudoalignment_algorithm::FULL_INTERSECTION:
                    // index.fetch_color_set_ids(record.seq, color_set_ids);
                    index.pseudoalign_full_intersection(query.cids, colors, tmp);
                    break;
                case pseudoalignment_algorithm::THRESHOLD_UNION:
                    index.pseudoalign_threshold_union(query.seq, colors, threshold);
                    cerr << "THRESHOLD_UNION not supported!! Needs Fixing!!" << endl;
                    exit(1);
                    break;
                default:
                    break;
            }

            if constexpr (std::is_same_v<util::preprocessed_query_reader, QueryReader>) {
                for (auto qid: query.ids) {
                    num_reads += 1;
                    output_buffer.write(qid, colors);
                    if (!colors.empty()) {
                        num_mapped_reads += 1;
                    }
                }
            } else {
                num_reads += 1;
                output_buffer.write(query.id, colors);
                if (!colors.empty()) {
                    num_mapped_reads += 1;
                }
            }

            colors.clear();

            if (verbose and num_reads > 0 and num_reads % 1000000 == 0) {
                iomut.lock();
                std::cout << "processed " << num_reads << " reads" << std::endl;
                iomut.unlock();
            }
            qg.next();
        }
    }
    return 0;
}

template <typename FulgorIndex, typename Formatter,  typename QueryReader>
int pseudoalign(FulgorIndex& index, QueryReader& query_reader,
                Formatter& formatter, uint64_t num_threads, double threshold,
                pseudoalignment_algorithm ps_alg, const bool verbose) {

    std::cerr << "query mode : " << to_string(ps_alg, threshold) << "\n";

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> t;
    t.start();

    std::atomic<uint64_t> num_mapped_reads{0};
    std::atomic<uint64_t> num_reads{0};

    assert(num_threads >= 2);

    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    std::mutex iomut;

    for (uint64_t i = 1; i != num_threads; ++i) {
        workers.push_back(std::thread([&index, &query_reader, &num_reads, &num_mapped_reads, ps_alg,
                                       threshold, &formatter, &iomut, verbose]() {
            pseudoalign(index, query_reader, num_reads, num_mapped_reads, ps_alg, threshold, formatter,
                        iomut, verbose);
        }));
    }

    for (auto& w : workers) w.join();

    t.stop();
    if (verbose) essentials::logger("DONE");

    if (verbose) {
        std::cout << "mapped " << num_reads << " reads" << std::endl;
        std::cout << "elapsed = " << t.elapsed() << " millisec / ";
        std::cout << t.elapsed() / 1000 << " sec / ";
        std::cout << t.elapsed() / 1000 / 60 << " min / ";
        std::cout << (t.elapsed() * 1000) / num_reads << " musec/read" << std::endl;
        std::cout << "num_mapped_reads " << num_mapped_reads << "/" << num_reads << " ("
                  << (num_mapped_reads * 100.0) / num_reads << "%)" << std::endl;
    }

    return 0;
}

template <typename FulgorIndex, typename Formatter>
void fetch_and_deduplicate_sets(const std::string& query_filename,
                                Formatter& output_formatter,
                                std::string& tmp_filename,
                                FulgorIndex& index,
                                uint64_t num_threads,
                                bool verbose) {
    auto output_buffer = output_formatter.buffer();
    if (verbose) essentials::logger("*** START: fetching color set ids");

    std::ofstream tmp_file(tmp_filename, std::ios::binary);
    auto query_filenames = std::vector({query_filename});
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(query_filenames, num_threads,
                                                             num_threads - 1);
    rparser.start();
    std::vector<std::thread> workers;
    std::mutex outfile_mut, iomut;

    constexpr int32_t buff_thresh = 50;
    std::atomic<uint64_t> num_reads = 0;
    auto fetch = [&rparser, &index, &tmp_file, &outfile_mut, &iomut, &num_reads] () {
        uint32_t buff_size = 0;
        std::vector<uint32_t> color_set_ids;
        std::stringstream ss;

        auto rg = rparser.getReadGroup();
        while (rparser.refill(rg)) {
            uint32_t read_id = rg.chunk_frag_offset().frag_idx;

            for (auto const& record: rg) {
                index.fetch_color_set_ids(record.seq, color_set_ids);

                buff_size += 1;
                num_reads += 1;

                ss.write(reinterpret_cast<char*>(&read_id), sizeof(read_id));
                uint32_t num_color_sets = static_cast<uint32_t>(color_set_ids.size());
                ss.write(reinterpret_cast<char*>(&num_color_sets), sizeof(num_color_sets));
                if (num_color_sets > 0) {
                    ss.write(reinterpret_cast<char*>(color_set_ids.data()), num_color_sets * sizeof(color_set_ids[0]));
                }

                color_set_ids.clear();
                if (num_reads > 0 and num_reads % 1000000 == 0) {
                    iomut.lock();
                    std::cout << "fetched " << num_reads << " reads" << std::endl;
                    iomut.unlock();
                }
                if (buff_size > buff_thresh) {
                    std::string outs = ss.str();
                    ss.str("");
                    outfile_mut.lock();
                    tmp_file.write(outs.data(), outs.size());
                    outfile_mut.unlock();
                    buff_size = 0;
                }
                ++read_id;
            }
        }
        if (buff_size > 0) {
            std::string outs = ss.str();
            ss.str("");
            outfile_mut.lock();
            tmp_file.write(outs.data(), outs.size());
            outfile_mut.unlock();
            buff_size = 0;
        }
    };

    for (uint64_t i = 0; i < num_threads; ++i) {
        workers.push_back(std::thread(fetch));
    }
    for (auto& w : workers) w.join();
    rparser.stop();
    tmp_file.close();

    if (verbose) essentials::logger("*** DONE: fetching color set ids");
    if (verbose) essentials::logger("*** START: deduplicating queries");

    std::ifstream ifile(tmp_filename, std::ios::binary);
    std::vector<std::vector<uint32_t>> queries;
    std::vector<uint32_t> tmp;
    uint32_t read_num = 0;
    while (ifile.read(reinterpret_cast<char*>(&read_num), sizeof(read_num))) {
        uint32_t num_colors = 0;
        ifile.read(reinterpret_cast<char*>(&num_colors), sizeof(num_colors));

        tmp.resize(num_colors + 1);
        tmp[0] = read_num;
        ifile.read(reinterpret_cast<char*>(&tmp[1]), num_colors * sizeof(num_colors));
        if (tmp.size() > 1) {
            queries.push_back(tmp);
            tmp.clear();
        } else {
            // just write out the unmapped reads here
            output_buffer.write(read_num, {});
        }
    }

    if (queries.empty()) { return; }

    std::sort(queries.begin(), queries.end(),
    [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) -> bool {
      return std::lexicographical_compare(a.begin() + 1, a.end(), b.begin() + 1,
                                          b.end());
    });

    auto curr = queries.begin();
    auto next = curr++;
    size_t identical_lists = 0;
    uint64_t identical_sizes = 0;
    while (next < queries.end()) {
        if ((curr->size() == next->size()) and
            std::equal(curr->begin() + 1, curr->end(), next->begin() + 1)) {
            identical_sizes += next->size() - 1;
            next->resize(1);  // retain only the id.
            ++identical_lists;
            ++next;
        } else {
            curr = next;
            ++next;
        }
    }

    std::cerr << "number of identical lists = " << identical_lists << " (ignoring "
              << identical_sizes << " set ids)\n";
    std::cerr << "sorted!\n";
    ifile.close();
    std::ofstream ofile(tmp_filename, std::ios::trunc | std::ios::binary);
    for (auto& query : queries) {
        uint32_t s = query.size();
        ofile.write(reinterpret_cast<char*>(&s), sizeof(s));
        if (s > 0) {
            ofile.write(reinterpret_cast<char*>(query.data()), sizeof(query[0]) * s);
        }
    }
    ofile.close();

    if (verbose) essentials::logger("*** DONE: deduplicating queries");
}

int pseudoalign(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("query_filename", "Query filename in FASTA/FASTQ format (optionally gzipped).", "-q",
               true);
    parser.add("output_filename",
               "File where output will be written. You can specify \"/dev/stdout\" to write "
               "output to stdout. In this case, it is also recommended to use the --verbose flag "
               "to avoid printing status messages to stdout.",
               "-o", true);
    parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    parser.add("verbose", "Verbose output during query (default is false).", "--verbose", false,
               true);
    parser.add("threshold",
               "Threshold for threshold_union algorithm. It must be a float in (0.0,1.0].", "-r",
               false);
    parser.add("deduplicate", "Removes duplicate queries before pseudoalignment (default is false)."
               " Only works on Full-Intersection. Creates a temporary file in the executable's directory.",
               "--deduplicate", false, true);
    parser.add("format", "Format of the output file. Must either ascii, binary, compressed"
               " (default is ascii).", "--format", false);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto query_filename = parser.get<std::string>("query_filename");
    auto output_filename = parser.get<std::string>("output_filename");

    bool deduplicate = parser.get<bool>("deduplicate");
    auto output_format = parser.parsed("format") ? parser.get<std::string>("format") : "ascii";

    uint64_t num_threads = 1;
    if (parser.parsed("num_threads")) num_threads = parser.get<uint64_t>("num_threads");
    if (num_threads == 1) {
        num_threads += 1;
        std::cerr
            << "1 thread was specified, but an additional thread will be allocated for parsing"
            << std::endl;
    }

    double threshold = constants::invalid_threshold;
    if (parser.parsed("threshold")) threshold = parser.get<double>("threshold");
    if (threshold == 0.0 or threshold > 1.0) {
        std::cerr << "threshold must be a float in (0.0,1.0]" << std::endl;
        return 1;
    }

    auto ps_alg = pseudoalignment_algorithm::FULL_INTERSECTION;
    if (threshold != constants::invalid_threshold) {
        if (deduplicate) {
            throw new std::runtime_error("Deduplication not available for threshold < 1.0. Remove --deduplicate flag.");
        }
        ps_alg = pseudoalignment_algorithm::THRESHOLD_UNION;
    }

    bool verbose = parser.get<bool>("verbose");
    if (verbose) util::print_cmd(argc, argv);

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
    }

    std::variant<std::monostate, util::psa_ascii_formatter, util::psa_binary_formatter, util::psa_compressed_formatter> formatter;
    if (output_format == "ascii") {
        formatter.emplace<util::psa_ascii_formatter>(output_filename);
    } else if (output_format == "binary") {
        formatter.emplace<util::psa_binary_formatter>(output_filename);
    } else if (output_format == "compressed") {
        formatter.emplace<util::psa_compressed_formatter>(output_filename);
    } else {
        throw std::runtime_error("Unknown output format. Supported formats: ascii, binary, compressed.");
    }

    std::string tmp_filename = "queries.tmp";

    std::visit([&index_filename, &query_filename, &output_filename, &tmp_filename,
                deduplicate, num_threads, threshold, ps_alg, verbose]
                      (auto&& index, auto&& formatter) {
        if (verbose) essentials::logger("*** START: loading the index");
        essentials::load(index, index_filename.c_str());
        if (verbose) essentials::logger("*** DONE: loading the index");

        if (verbose) essentials::logger("performing queries from file '" + query_filename + "'...");

        if constexpr (std::is_same_v<std::decay_t<decltype(formatter)>, util::psa_compressed_formatter>) {
            formatter.set_num_colors(index.num_colors());
        }
        if constexpr (!std::is_same_v<std::decay_t<decltype(formatter)>, std::monostate>) {
            std::ofstream out(output_filename);

            if (deduplicate) {
                fetch_and_deduplicate_sets(query_filename, formatter, tmp_filename, index, num_threads, verbose);
                util::preprocessed_query_reader query_reader(tmp_filename, num_threads);
                pseudoalign(index, query_reader, formatter, num_threads, threshold, ps_alg, verbose);
            } else {
                util::fastq_query_reader query_reader(query_filename, num_threads, index);
                pseudoalign(index, query_reader, formatter, num_threads, threshold, ps_alg, verbose);
            }
        }

    }, index, formatter);

    std::remove(tmp_filename.c_str());

    return 0;
}

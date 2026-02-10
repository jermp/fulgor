#include <iostream>
#include <fstream>
#include <sstream>

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

template <typename FulgorIndex>
int pseudoalign(FulgorIndex const& index, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser,
                std::atomic<uint64_t>& num_reads, std::atomic<uint64_t>& num_mapped_reads,
                pseudoalignment_algorithm algo, const double threshold, std::ofstream& out_file,
                std::mutex& iomut, std::mutex& ofile_mut, const bool verbose)  //
{
    std::vector<uint32_t> tmp, colors;  // result of pseudoalignment
    std::vector<uint32_t> color_set_ids;
    std::stringstream ss;
    uint64_t buff_size = 0;
    constexpr uint64_t buff_thresh = 50;

    auto rg = rparser.getReadGroup();
    while (rparser.refill(rg)) {
        for (auto const& record : rg) {
            switch (algo) {
                case pseudoalignment_algorithm::FULL_INTERSECTION:
                    index.fetch_color_set_ids(record.seq, color_set_ids);
                    index.pseudoalign_full_intersection(color_set_ids, colors, tmp);
                    break;
                case pseudoalignment_algorithm::THRESHOLD_UNION:
                    index.pseudoalign_threshold_union(record.seq, colors, threshold);
                    break;
                default:
                    break;
            }
            buff_size += 1;
            if (!colors.empty()) {
                num_mapped_reads += 1;
                ss << record.name << '\t' << colors.size();
                for (auto c : colors) { ss << "\t" << c; }
                ss << '\n';
            } else {
                ss << record.name << "\t0\n";
            }
            num_reads += 1;
            colors.clear();
            if (verbose and num_reads > 0 and num_reads % 1000000 == 0) {
                iomut.lock();
                std::cout << "mapped " << num_reads << " reads" << std::endl;
                iomut.unlock();
            }
            if (buff_size > buff_thresh) {
                std::string outs = ss.str();
                ss.str("");
                ofile_mut.lock();
                out_file.write(outs.data(), outs.size());
                ofile_mut.unlock();
                buff_size = 0;
            }
        }
    }

    // dump anything left in the buffer
    if (buff_size > 0) {
        std::string outs = ss.str();
        ss.str("");
        ofile_mut.lock();
        out_file.write(outs.data(), outs.size());
        ofile_mut.unlock();
        buff_size = 0;
    }

    return 0;
}

template <typename FulgorIndex>
int pseudoalign(FulgorIndex& index, std::string const& query_filename,
                std::string const& output_filename, uint64_t num_threads, double threshold,
                pseudoalignment_algorithm ps_alg, const bool verbose) {

    std::cerr << "query mode : " << to_string(ps_alg, threshold) << "\n";

    std::ifstream is(query_filename.c_str());
    if (!is.good()) {
        std::cerr << "error in opening the file '" + query_filename + "'" << std::endl;
        return 1;
    }

    if (verbose) essentials::logger("performing queries from file '" + query_filename + "'...");
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> t;
    t.start();

    std::atomic<uint64_t> num_mapped_reads{0};
    std::atomic<uint64_t> num_reads{0};

    auto query_filenames = std::vector<std::string>({query_filename});
    assert(num_threads >= 2);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(query_filenames, num_threads,
                                                             num_threads - 1);

    rparser.start();
    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    std::mutex iomut;
    std::mutex ofile_mut;

    std::ofstream out_file;
    out_file.open(output_filename, std::ios::out | std::ios::trunc);
    if (!out_file) {
        std::cerr << "could not open output file " + output_filename << std::endl;
        return 1;
    }

    for (uint64_t i = 1; i != num_threads; ++i) {
        workers.push_back(std::thread([&index, &rparser, &num_reads, &num_mapped_reads, ps_alg,
                                       threshold, &out_file, &iomut, &ofile_mut, verbose]() {
            pseudoalign(index, rparser, num_reads, num_mapped_reads, ps_alg, threshold, out_file,
                        iomut, ofile_mut, verbose);
        }));
    }

    for (auto& w : workers) w.join();
    rparser.stop();

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

template <typename FulgorIndex>
void fetch_and_deduplicate_sets(const std::string& query_filename,
                                std::ofstream& out_file,
                                std::string& tmp_filename,
                                FulgorIndex& index,
                                uint64_t num_threads) {
    essentials::logger("*** START: fetching color set ids");

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
    };

    for (uint64_t i = 0; i < num_threads; ++i) {
        workers.push_back(std::thread(fetch));
    }
    for (auto& w : workers) w.join();
    rparser.stop();
    tmp_file.close();

    essentials::logger("*** DONE: fetching color set ids");
    essentials::logger("*** START: deduplicating queries");

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
            out_file << read_num << "\t0\n";
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

    essentials::logger("*** DONE: deduplicating queries");
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
        ps_alg = pseudoalignment_algorithm::THRESHOLD_UNION;
    }

    bool verbose = parser.get<bool>("verbose");
    if (verbose) util::print_cmd(argc, argv);

    std::variant<hfur_index_t, mdfur_index_t, mfur_index_t, dfur_index_t> index;
    if (sshash::util::ends_with(index_filename,
                                constants::mdfur_filename_extension)) {
        index = mdfur_index_t();
    } else if (sshash::util::ends_with(index_filename,
                                       constants::mfur_filename_extension)) {
        index = mfur_index_t();
    } else if (sshash::util::ends_with(index_filename,
                                      constants::dfur_filename_extension)) {
        index = dfur_index_t();
    } else if (sshash::util::ends_with(index_filename, constants::hfur_filename_extension)) {
        index = hfur_index_t();
    } else {
        std::cerr << "Wrong index filename supplied." << std::endl;
    }

    std::string tmp_filename = "queries.tmp";

    std::visit([&index_filename, &query_filename, &output_filename, &tmp_filename,
                deduplicate, num_threads, threshold, ps_alg, verbose]
                      (auto&& index) {
        essentials::logger("*** START: loading the index");
        essentials::load(index, index_filename.c_str());
        essentials::logger("*** DONE: loading the index");

        std::ofstream out(output_filename);
        if (deduplicate) {
            fetch_and_deduplicate_sets(query_filename, out, tmp_filename, index, num_threads);
        }
        pseudoalign(index, query_filename, output_filename, num_threads, threshold, ps_alg, verbose);
    }, index);

    std::remove(tmp_filename.c_str());

    return 1;
}

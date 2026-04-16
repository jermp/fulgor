#include <iostream>
#include <fstream>
#include <sstream>

#include "src/kmer_conservation.cpp"

using namespace fulgor;

template <typename FulgorIndex>
void kmer_conservation(FulgorIndex const& index,
                       fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser,
                       std::ofstream& out_file, std::mutex& ofile_mut,
                       query_options& options)  //
{
    std::vector<kmer_conservation_triple> kmer_conservation_info;
    std::stringstream ss;
    uint64_t buff_size = 0;
    constexpr uint64_t buff_thresh = 50;

    auto rg = rparser.getReadGroup();
    while (rparser.refill(rg)) {
        for (auto const& record : rg) {
            assert(record.seq.length() < (uint64_t(1) << 32));
            index.kmer_conservation(record.seq, kmer_conservation_info);
            buff_size += 1;
            if (!kmer_conservation_info.empty()) {
                ss << record.name << '\t' << kmer_conservation_info.size();
                for (auto kct : kmer_conservation_info) {
                    ss << "\t(" << kct.start_pos_in_query << ' ' << kct.num_kmers << ' '
                       << kct.color_set_id << ')';
                }
                ss << '\n';
            } else {
                ss << record.name << "\t0\n";
            }
            kmer_conservation_info.clear();
            options.increment_processed_reads();
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
}

template <typename FulgorIndex>
int kmer_conservation(std::string const& index_filename, std::string const& query_filename,
                      std::string const& output_filename, query_options& options) {
    FulgorIndex index;
    if (options.verbose) essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    if (options.verbose) essentials::logger("DONE");

    std::ifstream is(query_filename.c_str());
    if (!is.good()) {
        std::cerr << "error in opening the file '" + query_filename + "'" << std::endl;
        return 1;
    }

    if (options.verbose) {
        essentials::logger("performing queries from file '" + query_filename + "'...");
    }

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> t;
    t.start();

    const uint64_t num_threads = options.num_threads;
    auto query_filenames = std::vector<std::string>({query_filename});
    assert(num_threads >= 2);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(query_filenames, num_threads,
                                                             num_threads - 1);

    rparser.start();
    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    std::mutex ofile_mut;

    std::ofstream out_file;
    out_file.open(output_filename, std::ios::out | std::ios::trunc);
    if (!out_file) {
        std::cerr << "could not open output file " + output_filename << std::endl;
        return 1;
    }

    for (uint64_t i = 1; i != num_threads; ++i) {
        workers.push_back(std::thread([&index, &rparser, &out_file, &ofile_mut, &options]() {
            kmer_conservation(index, rparser, out_file, ofile_mut, options);
        }));
    }

    for (auto& w : workers) w.join();
    rparser.stop();

    t.stop();
    if (options.verbose) essentials::logger("DONE");

    if (options.verbose) {
        std::cout << "processed " << options.num_reads << " reads" << std::endl;
        std::cout << "elapsed = " << t.elapsed() << " millisec / ";
        std::cout << t.elapsed() / 1000 << " sec / ";
        std::cout << t.elapsed() / 1000 / 60 << " min / ";
        std::cout << (t.elapsed() * 1000) / options.num_reads << " musec/read" << std::endl;
    }

    return 0;
}

int kmer_conservation(int argc, char** argv) {
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
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto query_filename = parser.get<std::string>("query_filename");
    auto output_filename = parser.get<std::string>("output_filename");

    uint64_t num_threads = 1;
    if (parser.parsed("num_threads")) num_threads = parser.get<uint64_t>("num_threads");
    if (num_threads == 1) {
        num_threads += 1;
        std::cerr
            << "1 thread was specified, but an additional thread will be allocated for parsing"
            << std::endl;
    }

    bool verbose = parser.get<bool>("verbose");
    if (verbose) util::print_cmd(argc, argv);

    query_options options(verbose, num_threads);

    if (is_meta_diff(index_filename)) {
        return kmer_conservation<mdfur_index_t>(index_filename, query_filename, output_filename,
                                                options);
    } else if (is_meta(index_filename)) {
        return kmer_conservation<mfur_index_t>(index_filename, query_filename, output_filename,
                                               options);
    } else if (is_diff(index_filename)) {
        return kmer_conservation<dfur_index_t>(index_filename, query_filename, output_filename,
                                               options);
    } else if (is_hybrid(index_filename)) {
        return kmer_conservation<hfur_index_t>(index_filename, query_filename, output_filename,
                                               options);
    }

    std::cerr << "Wrong index filename supplied." << std::endl;

    return 1;
}

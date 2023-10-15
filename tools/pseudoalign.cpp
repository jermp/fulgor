#include <iostream>
#include <fstream>
#include <sstream>

#include "../external/CLI11.hpp"
#include "../external/sshash/include/gz/zip_stream.hpp"
#include "../external/FQFeeder/include/FastxParser.hpp"
#include "../external/FQFeeder/src/FastxParser.cpp"

#include "../piscem_psa/hit_searcher.cpp"
#include "../kallisto_psa/psa.cpp"
#include "../src/psa/full_intersection.cpp"
#include "../src/psa/threshold_union.cpp"

using namespace fulgor;

enum class pseudoalignment_algorithm : uint8_t {
    FULL_INTERSECTION,
    THRESHOLD_UNION,
    SKIPPING,
    SKIPPING_KALLISTO
};

std::string to_string(pseudoalignment_algorithm algo, double threshold) {
    std::string o;
    switch (algo) {
        case pseudoalignment_algorithm::FULL_INTERSECTION:
            o = "full-intersection";
            break;
        case pseudoalignment_algorithm::THRESHOLD_UNION:
            o = "threshold-intersection (t = " + std::to_string(threshold) + ")";
            break;
        case pseudoalignment_algorithm::SKIPPING:
            o = "skipping-pseudoalignment";
            break;
        case pseudoalignment_algorithm::SKIPPING_KALLISTO:
            o = "kallisto-pseudoalignment";
            break;
    }
    return o;
}

template <typename FulgorIndex>
int do_map(FulgorIndex const& index, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser,
           std::atomic<uint64_t>& num_reads, std::atomic<uint64_t>& num_mapped_reads,
           pseudoalignment_algorithm algo, const double threshold, std::ofstream& out_file,
           std::mutex& iomut, std::mutex& ofile_mut) {
    std::vector<uint32_t> colors;  // result of pseudo-alignment
    std::stringstream ss;
    int32_t buff_size{0};
    constexpr int32_t buff_thresh{100};

    if ((algo == pseudoalignment_algorithm::SKIPPING or
         algo == pseudoalignment_algorithm::SKIPPING_KALLISTO) and
        (index.get_dict().canonicalized())) {
        std::vector<uint32_t> unitig_ids;                           // for use with skipping
        std::vector<std::pair<projected_hits, int>> kallisto_hits;  // for use with kallisto psa

        piscem_psa::hit_searcher<FulgorIndex> hs(&index);
        sshash::streaming_query_canonical_parsing qc(&index.get_dict());

        auto get_hits_piscem_psa = [&qc, &hs](const std::string& seq,
                                              std::vector<uint32_t>& unitig_ids) -> void {
            hs.clear();
            auto had_hits = hs.get_raw_hits_sketch(seq, qc, true, false);
            if (had_hits) {
                for (auto& h : hs.get_left_hits()) {
                    if (!h.second.empty()) { unitig_ids.push_back(h.second.contigIdx_); }
                }
            }
        };

        auto get_hits_kallisto_psa = [&index, &kallisto_hits](
                                         const std::string& seq,
                                         std::vector<uint32_t>& unitig_ids) -> void {
            kallisto_hits.clear();
            match(seq, seq.length(), &index, kallisto_hits);
            if (!kallisto_hits.empty()) {
                for (auto& h : kallisto_hits) { unitig_ids.push_back(h.first.contigIdx_); }
            }
        };
        // Get the read group by which this thread will
        // communicate with the parser (*once per-thread*)
        auto rg = rparser.getReadGroup();
        while (rparser.refill(rg)) {
            // Here, rg will contain a chunk of read pairs we can process.
            for (auto const& record : rg) {
                switch (algo) {
                    case pseudoalignment_algorithm::SKIPPING:
                        get_hits_piscem_psa(record.seq, unitig_ids);
                        break;
                    case pseudoalignment_algorithm::SKIPPING_KALLISTO:
                        get_hits_kallisto_psa(record.seq, unitig_ids);
                        break;
                    default:
                        break;
                }

                num_reads += 1;
                index.intersect_unitigs(unitig_ids, colors);
                if (!colors.empty()) {
                    num_mapped_reads += 1;
                    ss << record.name << "\t" << colors.size();
                    for (auto c : colors) { ss << "\t" << c; }
                    ss << "\n";
                } else {
                    ss << record.name << "\t0\n";
                }
                colors.clear();
                unitig_ids.clear();

                buff_size += 1;
                if (num_reads > 0 and num_reads % 1000000 == 0) {
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
    } else {
        auto rg = rparser.getReadGroup();
        while (rparser.refill(rg)) {
            for (auto const& record : rg) {
                switch (algo) {
                    case pseudoalignment_algorithm::FULL_INTERSECTION:
                        index.pseudoalign_full_intersection(record.seq, colors);
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
                    ss << record.name << "\t" << colors.size();
                    for (auto c : colors) { ss << "\t" << c; }
                    ss << "\n";
                } else {
                    ss << record.name << "\t0\n";
                }
                num_reads += 1;
                colors.clear();
                if (num_reads > 0 and num_reads % 1000000 == 0) {
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
int pseudoalign(std::string const& index_filename, std::string const& query_filename,
                std::string const& output_filename, uint64_t num_threads, double threshold,
                pseudoalignment_algorithm algo) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");

    // if not a skipping variant and no threshold set, then set the algorithm
    if ((algo == pseudoalignment_algorithm::FULL_INTERSECTION) and
        (threshold != constants::invalid_threshold)) {
        algo = pseudoalignment_algorithm::THRESHOLD_UNION;
    }

    std::cerr << "query mode : " << to_string(algo, threshold) << "\n";

    if (((algo == pseudoalignment_algorithm::SKIPPING) or
         (algo == pseudoalignment_algorithm::SKIPPING_KALLISTO)) and
        !(index.get_dict().canonicalized())) {
        std::cout << "==> Warning: skipping is only supported for canonicalized indexes. <=="
                  << std::endl;
    }

    std::ifstream is(query_filename.c_str());
    if (!is.good()) {
        std::cerr << "error in opening the file '" + query_filename + "'" << std::endl;
        return 1;
    }

    essentials::logger("performing queries from file '" + query_filename + "'...");
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> t;
    t.start();

    std::atomic<uint64_t> num_mapped_reads{0};
    std::atomic<uint64_t> num_reads{0};

    auto query_filenames = std::vector<std::string>({query_filename});
    if (num_threads == 1) {
        num_threads += 1;
        essentials::logger(
            "1 thread was specified, but an additional thread will be allocated for parsing");
    }
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(query_filenames, num_threads,
                                                             num_threads - 1);

    rparser.start();
    std::vector<std::thread> workers;
    std::mutex iomut;
    std::mutex ofile_mut;

    std::ofstream out_file;
    out_file.open(output_filename, std::ios::out | std::ios::trunc);
    if (!out_file) {
        essentials::logger("could not open output file " + output_filename);
        return 1;
    }

    for (uint64_t i = 1; i != num_threads; ++i) {
        workers.push_back(std::thread([&index, &rparser, &num_reads, &num_mapped_reads, algo,
                                       threshold, &out_file, &iomut, &ofile_mut]() {
            do_map(index, rparser, num_reads, num_mapped_reads, algo, threshold, out_file, iomut,
                   ofile_mut);
        }));
    }

    for (auto& w : workers) { w.join(); }
    rparser.stop();

    t.stop();
    essentials::logger("DONE");

    std::cout << "mapped " << num_reads << " reads" << std::endl;
    std::cout << "elapsed = " << t.elapsed() << " millisec / ";
    std::cout << t.elapsed() / 1000 << " sec / ";
    std::cout << t.elapsed() / 1000 / 60 << " min / ";
    std::cout << (t.elapsed() * 1000) / num_reads << " musec/read" << std::endl;
    std::cout << "num_mapped_reads " << num_mapped_reads << "/" << num_reads << " ("
              << (num_mapped_reads * 100.0) / num_reads << "%)" << std::endl;

    return 0;
}

int pseudoalign(std::string const& index_filename, std::string const& query_filename,
                std::string const& output_filename, uint64_t num_threads, double threshold,
                pseudoalignment_algorithm algo, bool meta) {
    // if (meta) {
    return pseudoalign<meta_index_type>(index_filename, query_filename, output_filename,
                                        num_threads, threshold, algo);
    // }
    // return pseudoalign<index_type>(index_filename, query_filename, output_filename, num_threads,
    //                                threshold, algo);
}

int pseudoalign(int argc, char** argv) {
    std::string index_filename;
    std::string query_filename;
    std::string output_filename;
    uint64_t num_threads = 1;
    double threshold = constants::invalid_threshold;
    pseudoalignment_algorithm algo = pseudoalignment_algorithm::FULL_INTERSECTION;
    bool meta = false;

    CLI::App app{"Perform (color-only) pseudoalignment to a Fulgor index."};
    app.add_option("-i,--index", index_filename, "Fulgor index filename")
        ->required()
        ->check(CLI::ExistingFile);

    app.add_flag_callback(
        "--meta", [&meta]() { meta = true; }, "Specify if the Fulgor index is meta-colored.");

    app.add_option("-q,--query", query_filename,
                   "Query filename inf FASTA/FASTQ format (optionally gzipped)")
        ->required();
    app.add_option("-o,--output", output_filename, "File where output should be written")
        ->required();
    app.add_option("-t,--threads", num_threads, "Number of threads")->default_val(1);
    app.add_option("--threshold", threshold, "Threshold for threshold_union algorithm.")
        ->check(CLI::Range(0.0, 1.0));
    auto skip_opt = app.add_flag_callback(
        "--skipping", [&algo]() { algo = pseudoalignment_algorithm::SKIPPING; },
        "Enable the skipping heuristic in pseudoalignment");
    app.add_flag_callback(
           "--skipping-kallisto",
           [&algo]() { algo = pseudoalignment_algorithm::SKIPPING_KALLISTO; },
           "Enable the kallisto skipping heuristic in pseudoalignment")
        ->excludes(skip_opt);
    CLI11_PARSE(app, argc, argv);

    util::print_cmd(argc, argv);

    return pseudoalign(index_filename, query_filename, output_filename, num_threads, threshold,
                       algo, meta);
}

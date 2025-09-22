#include <iostream>
#include <fstream>
#include <sstream>

#include "src/ps_full_intersection.cpp"
#include "src/ps_threshold_union.cpp"
#include "external/sketch/include/sketch/hll.h"

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

template<typename Container>
struct container_hash{
    std::size_t operator()(Container const& vec) const {
        std::size_t seed = vec.size();
        for(auto x : vec) {
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

typedef std::unordered_map<std::vector<uint32_t>, std::vector<uint32_t>, container_hash<std::vector<uint32_t>>> cache_t;

template <typename FulgorIndex>
int find_partial_queries(FulgorIndex const& index, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser,
           std::atomic<uint64_t>& num_reads, std::atomic<uint64_t>& num_mapped_reads,
           std::mutex& iomut, std::vector<sketch::hll_t>& sketches, std::vector<uint8_t>& counts, const bool verbose)
{
    std::vector<uint32_t> sets_ids;  // result of pseudoalignment

    typename sketch::hll_t::HashType hasher;
    auto rg = rparser.getReadGroup();
    while (rparser.refill(rg)) {
        auto rg_offset_info = rg.chunk_frag_offset();
        uint32_t file_idx = rg_offset_info.file_idx;
        // TODO: This is a horrible hack right now! 
        // we encode the file of origin (if there is > 1 input file)
        // in the top 2 bits, and the read id in that file in 
        // the lower 30 bits. Do this more robustly. It doesn't 
        // make sense with > 4 inputs.
        uint32_t file_mask = (0x2 & file_idx) << 30;
        uint32_t read_id = rg_offset_info.frag_idx;

        for (auto const& record : rg) {
            uint32_t record_info = file_mask | read_id;
            std::string sequence = record.seq;
            if (sequence.length() < index.k()) continue;
            sets_ids.clear();
            std::vector<uint64_t> unitig_ids;

            { /* stream through */
                sshash::streaming_query<kmer_type, true> query(&index.get_k2u());
                query.reset();
                const uint64_t num_kmers = sequence.length() - index.k() + 1;
                for (uint64_t i = 0, prev_unitig_id = -1; i != num_kmers; ++i) {
                    char const* kmer = sequence.data() + i;
                    auto answer = query.lookup_advanced(kmer);
                    if (answer.kmer_id != sshash::constants::invalid_uint64) {  // kmer is positive
                        if (answer.contig_id != prev_unitig_id) {
                            unitig_ids.push_back(answer.contig_id);
                            prev_unitig_id = answer.contig_id;
                        }
                    }
                }
            }

            /* here we use it to hold the color set ids;
               in meta_intersect we use it to hold the partition ids */
            std::vector<uint32_t> tmp;

            /* deduplicate unitig_ids */
            std::sort(unitig_ids.begin(), unitig_ids.end());
            auto end_unitigs = std::unique(unitig_ids.begin(), unitig_ids.end());
            tmp.reserve(end_unitigs - unitig_ids.begin());
            for (auto it = unitig_ids.begin(); it != end_unitigs; ++it) {
                uint32_t unitig_id = *it;
                uint32_t color_set_id = index.u2c(unitig_id);
                tmp.push_back(color_set_id);
            }

            /* deduplicate color set ids */
            std::sort(tmp.begin(), tmp.end());
            auto end_tmp = std::unique(tmp.begin(), tmp.end());
            sets_ids.reserve(end_tmp - tmp.begin());
            if (end_tmp - tmp.begin() > 1){
                for (auto it = tmp.begin(); it != end_tmp; ++it) {
                    sets_ids.push_back(*it);
                    sketches[*it].add(hasher.hash(record_info));
                    counts[*it] += (counts[*it] < std::numeric_limits<uint8_t>::max());
                }
            }

            // store data
            sets_ids.clear();
            if (verbose and num_reads > 0 and num_reads % 1000000 == 0) {
                iomut.lock();
                std::cout << "mapped " << num_reads << " reads" << std::endl;
                iomut.unlock();
            }
            ++read_id;
        }
    }

    return 0;
}

template <typename FulgorIndex>
int batched_pseudoalign(FulgorIndex const& index, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser,
                     std::atomic<uint64_t>& num_reads, std::atomic<uint64_t>& num_mapped_reads,
                     pseudoalignment_algorithm algo, const double threshold, std::ofstream& out_file,
                     std::mutex& iomut, std::mutex& ofile_mut, bits::bit_vector const& s2m, bits::rank9 const& s2m_index,
                     cache_t& cache, std::vector<uint32_t> const& permutation, const bool verbose)  //
{
    std::vector<uint32_t> colors;  // result of pseudoalignment
    std::stringstream ss;
    uint64_t buff_size = 0;
    constexpr uint64_t buff_thresh = 50;
    uint64_t cached = 0;
    uint64_t saving = 0;
    uint64_t total = 0;

    auto rg = rparser.getReadGroup();
    while (rparser.refill(rg)) {
        for (auto const& record : rg) {
            std::string sequence = record.seq;
            if (sequence.length() < index.k()) continue;
            std::vector<uint64_t> unitig_ids;

            { /* stream through */
                sshash::streaming_query<kmer_type, true> query(&index.get_k2u());
                query.reset();
                const uint64_t num_kmers = sequence.length() - index.k() + 1;
                for (uint64_t i = 0, prev_unitig_id = -1; i != num_kmers; ++i) {
                    char const* kmer = sequence.data() + i;
                    auto answer = query.lookup_advanced(kmer);
                    if (answer.kmer_id != sshash::constants::invalid_uint64) {  // kmer is positive
                        if (answer.contig_id != prev_unitig_id) {
                            unitig_ids.push_back(answer.contig_id);
                            prev_unitig_id = answer.contig_id;
                        }
                    }
                }
            }

            /* here we use it to hold the color set ids;
               in meta_intersect we use it to hold the partition ids */
            std::vector<uint32_t> tmp;

            /* deduplicate unitig_ids */
            std::sort(unitig_ids.begin(), unitig_ids.end());
            auto end_unitigs = std::unique(unitig_ids.begin(), unitig_ids.end());
            tmp.reserve(end_unitigs - unitig_ids.begin());
            for (auto it = unitig_ids.begin(); it != end_unitigs; ++it) {
                uint32_t unitig_id = *it;
                uint32_t color_set_id = index.u2c(unitig_id);
                tmp.push_back(color_set_id);
            }

            /* deduplicate color set ids */
            std::sort(tmp.begin(), tmp.end());
            auto end_tmp = std::unique(tmp.begin(), tmp.end());
            if (end_tmp - tmp.begin() > 0){
                std::vector<uint32_t> batch;
                total += end_tmp - tmp.begin();
                
                std::vector<uint32_t> permuted;
                permuted.reserve(end_tmp - tmp.begin());
                for (auto it = tmp.begin(); it != end_tmp; ++it) {
                    permuted.push_back(permutation[*it]);
                }
                std::sort(permuted.begin(), permuted.end());

                uint32_t group_id = s2m_index.rank1(s2m, permuted[0]);
                for (auto id : permuted) {
                    uint32_t curr_group_id = s2m_index.rank1(s2m, id);
                    if (curr_group_id != group_id){
                        if (batch.size() > 1){
                            iomut.lock();
                            if (cache.count(batch) != 0){
                                cached++;
                                saving += batch.size() - 1;
                            } else {
                                cache.emplace(batch, std::vector<uint32_t>(1, 0));
                            }
                            iomut.unlock();
                        }
                        group_id = curr_group_id;
                        batch.clear();
                    }
                    batch.push_back(id);
                }
                num_mapped_reads += 1;
            }


            buff_size += 1;
            if (!tmp.empty()) {
                ss << record.name << '\t' << tmp.size();
                for (auto c : tmp) { ss << "\t" << c; }
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

    {
        iomut.lock();
        cout << " [#] " << cached << "/" << total << " cached intersections were used (" << saving << " lists were not iterated)" << endl;
        iomut.unlock();
    }

    return 0;
}

template <typename FulgorIndex>
int pseudoalign(FulgorIndex const& index, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser,
                std::atomic<uint64_t>& num_reads, std::atomic<uint64_t>& num_mapped_reads,
                pseudoalignment_algorithm algo, const double threshold, std::ofstream& out_file,
                std::mutex& iomut, std::mutex& ofile_mut, const bool verbose)  //
{
    std::vector<uint32_t> colors;  // result of pseudoalignment
    std::stringstream ss;
    uint64_t buff_size = 0;
    constexpr uint64_t buff_thresh = 50;

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
int pseudoalign(std::string const& index_filename, std::string const& query_filename,
                std::string const& output_filename, uint64_t num_threads, double threshold,
                pseudoalignment_algorithm ps_alg, const bool verbose, const bool batched) {
    FulgorIndex index;
    if (verbose) essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    if (verbose) essentials::logger("DONE");

    std::cerr << "query mode : " << (batched ? "BATCHED-" : "") << to_string(ps_alg, threshold) << "\n";

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

    const uint64_t p = 7; // 6 is sufficiently fast
    if (batched) {
        const uint64_t num_color_sets = index.num_color_sets();
        std::vector<sketch::hll_t> sketches (num_color_sets, sketch::hll_t(p));
        std::vector<uint8_t> set_counts(num_color_sets);

        for (uint64_t i = 1; i != num_threads; ++i) {

            workers.push_back(std::thread([&index, &rparser, &num_reads, &num_mapped_reads, 
                                           &iomut, &sketches, &set_counts, verbose]() {
                find_partial_queries(index, rparser, num_reads, num_mapped_reads, iomut, sketches, set_counts, verbose);
            }));
        }
        for (auto& w : workers) w.join();

        {
            const uint8_t sketch_threshold = 2;
            // std::vector<sketch::hll_t> sketches;
            // sketches.reserve(num_color_sets);
            for (uint64_t set_id = 0; set_id < num_color_sets; set_id++){
                if (set_counts[set_id] < sketch_threshold) {
                   sketches[set_id] = sketch::hll_t(p);
                }
            }

            std::ofstream out("sketches.bin", std::ios::binary);
            if (!out.is_open()) throw std::runtime_error("cannot open file");
            const uint64_t num_bytes = 1ULL << p;
            out.write(reinterpret_cast<char const*>(&num_bytes), 8);
            out.write(reinterpret_cast<char const*>(&num_color_sets), 8);
            for (auto const& x : sketches) {
                assert(x.m() == num_bytes);
                assert(x.m() == x.core().size());
                uint8_t const* data = x.data();
                out.write(reinterpret_cast<char const*>(data), num_bytes);
            }
            out.close();
        }
        
        essentials::logger("step 3. clustering sketches");

        std::ifstream in("sketches.bin", std::ios::binary);
        if (!in.is_open()) throw std::runtime_error("error in opening file");

        std::vector<kmeans::point> points;
        uint64_t num_bytes_per_point = 0;
        uint64_t num_points = 0;
        in.read(reinterpret_cast<char*>(&num_bytes_per_point), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&num_points), sizeof(uint64_t));
        cout << num_points << endl;
        points.resize(num_points, kmeans::point(num_bytes_per_point));
        for (auto& point : points) {
            in.read(reinterpret_cast<char*>(point.data()), num_bytes_per_point);
        }
        in.close();

        std::remove("sketches.bin");

        kmeans::clustering_parameters params;

        /* kmeans_divisive */
        constexpr float min_delta = 0.001;
        constexpr float max_iteration = 3;
        constexpr uint64_t min_cluster_size = 1;
        constexpr uint64_t seed = 0;
        params.set_min_delta(min_delta);
        params.set_max_iteration(max_iteration);
        params.set_min_cluster_size(min_cluster_size);
        params.set_random_seed(seed);
        params.set_num_threads(num_threads);
        auto clustering_data = kmeans::kmeans_divisive(points.begin(), points.end(), params);

        std::cout << "CLUSTERING DONE!" << std::endl; 

        uint64_t num_partitions = clustering_data.num_clusters;
        uint64_t max_partition_size = 0;
        std::vector<uint32_t> partition_size;
        std::vector<uint32_t> permutation;

        partition_size.resize(num_partitions + 1, 0);
        for (auto c : clustering_data.clusters) partition_size[c] += 1;

        /* take prefix sums */
        uint64_t val = 0;
        uint64_t bigger_that_one = 0;
        bits::bit_vector::builder bvb;
        bvb.resize(num_color_sets);
        for (auto& size : partition_size) {
            if (size > max_partition_size) max_partition_size = size;
            bigger_that_one += size > 1;

            uint64_t tmp = size;
            size = val;
            val += tmp;
            bvb.set(val-1);
        }
        bits::bit_vector set2meta;
        bvb.build(set2meta);
        bits::rank9 rank1_s2m;
        rank1_s2m.build(set2meta);
        assert(rank1_s2m.num_ones() == num_partitions);

        /* build permutation */
        auto counts = partition_size;  // copy
        permutation.resize(num_color_sets);
        assert(clustering_data.clusters.size() == index.num_color_sets());
        for (uint64_t i = 0; i != num_color_sets; ++i) {
            uint32_t cluster_id = clustering_data.clusters[i];
            permutation[i] = counts[cluster_id];
            counts[cluster_id] += 1;
        }

        std::cout << "Computed " << num_partitions << " partitions (" << bigger_that_one << " of which > 1)" << std::endl;

        essentials::logger("step 3. pseudoalignment");
        fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser2(query_filenames, num_threads,
                                                                 num_threads - 1);
        rparser2.start();
        cache_t cache;
        workers.clear();
        for (uint64_t i = 1; i != num_threads; ++i) {
            workers.push_back(std::thread([&index, &rparser2, &num_reads, &num_mapped_reads, ps_alg,
                                           threshold, &out_file, &iomut, &ofile_mut, &set2meta, &rank1_s2m, &cache, &permutation, verbose]() {
                batched_pseudoalign(index, rparser2, num_reads, num_mapped_reads, ps_alg, threshold, out_file,
                            iomut, ofile_mut, set2meta, rank1_s2m, cache, permutation, verbose);
            }));
        }

        for (auto& w : workers) w.join();
        rparser2.stop();
    } else {
        for (uint64_t i = 1; i != num_threads; ++i) {
            workers.push_back(std::thread([&index, &rparser, &num_reads, &num_mapped_reads, ps_alg,
                                           threshold, &out_file, &iomut, &ofile_mut, verbose]() {
                pseudoalign(index, rparser, num_reads, num_mapped_reads, ps_alg, threshold, out_file,
                            iomut, ofile_mut, verbose);
            }));
        }

        for (auto& w : workers) w.join();

    }

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
    parser.add("batched", "WARNING, DEVELOPMENT FEATURE. Process queries in batches", "--batched", false,
               true);
    parser.add("threshold",
               "Threshold for threshold_union algorithm. It must be a float in (0.0,1.0].", "-r",
               false);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto query_filename = parser.get<std::string>("query_filename");
    auto output_filename = parser.get<std::string>("output_filename");
    bool batched = parser.get<bool>("batched");

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

    if (sshash::util::ends_with(index_filename,
                                constants::meta_diff_colored_fulgor_filename_extension)) {
        return pseudoalign<meta_differential_index_type>(index_filename, query_filename,
                                                         output_filename, num_threads, threshold,
                                                         ps_alg, verbose, batched);
    } else if (sshash::util::ends_with(index_filename,
                                       constants::meta_colored_fulgor_filename_extension)) {
        return pseudoalign<meta_index_type>(index_filename, query_filename, output_filename,
                                            num_threads, threshold, ps_alg, verbose, batched);
    } else if (sshash::util::ends_with(index_filename,
                                       constants::diff_colored_fulgor_filename_extension)) {
        return pseudoalign<differential_index_type>(index_filename, query_filename, output_filename,
                                                    num_threads, threshold, ps_alg, verbose, batched);
    } else if (sshash::util::ends_with(index_filename, constants::fulgor_filename_extension)) {
        return pseudoalign<index_type>(index_filename, query_filename, output_filename, num_threads,
                                       threshold, ps_alg, verbose, batched);
    }

    std::cerr << "Wrong index filename supplied." << std::endl;

    return 1;
}

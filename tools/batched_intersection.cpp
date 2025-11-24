#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>

#include "../external/FQFeeder/include/FastxParser.hpp"
#include "../external/FQFeeder/include/blockingconcurrentqueue.h"

#include "../external/unordered_dense/include/ankerl/unordered_dense.h"
#include "include/util.hpp"

void preprocess_pref_suff(index_type const& index, frequent_map_t& frequent_map, std::vector<query_t>* sbatch) {
    const uint64_t ps_thresh = 2;

    std::vector<uint32_t> pkey;
    std::vector<uint32_t> skey;
    query_t prev = (*sbatch)[0];
    std::vector<uint32_t> colors_tmp;
    // use a stack to improve cached intersection efficiency?
    auto cache_partial_intersection = [&](std::vector<uint32_t> const& key) {
        if (frequent_map.find(key) != frequent_map.end()) return;

        colors_tmp.clear();
        intersect_color_ids(index, key, colors_tmp);
        frequent_map.emplace(key, colors_tmp);
    };

    for (uint64_t qid = 1; qid < sbatch->size(); ++qid) {
        auto& query = (*sbatch)[qid];
        if (query.cids.empty()) continue;
        pkey.clear();
        skey.clear();

        auto ps_info = find_max_prefix_suffix(prev.cids, query.cids);
        if (ps_info.prefix_len >= ps_thresh) {
            pkey.assign(query.cids.begin(), query.cids.begin() + ps_info.prefix_len);
            cache_partial_intersection(pkey);
        }
        if (ps_info.suffix_len >= ps_thresh) {
            skey.assign(query.cids.rbegin(), query.cids.rbegin() + ps_info.suffix_len);
            cache_partial_intersection(skey);
        }
        prev = query;
    }
}

std::pair<uint64_t, uint64_t> match_pref_suff(frequent_map_t& frequent_map, query_t& query,
                                              std::vector<color_info>& uncompressed_res) {
    auto fit = query.cids.begin();
    auto rit = query.cids.rbegin();
    color_info cached_pref(query.cids.end(), query.cids.end(), query.cids.end());
    color_info cached_suff(query.cids.end(), query.cids.end(), query.cids.end());
    std::vector<uint32_t> pkey;
    std::vector<uint32_t> skey;
    uint64_t pkey_size = 0;
    uint64_t skey_size = 0;

    while (fit < query.cids.end()) {
        pkey.push_back(*fit);
        skey.push_back(*rit);
        if (frequent_map.find(pkey) != frequent_map.end()) {
            cached_pref = {frequent_map[pkey].begin(), frequent_map[pkey].end(),
                          frequent_map[pkey].begin()};
            pkey_size = pkey.size();
        }
        if (frequent_map.find(skey) != frequent_map.end()) {
            cached_suff = {frequent_map[skey].begin(), frequent_map[skey].end(),
                          frequent_map[skey].begin()};
            skey_size = skey.size();
        }
        ++fit;
        ++rit;
    }

    if (cached_pref.begin != query.cids.end()) {
        uncompressed_res.push_back(cached_pref);
    }
    if (cached_suff.begin != query.cids.end()) {
        uncompressed_res.push_back(cached_suff);
    }

    return {pkey_size, skey_size};
}

void preprocess_intersection(index_type const& index, frequent_map_t& frequent_map, std::vector<query_t>* sbatch) {
    const uint64_t ps_thresh = 3;

    query_t prev = (*sbatch)[0];
    uint64_t prev_id = 0;
    std::vector<uint32_t> colors_tmp, cached_int;
    // use a stack to improve cached intersection efficiency?
    auto cache_partial_intersection = [&](std::vector<uint32_t> const& key) {
        if (frequent_map.find(key) != frequent_map.end()) return;

        colors_tmp.clear();
        intersect_color_ids(index, key, colors_tmp);
        frequent_map.emplace(key, colors_tmp);
    };

    std::vector<uint32_t> common_ids;
    std::vector<std::vector<uint32_t>> cache;

    for (uint64_t qid = 1; qid < sbatch->size(); ++qid) {
        query_t& query = (*sbatch)[qid];
        if (query.cids.empty()) continue;

        std::vector<color_info> iterators = {
            {prev.cids.begin(), prev.cids.end(), prev.cids.begin()},
            {query.cids.begin(), query.cids.end(), query.cids.begin()}};
        common_ids.clear();
        intersect_uncompressed(iterators, UINT32_MAX - 1, common_ids);
        if (common_ids.size() >= ps_thresh) {
            cache.push_back(common_ids);
            (*sbatch)[qid].cached_intersection = common_ids;
            if ((*sbatch)[prev_id].cached_intersection.size() < common_ids.size()) {
                (*sbatch)[prev_id].cached_intersection = common_ids;
            }
        }

        prev = query;
        prev_id = qid;
    }
    std::sort(cache.begin(), cache.end(), [](std::vector<uint32_t> const& a, std::vector<uint32_t> const& b) {
        return a.size() < b.size();
    });
    for (auto const& c : cache) {
        cache_partial_intersection(c);
    }
}

void match_intersection(frequent_map_t& frequent_map, query_t& query, std::vector<color_info>& uncompressed_res) {
    auto end_it = std::set_difference(query.cids.begin(), query.cids.end(),
                        query.cached_intersection.begin(), query.cached_intersection.end(),
                       query.cids.begin());
    query.cids.erase(end_it, query.cids.end());
    if (frequent_map.find(query.cached_intersection) != frequent_map.end()) {
        uncompressed_res.emplace_back(frequent_map[query.cached_intersection].begin(), frequent_map[query.cached_intersection].end(), frequent_map[query.cached_intersection].begin());
    }
}

void vertical_partition_queries(std::vector<query_t>* batch, std::vector<uint32_t>& cluster_sizes, std::vector<uint32_t>& permutation) {
    assert(cluster_sizes.empty());
    const uint64_t p = 8;
    const uint64_t num_bytes_per_point = 1ULL << p;
    std::unordered_map<uint32_t, sketch::hll_t> sketches; // map since set ids are more or less random
    typename sketch::hll_t::HashType hasher;

    for (auto& query: *batch) {
        for (auto it = query.cids.begin(); it != query.cids.end(); ++it) {
            if (sketches.find(*it) == sketches.end()) {
                sketches[*it] = sketch::hll_t(p);
            }
            sketches[*it].add(hasher.hash(query.id));
        }
    }
    const uint64_t num_points = sketches.size();
    if (num_points == 0) {
        cout << "wow"<< endl;
        return ;
    };

    std::vector<kmeans::point> points(num_points, kmeans::point(num_bytes_per_point));
    std::vector<uint32_t> ids;
    ids.reserve(num_points);
    auto points_it = points.begin();
    for (auto& [id, sketch] : sketches) {
        ids.push_back(id);
        const char* raw = reinterpret_cast<const char*>(sketch.data());
        std::memcpy(points_it->data(), raw, num_bytes_per_point);
        ++points_it;
    }
    assert(ids.size() == num_points);

    /* kmeans_divisive */
    kmeans::clustering_parameters params;
    constexpr float min_delta = 0.001;
    constexpr float max_iteration = 10;
    constexpr uint64_t min_cluster_size = 1;
    constexpr uint64_t seed = 0;
    params.set_min_delta(min_delta);
    params.set_max_iteration(max_iteration);
    params.set_min_cluster_size(min_cluster_size);
    params.set_random_seed(seed);
    params.set_num_threads(1);

    auto clustering_data = kmeans::kmeans_divisive(points.begin(), points.end(), params);
    cluster_sizes.resize(clustering_data.num_clusters + 1);
    for (auto cluster_id : clustering_data.clusters) {
        cluster_sizes[cluster_id]++;
    }
    uint64_t val = 0;
    for (auto& size : cluster_sizes) {
        uint64_t tmp = size;
        size = val;
        val += tmp;
    }

    std::unordered_map<uint32_t, uint32_t> id_to_meta;
    for (uint64_t i = 0; i < num_points; ++i) {
        id_to_meta[ids[i]] = cluster_sizes[clustering_data.clusters[i]]++;
    }
    assert(id_to_meta.size() == num_points);
    for (auto& query: *batch) {
        for (uint64_t i = 0; i < query.cids.size(); ++i) {
            query.cids[i] = id_to_meta[query.cids[i]];
        }
        std::sort(query.cids.begin(), query.cids.end());
    }
    permutation.resize(num_points);
    for (auto& [original_id, meta_id] : id_to_meta) {
        permutation[meta_id] = original_id;
    }
}

using namespace fulgor;

void process_lines_batched(index_type& index,
                        moodycamel::BlockingConcurrentQueue<std::vector<query_t>*>& q,
                        std::atomic<bool>& done, std::atomic<uint32_t>& in_flight,
                        std::mutex& ofile_mut, std::ofstream& out_file,
                        std::atomic<uint64_t>& num_mapped_reads) {
    std::stringstream ss;
    int32_t buff_size{0};
    constexpr int32_t buff_thresh{50};
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> compressed_timer;
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> uncompressed_timer;
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> preprocess_timer;
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> mapping_timer;
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> working_timer;
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> writing_timer;
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> other_timer;
    uint64_t num_compressed_int = 0, num_uncompressed_int = 0;

    std::vector<uint32_t> colors;
    std::vector<uint32_t> colors_tmp;
    std::vector<color_info> uncompressed_res;
    std::vector<query_t>* sbatch = nullptr;
    stringstream intersection_ss;

    while (q.try_dequeue(sbatch) or !done or (in_flight > 0)) {
        if (sbatch != nullptr) {
            frequent_map_t frequent_map;
            working_timer.start();
            in_flight -= 1;

            preprocess_timer.start();
            // preprocess_pref_suff(index, frequent_map, sbatch);
            // preprocess_intersection(index, frequent_map, sbatch);
            std::vector<uint32_t> cluster_begins, permutation;
            vertical_partition_queries(sbatch, cluster_begins, permutation);
            preprocess_timer.stop();
            uint32_t query_id;

            for (auto& query : *sbatch) {
                assert(query_id != query.id);
                query_id = query.id;
                if (!query.cids.empty()) {
                    // interesect the color ids to get the colors
                    colors.clear();
                    colors_tmp.clear();
                    uncompressed_res.clear();

                    mapping_timer.start();
                    // auto [pkey_size, skey_size] =  match_pref_suff(frequent_map, query, uncompressed_res);
                    // std::vector<uint32_t> uncached_cids = {query.cids.begin()+pkey_size, max(query.cids.end()-skey_size, query.cids.begin()+pkey_size)};
                    // match_intersection(frequent_map, query, uncompressed_res);
                    std::vector<uint32_t> partial_set;
                    std::vector<std::vector<uint32_t>> partial_sets;
                    uint64_t curr_partition = 0;
                    while (query.cids.front() >= cluster_begins[curr_partition+1]) {
                        ++curr_partition;
                    }
                    for (auto c : query.cids) {
                        if (c < cluster_begins[curr_partition + 1]) {
                            partial_set.push_back(permutation[c]); // insert original id
                        } else {
                            assert(!partial_set.empty());
                            if (frequent_map.find(partial_set) == frequent_map.end()) {
                                compressed_timer.start();
                                num_compressed_int++;
                                intersect_color_ids(index, partial_set, colors_tmp);
                                compressed_timer.stop();
                                frequent_map.emplace(partial_set, colors_tmp);
                                partial_sets.push_back(partial_set);
                            }
                            uncompressed_res.emplace_back(frequent_map[partial_set]);

                            partial_set.clear();
                            partial_set.push_back(c);
                            while (c >= cluster_begins[curr_partition + 1]) {
                                ++curr_partition;
                            }
                            colors_tmp.clear();
                        }
                    }
                    if (frequent_map.find(partial_set) == frequent_map.end()) {
                        compressed_timer.start();
                        num_compressed_int++;
                        intersect_color_ids(index, partial_set, colors_tmp);
                        compressed_timer.stop();
                        frequent_map.emplace(partial_set, colors_tmp);
                        partial_sets.push_back(partial_set);
                    }
                    uncompressed_timer.start();
                    num_uncompressed_int++;
                    uncompressed_res.emplace_back(frequent_map[partial_set]);
                    uncompressed_timer.stop();
                    colors_tmp.clear();
                    partial_set.clear();

                    mapping_timer.stop();

                    intersect_uncompressed(uncompressed_res, index.num_colors(), colors);

                    /*
                    if (!uncompressed_res.empty()){
                        colors_tmp.clear();
                        if (!query.cids.empty()) {
                            compressed_timer.start();
                            num_compressed_int++;
                            intersect_color_ids(index, query.cids, colors_tmp);
                            compressed_timer.stop();
                        }
                        uncompressed_timer.start();
                        num_uncompressed_int++;
                        intersect_uncompressed(uncompressed_res, index.num_colors(), colors);
                        uncompressed_timer.stop();
                    } else {
                        compressed_timer.start();
                        num_compressed_int++;
                        intersect_color_ids(index, query.cids, colors);
                        compressed_timer.stop();
                    }
                }
                */
                }

                other_timer.start();
                std::string tmpstr;
                if (!colors.empty()) {
                    ss << query_id << "\t" << colors.size();
                    //for (auto c : colors) {}
                    util::vec_to_tsv(colors, tmpstr);
                    ss << "\t" << tmpstr;
                    ss << "\n";
                    num_mapped_reads++;
                } else {
                    // num_mapped_reads -= 1;
                    ss << query_id << "\t0\n";
                }
                buff_size += 1;
                other_timer.stop();


                if (buff_size > buff_thresh) {
                    writing_timer.start();
                    std::string outs = ss.str();
                    ss.str("");
                    ofile_mut.lock();
                    out_file.write(outs.data(), outs.size());
                    ofile_mut.unlock();
                    buff_size = 0;
                    writing_timer.stop();
                }
            }

            delete sbatch;
            sbatch = nullptr;
            working_timer.stop();
        }
    }

    // dump anything left in the buffer
    if (buff_size > 0) {
        writing_timer.start();
        std::string outs = ss.str();
        ss.str("");
        ofile_mut.lock();
        out_file.write(outs.data(), outs.size());
        ofile_mut.unlock();
        buff_size = 0;
        writing_timer.stop();
    }
    // std::cerr << "(thread_index : " << thread_index << ") total_aln: " << total_aln << ",
    // skipped_aln: " << skipped_aln << ", batch_ctr: " << batch_ctr << "\n";
    ofile_mut.lock();
    cout << "Thread worked for: " << working_timer.elapsed()/1000 << "ms, ";
    cout << "Other stuff took: " << other_timer.elapsed()/1000 << "ms, ";
    cout << "Compressed intersections took: " << compressed_timer.elapsed()/1000 << "ms ("<< compressed_timer.elapsed() / num_compressed_int << " musec per intersection), ";
    cout << "Uncompressed intersections took: " << uncompressed_timer.elapsed()/1000 << "ms ("<< uncompressed_timer.elapsed() / num_uncompressed_int << " musec per intersection), ";
    cout << "Writing took: " << writing_timer.elapsed()/1000 << "ms, ";
    cout << "Preprocessing took: " << preprocess_timer.elapsed()/1000 << "ms, ";
    cout << "Mapping took: " << mapping_timer.elapsed()/1000 << "ms, ";
    cout << endl;
    ofile_mut.unlock();
}

void do_intersection_batched(index_type& index, size_t num_threads_in,
                             const std::string& query_filename, std::ofstream& output,
                             std::atomic<uint64_t>& num_mapped_reads) {
    std::ifstream ifile(query_filename, std::ios::binary);

    std::atomic<bool> done{false};
    std::atomic<uint32_t> in_flight{0};
    const uint64_t batch_thresh = 1000;

    moodycamel::BlockingConcurrentQueue<std::vector<query_t>*> q(3 * num_threads_in);

    std::thread producer([&q, &ifile, &done, &in_flight]() -> void {
        auto* batch = new std::vector<query_t>();
        batch->reserve(batch_thresh);
        uint32_t list_len = 0;
        size_t nbatches = 0;

        while (ifile.read(reinterpret_cast<char*>(&list_len), sizeof(list_len))) {
            if (list_len > 0) {
                // if (list_len > 100) { std::cerr << "list_len " << list_len << "\n"; }
                uint32_t query_id;
                ifile.read(reinterpret_cast<char*>(&query_id), sizeof(query_id));
                query_t query(query_id, list_len -1); // -1 since the first one is the id
                ifile.read(reinterpret_cast<char*>(query.cids.data()),
                           (list_len - 1) * sizeof(list_len));

                // if the current batch is big enough that we want to push it
                // make sure the current element isn't a duplicate (i.e. size 1)
                // and push the batch *before* we add the next element, which
                // may be the start of a new duplicate run.
                if ((batch->size() >= batch_thresh) and (!query.cids.empty())) {
                    while (!q.try_enqueue(batch)) {}
                    batch = new std::vector<query_t>();
                    in_flight += 1;
                    ++nbatches;
                }
                // always push back the current element. If it's in the middle of
                // a duplicate run, we keep growing the previous batch, otherwise
                // this is the first element of the enw batch.
                batch->push_back(query);
            } else {
                std::cerr << "should not happen\n";
            }
        }

        if (!batch->empty()) {
            q.enqueue(batch);
            in_flight += 1;
        }
        done = true;
    });

    std::mutex ofmut;
    std::vector<std::thread> workers;
    workers.reserve(num_threads_in);
    num_mapped_reads = 0;
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> timer;
    timer.start();
    for (size_t i = 0; i < num_threads_in; ++i) {
        workers.push_back(std::thread(
            [&index, &q, &done, &in_flight, &ofmut, &output, &num_mapped_reads]() -> void {
                process_lines_batched(index, q, done, in_flight, ofmut, output, num_mapped_reads); //, i, num_mapped_reads);
            }));
    }
    producer.join();
    for (auto& t : workers) { t.join(); }

    timer.stop();

    std::cout << "Only intersection took = " << timer.elapsed() << " millisec / ";
    std::cout << timer.elapsed() / 1000 << " sec / ";
    std::cout << timer.elapsed() / 1000 / 60 << " min / " << std::endl;
}

void sort_file_by_similarity(const std::string& tmp_outname,
                                             std::ofstream& output, const uint64_t num_threads_in) {
    std::vector<std::vector<uint32_t>> queries;
    std::ifstream ifile(tmp_outname, std::ios::binary);

    std::vector<uint32_t> vals;
    uint32_t read_num = 0;
    while (ifile.read(reinterpret_cast<char*>(&read_num), sizeof(read_num))) {
        uint32_t num_colors = 0;
        ifile.read(reinterpret_cast<char*>(&num_colors), sizeof(num_colors));

        vals.resize(num_colors + 1);
        vals[0] = read_num;
        ifile.read(reinterpret_cast<char*>(&vals[1]), num_colors * sizeof(num_colors));
        if (vals.size() > 1) {
            queries.push_back(vals);
            vals.clear();
        } else {
            // just write out the unmapped reads here
            output << read_num << "\t0\n";
        }
    }

    if (queries.empty()) {
        return;
    }

    essentials::logger("Begin sort");

    std::sort(queries.begin(), queries.end(),
              [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) -> bool {
                  return std::lexicographical_compare(a.begin() + 1, a.end(), b.begin() + 1,
                                                      b.end());
              });

    essentials::logger("Done sort");
    // probably should deduplicate here in some way, to remove noise for the clustering
    // store a vector of vectors of query ids, each inner vector is a set of equal queries
    // when clustering, compute the sketch for only one

    ifile.close();

    const uint64_t num_queries = queries.size(), batch_size = 250000;
    const uint64_t p = 8;
    uint64_t curr_batch = 0;
    typename sketch::hll_t::HashType hasher;
    sketch::hll_t sketch(p);
    const uint64_t num_bytes_per_point = 1ULL << p;
    std::vector<kmeans::point> points(batch_size, kmeans::point(num_bytes_per_point));

    /* kmeans_divisive */
    kmeans::clustering_parameters params;
    constexpr float min_delta = 0.001; // TODO: try bigger min_delta
    constexpr float max_iteration = 10;
    constexpr uint64_t min_cluster_size = 1;
    constexpr uint64_t seed = 0;
    params.set_min_delta(min_delta);
    params.set_max_iteration(max_iteration);
    params.set_min_cluster_size(min_cluster_size);
    params.set_random_seed(seed);
    params.set_num_threads(num_threads_in);

    std::ofstream ofile(tmp_outname, std::ios::trunc | std::ios::binary);

    essentials::logger("Begin batch clustering");
    while (curr_batch * batch_size < num_queries) {
        stringstream ss;
        ss << "Clustering batch " << curr_batch + 1 << " / " << num_queries/batch_size + 1;
        essentials::logger(ss.str());

        const uint64_t batch_start = curr_batch * batch_size;
        const uint64_t batch_end = std::min((curr_batch + 1) * batch_size, num_queries);
        points.resize(batch_end - batch_start);
        for (uint64_t i = batch_start; i < batch_end; ++i) {
            auto& cids = queries[i];
            sketch.clear();
            for (auto it = cids.begin()+1; it != cids.end(); ++it) {
                sketch.add(hasher.hash(*it));
            }
            const char* raw = reinterpret_cast<const char*>(sketch.data());
            std::memcpy(points[i-batch_start].data(), raw, num_bytes_per_point);
        }
        auto clustering_data = kmeans::kmeans_divisive(points.begin(), points.end(), params);
        std::vector<uint32_t> cluster_sizes(clustering_data.num_clusters + 1);
        {
            for (auto cluster_id : clustering_data.clusters) {
                cluster_sizes[cluster_id]++;
            }
            uint64_t val = 0;
            for (auto& size : cluster_sizes) {
                uint64_t tmp = size;
                size = val;
                val += tmp;
            }
        }

        // build permutation
        assert(clustering_data.clusters.size() == batch_end - batch_start);
        std::vector<uint32_t> permutation(batch_end - batch_start);
        std::vector<uint32_t> prev ;
        for (uint32_t i = 0; i < batch_end - batch_start; ++i) {
            uint32_t cluster_id = clustering_data.clusters[i];
            permutation[cluster_sizes[cluster_id]++] = i;
        }
        for (auto query_id : permutation) {
            query_id += batch_start;
            auto& query = queries[query_id];
            if (prev.size() == query.size() &&
                std::equal(prev.begin() + 1, prev.end(), query.begin() + 1)) {
                query.resize(1);
            } else {
                prev = query;
            }
            uint32_t s = query.size();
            ofile.write(reinterpret_cast<char*>(&s), sizeof(s));
            ofile.write(reinterpret_cast<char*>(query.data()), sizeof(s) * s);

        }

        ++curr_batch;
    }
    essentials::logger("Done batch clustering");

    ofile.close();
}

int intersect_colors_batched(int argc, char** argv) {
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

    index_type index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");

    // first map to get the color lists
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

    char tmp_filename[] = "tmp_color_listXXXXXX";
    mkstemp(tmp_filename);

    std::string tmp_outname(tmp_filename);
    essentials::logger("writing temporary color lists to " + tmp_outname);
    std::ofstream tmp_file;
    tmp_file.open(tmp_outname, std::ios::out | std::ios::trunc);
    if (!tmp_file) {
        essentials::logger("could not open output file " + tmp_outname);
        return 1;
    }

    for (size_t i = 1; i < num_threads; ++i) {
        workers.push_back(std::thread([&index, &rparser, &num_reads, &num_mapped_reads, &tmp_file,
                                       &iomut, &ofile_mut]() {
            do_color_map(index, rparser, num_reads, num_mapped_reads, tmp_file, iomut, ofile_mut);
        }));
    }

    for (auto& w : workers) { w.join(); }
    rparser.stop();

    t.stop();
    essentials::logger("DONE");
    tmp_file.close();

    std::cout << "mapped " << num_reads << " reads" << std::endl;
    std::cout << "elapsed = " << t.elapsed() << " millisec / ";
    std::cout << t.elapsed() / 1000 << " sec / ";
    std::cout << t.elapsed() / 1000 / 60 << " min / ";
    std::cout << (t.elapsed() * 1000) / num_reads << " musec/read" << std::endl;
    std::cout << "num_mapped_reads " << num_mapped_reads << "/" << num_reads << " ("
              << (num_mapped_reads * 100.0) / num_reads << "%)" << std::endl;

    std::ofstream output(output_filename);
    sort_file_by_similarity(tmp_outname, output, num_threads);

    do_intersection_batched(index, num_threads, tmp_outname, output, num_mapped_reads);
    std::remove(tmp_outname.c_str());

    std::cout << "num_mapped_reads " << num_mapped_reads << "/" << num_reads << " ("
              << (num_mapped_reads * 100.0) / num_reads << "%)" << std::endl;
    return 0;
}

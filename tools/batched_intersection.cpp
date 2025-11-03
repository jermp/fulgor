#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>

#include "../external/FQFeeder/include/FastxParser.hpp"
#include "../external/FQFeeder/include/blockingconcurrentqueue.h"

#include "../external/unordered_dense/include/ankerl/unordered_dense.h"

/* Credits to Jarno Alanko for the suggestion. */
void manual_conversion(std::vector<uint32_t> const& vec, std::string& s) {
    s.clear();
    s.reserve(vec.size() * 12);
    char buffer[32];
    buffer[31] = '\t';
    uint32_t tmp;
    for (uint32_t x : vec) {
        int len = 0;
        do {
            tmp = x / 10;
            buffer[30 - len++] = '0' + (x - tmp * 10);
            x = tmp;
        } while (x > 0);
        s.append(buffer + 31 - len, len + 1);
    }
    s.pop_back();
}

using namespace fulgor;

void process_lines_meta(index_type& index,
                        moodycamel::BlockingConcurrentQueue<std::vector<query_t>*>& q,
                        std::atomic<bool>& done, std::atomic<uint32_t>& in_flight,
                        std::mutex& ofile_mut, std::ofstream& out_file) {
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

    std::vector<uint32_t> pkey;
    std::vector<uint32_t> skey;
    std::vector<uint32_t> colors;
    std::vector<uint32_t> colors_tmp;
    std::vector<color_info> uncompressed_res;
    std::vector<query_t>* sbatch = nullptr;
    const uint64_t ps_thresh = 3;
    stringstream intersection_ss;


    while (q.try_dequeue(sbatch) or !done or (in_flight > 0)) {
        if (sbatch != nullptr) {
            frequent_map_t frequent_map;
            working_timer.start();
            in_flight -= 1;

            query_t prev = (*sbatch)[0];
            // use a stack to improve cached intersection efficiency?
            auto cache_partial_intersection = [&](std::vector<uint32_t> const& key) {
                if (frequent_map.find(key) != frequent_map.end()) return;

                colors_tmp.clear();
                intersect_color_ids(index, key, colors_tmp);
                frequent_map.emplace(key, colors_tmp);
            };

            preprocess_timer.start();
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
            preprocess_timer.stop();
            uint32_t query_id;

            for (auto& query : *sbatch) {
                // if the only thing in the cids vector is
                // 1 element (the read_id), then the color
                // output should be exactly the same as the
                // last call to index.instersect_color_ids
                // and so we don't recompute it here.
                assert(query_id != query.id);
                query_id = query.id;
                if (!query.cids.empty()) {
                    // interesect the color ids to get the colors
                    colors.clear();
                    uncompressed_res.clear();
                    pkey.clear();
                    skey.clear();

                    auto fit = query.cids.begin();
                    auto rit = query.cids.rbegin();
                    uint64_t pkey_size = 0;
                    uint64_t skey_size = 0;
                    color_info cached_pref(query.cids.end(), query.cids.end(), query.cids.end());
                    color_info cached_suff(query.cids.end(), query.cids.end(), query.cids.end());

                    mapping_timer.start();
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
                    mapping_timer.stop();

                    if (cached_pref.begin != query.cids.end()) {
                        uncompressed_res.push_back(cached_pref);
                        //cached_size += pkey_size;
                        //cached++;
                    }
                    if (cached_suff.begin != query.cids.end()) {
                        uncompressed_res.push_back(cached_suff);
                        //cached_size += pkey_size;
                        //cached++;
                    }
                    std::vector<uint32_t> uncached_cids = {query.cids.begin()+pkey_size, max(query.cids.end()-skey_size, query.cids.begin()+pkey_size)};

                    if (!uncompressed_res.empty()){
                        colors_tmp.clear();
                        if (!uncached_cids.empty()) {
                            compressed_timer.start();
                            intersect_color_ids(index, uncached_cids, colors_tmp); // if all cids are intersected, the result is correct
                            uncompressed_res.emplace_back(colors_tmp.begin(), colors_tmp.end(),
                                                          colors_tmp.begin());
                            compressed_timer.stop();
                        }
                        uncompressed_timer.start();
                        intersect_uncompressed(uncompressed_res, index.num_colors(), colors);
                        uncompressed_timer.stop();
                    } else {
                        compressed_timer.start();
                        intersect_color_ids(index, query.cids, colors);
                        compressed_timer.stop();
                    }
                }

                other_timer.start();
                std::string tmpstr;
                if (!colors.empty()) {
                    ss << query_id << "\t" << colors.size();
                    //for (auto c : colors) {}
                    manual_conversion(colors, tmpstr);
                    ss << "\t" << tmpstr;
                    ss << "\n";
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
            //tmpss << cached << " " << cached_size << endl;
            //cout << tmpss.str() << flush;

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
    cout << "Compressed intersections took: " << compressed_timer.elapsed()/1000 << "ms, ";
    cout << "Uncompressed intersections took: " << uncompressed_timer.elapsed()/1000 << "ms, ";
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
    const uint64_t batch_thresh = 10000;

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
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> timer;
    timer.start();
    for (size_t i = 0; i < num_threads_in; ++i) {
        workers.push_back(std::thread(
            [&index, &q, &done, &in_flight, &ofmut, &output]() -> void {
                process_lines_meta(index, q, done, in_flight, ofmut, output); //, i, num_mapped_reads);
            }));
    }
    producer.join();
    for (auto& t : workers) { t.join(); }

    timer.stop();
    std::cout << "Only intersection took = " << timer.elapsed() << " millisec / ";
    std::cout << timer.elapsed() / 1000 << " sec / ";
    std::cout << timer.elapsed() / 1000 / 60 << " min / " << std::endl;
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
    auto sorted = sort_file(tmp_outname, output);

    do_intersection_batched(index, num_threads, tmp_outname, output, num_mapped_reads);
    std::remove(tmp_outname.c_str());

    std::cout << "num_mapped_reads " << num_mapped_reads << "/" << num_reads << " ("
              << (num_mapped_reads * 100.0) / num_reads << "%)" << std::endl;
    return 0;
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>

#include "../external/FQFeeder/include/FastxParser.hpp"
#include "../external/FQFeeder/include/blockingconcurrentqueue.h"

#include "../external/unordered_dense/include/ankerl/unordered_dense.h"


using namespace fulgor;

void process_lines_meta(index_type& index,
                        moodycamel::BlockingConcurrentQueue<std::vector<query_t>*>& q,
                        std::atomic<bool>& done, std::atomic<uint32_t>& in_flight,
                        std::mutex& ofile_mut, std::ofstream& out_file) {
    std::stringstream ss;
    int32_t buff_size{0};
    constexpr int32_t buff_thresh{50};

    frequent_map_t frequent_map;
    std::vector<uint32_t> pkey;
    std::vector<uint32_t> skey;
    std::vector<uint32_t> colors;
    std::vector<uint32_t> colors_tmp;
    std::vector<color_info> uncompressed_res;
    std::vector<query_t>* sbatch = nullptr;
    const uint64_t ps_thresh = 3;

    while (q.try_dequeue(sbatch) or !done or (in_flight > 0)) {
        if (sbatch != nullptr) {
            in_flight -= 1;

            query_t prev = (*sbatch)[0];
            // use a stack to improve cached intersection efficiency?
            uint64_t cached = 0;
            uint64_t cached_size = 0;
            auto cache_partial_intersection = [&](std::vector<uint32_t> const& key) {
                cached++;
                cached_size += key.size();
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
                    // cache_partial_intersection(skey);
                }

                prev = query;
            }
            stringstream tmpss;
            tmpss << " > " << cached << " " << cached_size << endl;
            cout << tmpss.str() << flush;
            tmpss.str("");
            cached_size = 0;
            cached = 0;
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

                    auto it = query.cids.begin();
                    uint64_t pkey_size = 0;
                    color_info cached_res(query.cids.end(), query.cids.end(), query.cids.end());

                    while (it < query.cids.end()) {
                        pkey.push_back(*it);
                        if (frequent_map.find(pkey) != frequent_map.end()) {
                            cached_res = {frequent_map[pkey].begin(), frequent_map[pkey].end(),
                                          frequent_map[pkey].begin()};
                            pkey_size = pkey.size();
                        }
                        ++it;
                    }

                    if (cached_res.begin != query.cids.end()) {
                        uncompressed_res.push_back(cached_res);
                        std::vector<uint32_t> tmp_cids = {query.cids.begin()+pkey_size, query.cids.end()};
                        //cached_size += pkey_size;
                        //cached++;

                        colors_tmp.clear();
                        if (!tmp_cids.empty()) {
                            intersect_color_ids(index, tmp_cids, colors_tmp); // if all cids are intersected, the result is correct
                            uncompressed_res.emplace_back(colors_tmp.begin(), colors_tmp.end(),
                                                          colors_tmp.begin());
                        }
                        intersect_uncompressed(uncompressed_res, index.num_colors(), colors);
                    } else {
                        intersect_color_ids(index, query.cids, colors);
                    }
                }

                if (!colors.empty()) {
                    ss << query_id << "\t" << colors.size();
                    for (auto c : colors) { ss << "\t" << c; }
                    ss << "\n";
                } else {
                    // num_mapped_reads -= 1;
                    ss << query_id << "\t0\n";
                }
                buff_size += 1;

                if (buff_size > buff_thresh) {
                    std::string outs = ss.str();
                    ss.str("");
                    ofile_mut.lock();
                    out_file.write(outs.data(), outs.size());
                    ofile_mut.unlock();
                    buff_size = 0;
                }
            }
            //tmpss << cached << " " << cached_size << endl;
            //cout << tmpss.str() << flush;

            delete sbatch;
            sbatch = nullptr;
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
    // std::cerr << "(thread_index : " << thread_index << ") total_aln: " << total_aln << ",
    // skipped_aln: " << skipped_aln << ", batch_ctr: " << batch_ctr << "\n";
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
                if ((batch->size() >= 10000) and (!query.cids.empty())) {
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
            [&index, &q, &done, &in_flight, &ofmut, &output, &num_mapped_reads, i]() -> void {
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

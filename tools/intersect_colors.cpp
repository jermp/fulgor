#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>

// #include "../external/sshash/include/gz/zip_stream.hpp"
#include "../external/FQFeeder/include/FastxParser.hpp"
#include "../external/FQFeeder/include/blockingconcurrentqueue.h"

#include "../external/unordered_dense/include/ankerl/unordered_dense.h"

using namespace fulgor;

struct color_info {
    color_info(std::vector<uint32_t>::iterator begin, std::vector<uint32_t>::iterator end, std::vector<uint32_t>::iterator curr) 
        : begin(begin), end(end), curr(curr) {}
    std::vector<uint32_t>::iterator begin;
    std::vector<uint32_t>::iterator end;
    std::vector<uint32_t>::iterator curr;
    size_t size() const { return std::distance(begin, end); }
    bool is_exhausted() const { return curr >= end; }
};

struct custom_vec_hash {
    using is_avalanching = void;

    [[nodiscard]] auto operator()(std::vector<uint32_t> const& f) const noexcept -> uint64_t {
        return ankerl::unordered_dense::detail::wyhash::hash(f.data(), sizeof(uint32_t) * f.size());
    }
};

using frequent_map_t =
    ankerl::unordered_dense::map<std::vector<uint32_t>, std::vector<uint32_t>, custom_vec_hash>;

struct meta_sets {
    bits::bit_vector s2m;
    bits::rank9 s2m_index;
    std::vector<uint32_t> permutation;
    cache_t cache;
};

void next_geq(color_info& info, uint32_t tgt) {
    info.curr = std::lower_bound(info.curr, info.end, tgt);
}

template <typename ColorSets>
void intersect_color_ids(fulgor::index<ColorSets> const& index,
                         const std::vector<uint32_t>& color_ids, std::vector<uint32_t>& colors) {
    std::vector<typename ColorSets::iterator_type> iterators;
    iterators.reserve(color_ids.size()-1);
    for (auto it = color_ids.begin() + 1; it != color_ids.end(); ++it) {
        uint64_t color_set_id = *it;
        auto fwd_it = index.color_set(color_set_id);
        iterators.push_back(fwd_it);
    }
    fulgor::next_geq_intersect(iterators.begin(), iterators.end(), colors, index.num_colors());
}

void intersect_uncompressed(std::vector<color_info>& iterators, uint32_t num_colors,
                            std::vector<uint32_t>& colors) {
    assert(colors.empty());

    if (iterators.empty()) return;

    std::sort(iterators.begin(), iterators.end(),
              [](auto const& x, auto const& y) { return x.size() < y.size(); });

    if (iterators[0].is_exhausted()) {
        // std::cerr << "shouldn't happen";
        return;
    }

    uint32_t candidate = *(iterators[0].curr);
    uint64_t i = 1;
    while (candidate < num_colors) {
        for (; i < iterators.size(); ++i) {
            next_geq(iterators[i], candidate);
            uint32_t val = iterators[i].is_exhausted() ? num_colors : *(iterators[i].curr);
            if (val != candidate) {
                candidate = val;
                i = 0;
                break;
            }
        }
        if (i == iterators.size()) {
            colors.push_back(candidate);
            iterators[0].curr++;
            candidate = iterators[0].is_exhausted() ? num_colors : *(iterators[0].curr);
            i = 1;
        }
    }
}

void analyze_batch(std::vector<std::vector<uint32_t>>* batch, index_type& index,
                   frequent_map_t& frequent_map) {
    ankerl::unordered_dense::map<std::vector<uint32_t>, uint32_t, custom_vec_hash> count_map;
    std::vector<uint32_t> key;
    size_t k = 3;
    // uint32_t max_elem = 0;
    for (auto& cids : *batch) {
        if (cids.size() < k) { continue; }
        for (size_t i = 0; i < cids.size() - k; i++) {
            key.assign(cids.begin() + i, cids.begin() + i + k);
            auto& val = count_map[key];
            val += 1;
            // max_elem = std::max(val, max_elem);
        }
    }

    std::vector<uint32_t> colors;
    uint32_t thresh = 2;  // std::max(static_cast<uint32_t>(2), static_cast<uint32_t>(0.1 * max_elem) + 1);
    for (auto const& [key, val] : count_map) {
        auto kv = frequent_map.find(key);
        if (kv != frequent_map.end() and (val >= thresh)) {
            intersect_color_ids(index, key, colors);
            frequent_map[key] = colors;
            colors.clear();
        }
    }
}

// extract colors that are in the freuqent_map from
// the colors vector
std::vector<color_info> extract_frequent_colors(std::vector<uint32_t>& cids,
                                                frequent_map_t& frequent_map) {
    size_t k = 3;
    if (cids.size() < k) { return {}; }

    std::vector<uint32_t> cids_out;
    cids_out.reserve(cids.size());

    std::vector<color_info> res;
    std::vector<uint32_t> key;
    for (size_t i = 0; i < cids.size(); i++) {
        bool not_found = true;
        if (i < cids.size() - k) {
            key.assign(cids.begin() + i, cids.begin() + i + k);
            auto kv = frequent_map.find(key);
            if (kv != frequent_map.end()) {
                not_found = false;
                res.push_back({kv->second.begin(), kv->second.end(), kv->second.begin()});
                i += k;
            }
        }
        if (not_found) { cids_out.push_back(cids[i]); }
    }
    std::swap(cids, cids_out);
    return res;
}

struct pref_suf_bounds {
    uint32_t prefix_len{0};
    uint32_t suffix_len{0};
};

pref_suf_bounds find_max_prefix_suffix(std::vector<uint32_t>& prev_cid,
                                       std::vector<uint32_t>& new_cid) {
    if (prev_cid.empty() or new_cid.empty()) { return {0, 0}; }
    uint32_t pctr{0};
    uint32_t sctr{0};
    {
        auto pit = prev_cid.begin();
        auto cit = new_cid.begin();
        // count the length of the LCP
        for (; (*pit == *cit) and (pit < prev_cid.end()) and (cit < new_cid.end());
             ++pit, ++cit, ++pctr) {}
    }

    {
        auto pit = prev_cid.rbegin();
        auto cit = new_cid.rbegin();
        // count the length of the LCS
        for (; (*pit == *cit) and (pit < prev_cid.rend()) and (cit < new_cid.rend());
             ++pit, ++cit, ++sctr) {}
    }

    return {pctr, sctr};
}

void process_lines(index_type& index,
                   moodycamel::BlockingConcurrentQueue<std::vector<std::vector<uint32_t>>*>& q,
                   std::atomic<bool>& done, std::atomic<uint32_t>& in_flight, std::mutex& ofile_mut,
                   std::ofstream& out_file) {
    std::stringstream ss;
    int32_t buff_size{0};
    constexpr int32_t buff_thresh{100};

    frequent_map_t frequent_map;
    std::vector<uint32_t> pkey;
    std::vector<uint32_t> skey;
    std::vector<uint32_t> colors;
    std::vector<uint32_t> colors_tmp;

    std::vector<std::vector<uint32_t>>* sbatch;

    size_t total_cid_len = 0;
    size_t total_ps_len = 0;

    while (!done or in_flight > 0) {
        if (q.try_dequeue(sbatch)) {
            in_flight -= 1;
            frequent_map.clear();
            /*
               if (batch_ctr >= 5) {
               frequent_map.clear();
               batch_ctr = 0;
               }
               */
            // analyze_batch(sbatch, index, frequent_map);

            std::vector<uint32_t> prev_cids{};
            for (auto& cids : *sbatch) {
                pkey.clear();
                skey.clear();

                auto ps_info = find_max_prefix_suffix(prev_cids, cids);
                if (ps_info.prefix_len > 0) {
                    pkey.assign(cids.begin(), cids.begin() + ps_info.prefix_len);
                }
                if (ps_info.suffix_len > 0) {
                    skey.assign(cids.rbegin(), cids.rbegin() + ps_info.suffix_len);
                }

                std::vector<color_info> uncompressed_res;
                auto pkey_bak = pkey;
                auto skey_bak = skey;

                auto match_key = [&](auto& key, auto& key_bak) -> auto {
                    if (key.empty()) { return frequent_map.end(); }
                    auto pref_it = frequent_map.find(key);
                    while (!key.empty()) {
                        if (pref_it == frequent_map.end()) {
                            key.pop_back();
                            pref_it = frequent_map.find(key);
                        } else {
                            uncompressed_res.push_back({pref_it->second.begin(),
                                                        pref_it->second.end(),
                                                        pref_it->second.begin()});
                            break;
                        }
                    }

                    // if we didn't match the whole prefix, then
                    // compute the result for this prefix and
                    // put it in our hash
                    if ((key.size() >= 4) or (key.size() == key_bak.size())) {
                        std::swap(key, key_bak);
                    } else {
                        colors_tmp.clear();
                        std::vector<uint32_t> rem(key_bak.begin() + key.size(), key_bak.end());
                        intersect_color_ids(index, rem, colors_tmp);
                        uncompressed_res.push_back(
                            {colors_tmp.begin(), colors_tmp.end(), colors_tmp.begin()});
                        colors.clear();
                        intersect_uncompressed(uncompressed_res, index.num_colors(), colors);
                        frequent_map[key_bak] = colors;
                        colors.clear();
                    }
                    uncompressed_res.clear();
                    return frequent_map.find(key_bak);
                };

                auto pit = match_key(pkey, pkey_bak);
                auto sit = match_key(skey, skey_bak);

                if (pit != frequent_map.end()) {
                    uncompressed_res.push_back(
                        {pit->second.begin(), pit->second.end(), pit->second.begin()});
                }
                if (sit != frequent_map.end()) {
                    uncompressed_res.push_back(
                        {sit->second.begin(), sit->second.end(), sit->second.begin()});
                }

                total_cid_len += cids.size();
                total_ps_len += pkey_bak.size() + skey_bak.size();
                prev_cids = cids;
                std::vector<uint32_t>(cids.begin() + pkey_bak.size(), cids.end() - skey_bak.size())
                    .swap(cids);
                if (!cids.empty()) {
                    colors_tmp.clear();
                    intersect_color_ids(index, cids, colors_tmp);
                    uncompressed_res.push_back(
                        {colors_tmp.begin(), colors_tmp.end(), colors_tmp.begin()});
                }
                colors.clear();
                intersect_uncompressed(uncompressed_res, index.num_colors(), colors);
                /*
                // extract colors that are in the freuqent_map from
                // the colors vector
                std::vector<color_info> frequent_res = extract_frequent_colors(cids, frequent_map);

                // intersect what remains
                colors_tmp.clear();
                if (!cids.empty()) {
                index.intersect_color_ids(cids, colors_tmp);
                }
                // add the color vector we just computed
                if (!colors_tmp.empty()) {
                frequent_res.push_back( {colors_tmp.begin(), colors_tmp.end(), colors_tmp.begin()}
                );
                }
                // intersect_uncompressed
                intersect_uncompressed(frequent_res, index.num_colors(), colors);
                */

                if (!colors.empty()) {
                    // num_mapped_reads += 1;
                    ss << "read_name" << "\t" << colors.size();
                    for (auto c : colors) { ss << "\t" << c; }
                    ss << "\n";
                } else {
                    ss << "read_name" << "\t0\n";
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

                colors.clear();
            }
            delete sbatch;
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

    std::cerr << "total cid len: " << total_cid_len << ", total_ps_len: " << total_ps_len
              << " ratio " << static_cast<double>(total_ps_len) / total_cid_len << "\n";
}

void sketch_lines(
    moodycamel::BlockingConcurrentQueue<std::vector<std::vector<uint32_t>>*>& q,
    std::atomic<bool>& done, std::atomic<uint32_t>& in_flight, 
    std::vector<sketch::hll_t>& sketches, std::vector<uint8_t>& counts)
{
    typename sketch::hll_t::HashType hasher;
    std::vector<std::vector<uint32_t>>* sbatch = nullptr;

    while (q.try_dequeue(sbatch) or !done or (in_flight > 0)) {
        if (sbatch != nullptr) {
            in_flight -= 1;

            for (auto& cids : *sbatch) {
                // if the only thing in the cids vector is
                // 1 element (the read_id), then the color
                // output should be exactly the same as the
                // last call to index.instersect_color_ids
                // and so we don't recompute it here.
                if (cids.size() > 1) {
                    // interesect the color ids to get the colors
                    for (auto it = cids.begin() + 1; it != cids.end(); ++it){
                        sketches[*it].add(hasher.hash(cids[0]));
                        counts[*it] += (counts[*it] < std::numeric_limits<uint8_t>::max());
                    }
                } 
            }
            delete sbatch;
            sbatch = nullptr;
        }
    }
}

void process_lines_meta(
    index_type& index, moodycamel::BlockingConcurrentQueue<std::vector<std::vector<uint32_t>>*>& q,
    std::atomic<bool>& done, std::atomic<uint32_t>& in_flight, std::mutex& ofile_mut,
    std::ofstream& out_file, size_t thread_index, std::atomic<uint64_t>& num_mapped_reads,
    meta_sets& ms, std::mutex& ms_mut) {
    std::stringstream ss;
    int32_t buff_size{0};
    constexpr int32_t buff_thresh{500};

    uint64_t cached = 0;
    uint64_t saving = 0;
    uint64_t total = 0;

    std::vector<uint32_t> colors;
    std::vector<std::vector<uint32_t>>* sbatch = nullptr;

    size_t skipped_aln = 0;
    size_t total_aln = 0;
    size_t batch_ctr = 0;
    while (q.try_dequeue(sbatch) or !done or (in_flight > 0)) {
        if (sbatch != nullptr) {
            in_flight -= 1;
            batch_ctr++;

            // initially for each batch, the output
            // should be empty.
            colors.clear();

            for (auto& cids : *sbatch) {
                ++total_aln;
                // if the only thing in the cids vector is
                // 1 element (the read_id), then the color
                // output should be exactly the same as the
                // last call to index.instersect_color_ids
                // and so we don't recompute it here.
                if (cids.size() > 1) {
                    ss << "* ";
                    // interesect the color ids to get the colors
                    colors.clear();
                    std::vector<uint32_t> partial_set;
                    std::sort(cids.begin()+1, cids.end(), [&ms](uint32_t a, uint32_t b){
                            return ms.permutation[a] < ms.permutation[b];    
                        });

                    uint32_t group_id = ms.s2m_index.rank1(ms.s2m, ms.permutation[cids[1]]);
                    std::vector<uint32_t> iterators = {0}; // for the offset in intersect_color_ids
                    std::vector<color_info> cached_iterators;
                    partial_set.push_back(cids[1]);
                    for (auto it = cids.begin() + 1; it != cids.end(); ++it) {
                        uint32_t pid = ms.permutation[*it];
                        uint32_t curr_group_id = ms.s2m_index.rank1(ms.s2m, pid);
                        if (curr_group_id != group_id){
                            if (partial_set.size() > 1){
                                std::vector<uint32_t> partial_intersection;
                                ms_mut.lock();
                                if (ms.cache.count(partial_set) != 0){
                                    cached++;
                                    saving += partial_set.size() - 1;
                                } else {
                                    //ms_mut.unlock();
                                    //ms.cache[partial_set] = {};
                                    intersect_color_ids(index, partial_set, ms.cache[partial_set]);
                                    //ms_mut.lock();
                                }
                                cached_iterators.emplace_back(ms.cache[partial_set].begin(), ms.cache[partial_set].end(), ms.cache[partial_set].begin());
                                ms_mut.unlock();
                            } else {
                                assert(partial_set.size() == 1);
                                iterators.push_back(partial_set.front());
                            }
                            group_id = curr_group_id;
                            partial_set.clear();
                        }
                        partial_set.push_back(*it);
                    }

                    if (partial_set.size() > 1){
                        ms_mut.lock();
                        std::vector<uint32_t> partial_intersection;
                        if (ms.cache.count(partial_set) != 0){
                            cached_iterators.emplace_back(ms.cache[partial_set].begin(), ms.cache[partial_set].end(), ms.cache[partial_set].begin());
                            cached++;
                            saving += partial_set.size() - 1;
                        } else {
                            ms_mut.unlock();
                            intersect_color_ids(index, partial_set, partial_intersection);
                            ms_mut.lock();
                            ms.cache.emplace(partial_set, partial_intersection);
                            cached_iterators.emplace_back(ms.cache[partial_set].begin(), ms.cache[partial_set].end(), ms.cache[partial_set].begin());
                        }
                        ms_mut.unlock();
                    } else {
                        iterators.push_back(partial_set.front());
                    }

                    std::vector<uint32_t> tmp;
                    if (iterators.size() > 1){
                        intersect_color_ids(index, iterators, tmp);
                        // next_geq_intersect(iterators.begin(), iterators.end(), tmp, index.num_colors());
                        cached_iterators.emplace_back(tmp.begin(), tmp.end(), tmp.begin());
                    }
                    intersect_uncompressed(cached_iterators, index.num_colors(), colors);

                    // intersect_color_ids(index, cids, colors);
                } else {
                    ++skipped_aln;
                    ss << "# ";
                }

                if (!colors.empty()) {
                    ss << cids.front()+1 << "\t" << colors.size();
                    for (auto c : colors) { ss << "\t" << c; }
                    ss << "\n";
                } else {
                    num_mapped_reads -= 1;
                    ss << cids.front()+1 << "\t0\n";
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

void process_lines_simple(
    index_type& index, moodycamel::BlockingConcurrentQueue<std::vector<std::vector<uint32_t>>*>& q,
    std::atomic<bool>& done, std::atomic<uint32_t>& in_flight, std::mutex& ofile_mut,
    std::ofstream& out_file, size_t thread_index, std::atomic<uint64_t>& num_mapped_reads) {
    std::stringstream ss;
    int32_t buff_size{0};
    constexpr int32_t buff_thresh{500};

    std::vector<uint32_t> colors;
    std::vector<std::vector<uint32_t>>* sbatch = nullptr;

    size_t skipped_aln = 0;
    size_t total_aln = 0;
    size_t batch_ctr = 0;
    while (q.try_dequeue(sbatch) or !done or (in_flight > 0)) {
        if (sbatch != nullptr) {
            in_flight -= 1;
            batch_ctr++;

            // initially for each batch, the output
            // should be empty.
            colors.clear();

            for (auto& cids : *sbatch) {
                ++total_aln;
                // if the only thing in the cids vector is
                // 1 element (the read_id), then the color
                // output should be exactly the same as the
                // last call to index.instersect_color_ids
                // and so we don't recompute it here.
                if (cids.size() > 1) {
                    // interesect the color ids to get the colors
                    colors.clear();
                    intersect_color_ids(index, cids, colors);
                } else {
                    ++skipped_aln;
                }

                if (!colors.empty()) {
                    ss << cids.front() << "\t" << colors.size();
                    for (auto c : colors) { ss << "\t" << c; }
                    ss << "\n";
                } else {
                    num_mapped_reads -= 1;
                    ss << cids.front() << "\t0\n";
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

int do_color_map(index_type const& index, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser,
                 std::atomic<uint64_t>& num_reads, std::atomic<uint64_t>& num_mapped_reads,
                 std::ofstream& out_file, std::mutex& iomut, std::mutex& ofile_mut) {
    std::vector<uint32_t> colors;  // result of pseudo-alignment
    std::stringstream ss;
    int32_t buff_size = 0;
    constexpr int32_t buff_thresh = 50;

    auto rg = rparser.getReadGroup();
    while (rparser.refill(rg)) {
        // get info about the first read in this chunk
        // and pack it.
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
            index.pseudoalign_query_color_sets_ids(record.seq, colors);

            buff_size += 1;
            num_reads += 1;
            uint32_t record_info = file_mask | read_id;
            if (!colors.empty()) num_mapped_reads += 1;

            ss.write(reinterpret_cast<char*>(&record_info), sizeof(record_info));
            uint32_t num_colors = static_cast<uint32_t>(colors.size());
            ss.write(reinterpret_cast<char*>(&num_colors), sizeof(num_colors));

            for (auto c : colors) { ss.write(reinterpret_cast<char*>(&c), sizeof(c)); }

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
            ++read_id;
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

void sort_file(const std::string& tmp_outname, std::ofstream& output) {
    std::vector<std::vector<uint32_t>> v;
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
            v.push_back(vals);
            vals.clear();
        } else {
            // just write out the unmapped reads here
            output << read_num << "\t0\n";
        }
    }

    if (v.empty()) { return; }

    std::cerr << "begin sort\n";

    std::sort(v.begin(), v.end(),
              [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) -> bool {
                  return std::lexicographical_compare(a.begin() + 1, a.end(), b.begin() + 1,
                                                      b.end());
              });

    std::cerr << "done sort\n";

    auto curr = v.begin();
    auto next = curr++;
    size_t identical_lists = 0;
    while (next < v.end()) {
        if ((curr->size() == next->size()) and
            std::equal(curr->begin() + 1, curr->end(), next->begin() + 1)) {
            next->resize(1);  // retain only the id.
            ++identical_lists;
            ++next;
        } else {
            curr = next;
            ++next;
        }
    }

    std::cerr << "number of identical lists = " << identical_lists << "\n";
    std::cerr << "sorted!\n";
    ifile.close();
    std::ofstream ofile(tmp_outname, std::ios::trunc | std::ios::binary);
    for (auto vec_it = v.begin(); vec_it != v.end(); ++vec_it) {
        uint32_t s = vec_it->size();
        // if (s > 100) { std::cerr << "on output, size of v is " << s << "\n"; }
        ofile.write(reinterpret_cast<char*>(&s), sizeof(s));
        ofile.write(reinterpret_cast<char*>(vec_it->data()), sizeof(s) * s);
    }
    ofile.close();
}

void do_sketching(index_type& index, size_t num_threads_in, const std::string& query_filename ){
    const uint64_t p = 7;
    const uint64_t num_color_sets = index.num_color_sets();
    std::vector<sketch::hll_t> sketches (num_color_sets, sketch::hll_t(p));
    std::vector<uint8_t> counts(num_color_sets);

    std::ifstream ifile(query_filename, std::ios::binary);
    std::atomic<bool> done{false};
    std::atomic<uint32_t> in_flight{0};

    moodycamel::BlockingConcurrentQueue<std::vector<std::vector<uint32_t>>*> q(3 * num_threads_in);

    std::thread producer([&q, &ifile, &done, &in_flight]() -> void {
        std::vector<std::vector<uint32_t>>* batch = new std::vector<std::vector<uint32_t>>();
        uint32_t list_len = 0;
        size_t nbatches = 0;
        while (ifile.read(reinterpret_cast<char*>(&list_len), sizeof(list_len))) {
            if (list_len > 0) {
                // if (list_len > 100) { std::cerr << "list_len " << list_len << "\n"; }
                std::vector<uint32_t> working_vec;
                working_vec.resize(list_len);
                ifile.read(reinterpret_cast<char*>(working_vec.data()),
                           list_len * sizeof(list_len));

                // if the current batch is big enough that we want to push it
                // make sure the current element isn't a duplicate (i.e. size 1)
                // and push the batch *before* we add the next element, which
                // may be the start of a new duplicate run.
                if ((batch->size() >= 10000) and (working_vec.size() > 1)) {
                    while (!q.try_enqueue(batch)) {}
                    batch = new std::vector<std::vector<uint32_t>>();
                    in_flight += 1;
                    ++nbatches;
                }
                // always push back the current element. If it's in the middle of
                // a duplicate run, we keep growing the previous batch, otherwise
                // this is the first element of the enw batch.
                batch->push_back(working_vec);
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
    for (size_t i = 0; i < num_threads_in; ++i) {
        workers.push_back(std::thread(
            [&q, &done, &in_flight, &sketches, &counts]() -> void {
                sketch_lines(q, done, in_flight, sketches, counts);
            }));
    }

    producer.join();
    for (auto& t : workers) { t.join(); }

    {
        const uint8_t sketch_threshold = 2;
        // std::vector<sketch::hll_t> sketches;
        // sketches.reserve(num_color_sets);
        for (uint64_t set_id = 0; set_id < num_color_sets; set_id++){
            if (counts[set_id] < sketch_threshold) {
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
}

void do_intersection(index_type& index, size_t num_threads_in, const std::string& query_filename,
                     std::ofstream& output, std::atomic<uint64_t>& num_mapped_reads, 
                     meta_sets& ms) {
    std::ifstream ifile(query_filename, std::ios::binary);

    std::atomic<bool> done{false};
    std::atomic<uint32_t> in_flight{0};

    moodycamel::BlockingConcurrentQueue<std::vector<std::vector<uint32_t>>*> q(3 * num_threads_in);

    std::thread producer([&q, &ifile, &done, &in_flight]() -> void {
        std::vector<std::vector<uint32_t>>* batch = new std::vector<std::vector<uint32_t>>();
        uint32_t list_len = 0;
        size_t nbatches = 0;
        while (ifile.read(reinterpret_cast<char*>(&list_len), sizeof(list_len))) {
            if (list_len > 0) {
                // if (list_len > 100) { std::cerr << "list_len " << list_len << "\n"; }
                std::vector<uint32_t> working_vec;
                working_vec.resize(list_len);
                ifile.read(reinterpret_cast<char*>(working_vec.data()),
                           list_len * sizeof(list_len));

                // if the current batch is big enough that we want to push it
                // make sure the current element isn't a duplicate (i.e. size 1)
                // and push the batch *before* we add the next element, which
                // may be the start of a new duplicate run.
                if ((batch->size() >= 10000) and (working_vec.size() > 1)) {
                    while (!q.try_enqueue(batch)) {}
                    batch = new std::vector<std::vector<uint32_t>>();
                    in_flight += 1;
                    ++nbatches;
                }
                // always push back the current element. If it's in the middle of
                // a duplicate run, we keep growing the previous batch, otherwise
                // this is the first element of the enw batch.
                batch->push_back(working_vec);
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

    std::mutex ofmut, ms_mut;
    std::vector<std::thread> workers;
    workers.reserve(num_threads_in);
    for (size_t i = 0; i < num_threads_in; ++i) {
        workers.push_back(std::thread(
            [&index, &q, &done, &in_flight, &ofmut, &output, &num_mapped_reads, &ms, &ms_mut, i]() -> void {
                // process_lines_simple(index, q, done, in_flight, ofmut, output, i, num_mapped_reads);
                process_lines_meta(index, q, done, in_flight, ofmut, output, i, num_mapped_reads, ms, ms_mut);
            }));
    }

    producer.join();
    for (auto& t : workers) { t.join(); }
}

int intersect_colors(int argc, char** argv) {
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
    sort_file(tmp_outname, output);

    std::cout << "sketching" << std::endl;
    const uint64_t num_color_sets = index.num_color_sets();
    do_sketching(index, num_threads, tmp_outname);
    std::cout << "done sketching!" << std::endl;

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

    meta_sets ms;
    ms.s2m_index = rank1_s2m;
    ms.s2m = set2meta;
    ms.permutation.swap(permutation);

    do_intersection(index, num_threads, tmp_outname, output, num_mapped_reads, ms);
    std::remove(tmp_outname.c_str());

    std::cout << "num_mapped_reads " << num_mapped_reads << "/" << num_reads << " ("
              << (num_mapped_reads * 100.0) / num_reads << "%)" << std::endl;
    return 0;
}

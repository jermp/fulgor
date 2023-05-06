#pragma once

#include <numeric>
#include <vector>
#include <fstream>

#ifdef __linux__
#include <parallel/algorithm>
#endif

#include "util.hpp"
#include "../external/sshash/external/pthash/external/essentials/include/essentials.hpp"
#include "../external/sshash/external/pthash/external/mm_file/include/mm_file/mm_file.hpp"
#include "inverted_index.hpp"

namespace fulgor {

struct list_sorter {
    list_sorter(build_configuration const& build_config, uint32_t num_docs)
        : m_build_config(build_config)
        , m_run_identifier(std::chrono::high_resolution_clock::now().time_since_epoch().count())
        , m_num_docs(num_docs)
        , m_cur_offset(0)
        , m_num_bytes(0)
        , m_num_sorted_runs(0)
        , m_unique_color_classes_output_filename(build_config.file_base_name +
                                                 ".unique_color_classes")
        , m_map_output_filename(build_config.file_base_name + ".map") {
        m_num_bytes = build_config.ram_limit_in_GiB * essentials::GiB;
        m_num_bytes /= 2;
#ifdef __linux__
        m_num_bytes /= 2; /* since __gnu::parallel_sort is not in-place */
#endif
        m_offsets.reserve(m_num_bytes / 2 / (2 * sizeof(uint64_t)));
        m_lists.reserve(m_num_bytes / 2 / sizeof(uint32_t));
    }

    template <typename Iterator>
    void add(seg_id_t seg_id, Iterator begin, Iterator end, uint32_t size) {
        m_offsets.push_back({seg_id, m_cur_offset});
        m_lists.push_back(size);
        std::copy(begin, end, std::back_inserter(m_lists));
        m_cur_offset += size + 1;
        if (m_lists.size() * sizeof(uint32_t) + m_offsets.size() * 2 * sizeof(uint64_t) >=
            m_num_bytes) {
            sort_and_write_to_disk();
        }
    }

    void finalize() {
        sort_and_write_to_disk();  // last run
        merge_runs();
    }

    void remove_tmp_files() {
        for (uint64_t i = 0; i != m_num_sorted_runs; ++i) {
            std::string index_filename =
                util::get_tmp_filename(m_build_config.tmp_dirname, m_run_identifier, i);
            std::remove(index_filename.c_str());
        }
    }

private:
    build_configuration const& m_build_config;
    uint64_t m_run_identifier;
    uint32_t m_num_docs;
    uint64_t m_cur_offset;
    uint64_t m_num_bytes;
    uint64_t m_num_sorted_runs;
    std::vector<std::pair<seg_id_t, uint64_t>> m_offsets;  // (seg_id, offset into m_lists)
    std::vector<uint32_t> m_lists;

    /* dictionary of distinct color classes (inverted lists) */
    std::string m_unique_color_classes_output_filename;

    /* map from terms of the index (unitig ids) to dictionary entries */
    std::string m_map_output_filename;

    void sort_and_write_to_disk() {
        if (m_offsets.empty()) return;

        if (m_build_config.verbose) std::cout << " == sorting lists..." << std::endl;

#ifdef __linux__
        __gnu_parallel::sort
#else
        std::sort
#endif
            (m_offsets.begin(), m_offsets.end(), [&](auto offset_x, auto offset_y) {
                uint32_t listx_size = m_lists[offset_x.second];
                uint32_t listy_size = m_lists[offset_y.second];
                if (listx_size != listy_size) {
                    return listx_size < listy_size;  // non-decreasing length
                }
                /* lex order */
                uint32_t const* listx = m_lists.data() + offset_x.second + 1;
                uint32_t const* listy = m_lists.data() + offset_y.second + 1;
                for (uint64_t i = 0; i != listx_size; ++i, ++listx, ++listy) {
                    uint32_t x = *listx;
                    uint32_t y = *listy;
                    if (x == y) continue;
                    if (x < y) return true;
                    break;
                }
                return false;
            });

        std::string output_filename =
            util::get_tmp_filename(m_build_config.tmp_dirname, m_run_identifier, m_num_sorted_runs);
        inverted_index::builder ii_builder(m_num_docs, output_filename);

        if (m_build_config.verbose) {
            std::cout << " == writing to file '" << output_filename << "'..." << std::endl;
        }

        uint64_t num_ints = 0;
        for (auto offset : m_offsets) {
            uint32_t const* list = m_lists.data() + offset.second;
            uint32_t list_size = *list;
            list += 1;
            num_ints += list_size;
            seg_id_t seg_id = offset.first;
            ii_builder.write_list(seg_id, list, list_size);
        }

        m_lists.clear();
        m_offsets.clear();
        m_cur_offset = 0;
        m_num_sorted_runs += 1;

        ii_builder.finalize(num_ints);
    }

    void merge_runs() {
        if (m_num_sorted_runs == 0) return; /* all empty files */

        assert(m_lists.empty() and m_offsets.empty());

        std::ofstream out_map(m_map_output_filename.c_str(), std::ofstream::binary);
        if (!out_map.is_open()) throw std::runtime_error("cannot open file");

        inverted_index::builder ii_builder(m_num_docs, m_unique_color_classes_output_filename);

        /* input iterators and heap */
        constexpr bool with_cache = true;
        std::vector<inverted_index::iterator<with_cache>> iterators;
        std::vector<uint32_t> idx_heap;
        iterators.reserve(m_num_sorted_runs);
        idx_heap.reserve(m_num_sorted_runs);

        /* heap functions */
        auto stdheap_idx_comparator = [&](uint32_t i, uint32_t j) {
            auto list_i = iterators[i].list();
            auto list_j = iterators[j].list();
            if (list_i.size() != list_j.size()) return list_i.size() >= list_j.size();
            auto cache_i_begin = list_i.cache().begin();
            auto cache_j_begin = list_j.cache().begin();
            for (uint64_t i = 0; i != list_i.size(); ++i, ++cache_i_begin, ++cache_j_begin) {
                uint32_t x = *cache_i_begin;
                uint32_t y = *cache_j_begin;
                if (x == y) continue;
                if (x < y) return false;
                break;
            }
            return true;
        };
        auto advance_heap_head = [&]() {
            auto idx = idx_heap.front();
            iterators[idx].advance_to_next_list();
            if (iterators[idx].has_next()) {
                // percolate down the head
                uint64_t pos = 0;
                uint64_t size = idx_heap.size();
                while (2 * pos + 1 < size) {
                    uint64_t i = 2 * pos + 1;
                    if (i + 1 < size and stdheap_idx_comparator(idx_heap[i], idx_heap[i + 1])) ++i;
                    if (stdheap_idx_comparator(idx_heap[i], idx_heap[pos])) break;
                    std::swap(idx_heap[pos], idx_heap[i]);
                    pos = i;
                }
            } else {
                std::pop_heap(idx_heap.begin(), idx_heap.end(), stdheap_idx_comparator);
                idx_heap.pop_back();
            }
        };

        /* create the input iterators and the heap */
        std::vector<mm::file_source<uint64_t>> mm_files(m_num_sorted_runs);
        for (uint64_t i = 0; i != m_num_sorted_runs; ++i) {
            std::string index_filename =
                util::get_tmp_filename(m_build_config.tmp_dirname, m_run_identifier, i);
            mm_files[i].open(index_filename, mm::advice::sequential);
            iterators.emplace_back(mm_files[i].data(), mm_files[i].size());
            idx_heap.push_back(i);
        }
        std::make_heap(idx_heap.begin(), idx_heap.end(), stdheap_idx_comparator);

        uint64_t num_distinct_lists = 0;
        uint64_t num_processed_lists = 0;
        uint64_t num_ints_in_color_classes = 0;

        auto write_list = [&](uint32_t const* list, uint32_t list_size) {
            constexpr seg_id_t fake_seg_id = seg_id_t(-1);
            ii_builder.write_list(fake_seg_id, list, list_size);
            num_distinct_lists += 1;
            num_ints_in_color_classes += list_size;
        };

        auto write_dictionary_entry = [&](seg_id_t seg_id) {
            assert(num_distinct_lists > 0);
            uint64_t dictionary_entry = num_distinct_lists - 1;
            out_map.write(reinterpret_cast<char const*>(&seg_id), sizeof(seg_id_t));
            out_map.write(reinterpret_cast<char const*>(&dictionary_entry), sizeof(uint64_t));
        };

        std::cout << "processing lists..." << std::endl;

        /* write first list */
        auto& first_list = iterators[idx_heap.front()].list();
        auto prev_cache = first_list.cache();
        write_list(prev_cache.data(), prev_cache.size());
        write_dictionary_entry(first_list.seg_id());
        advance_heap_head();
        num_processed_lists += 1;

        while (!idx_heap.empty()) {
            auto& list = iterators[idx_heap.front()].list();
            auto curr_cache = list.cache();
            if (curr_cache.size() != prev_cache.size() or
                !std::equal(curr_cache.begin(), curr_cache.end(), prev_cache.begin())) {
                write_list(curr_cache.data(), curr_cache.size());
            }
            write_dictionary_entry(list.seg_id());
            advance_heap_head();
            prev_cache.swap(curr_cache);
            num_processed_lists += 1;
            if (num_processed_lists % 1000000 == 0 and m_build_config.verbose) {
                std::cout << " == processed " << num_processed_lists << " lists" << std::endl;
            }
        }
        if (m_build_config.verbose) {
            std::cout << " == processed " << num_processed_lists << " lists" << std::endl;
            std::cout << " == num_distinct_lists " << num_distinct_lists << std::endl;
            std::cout << " == num_ints_in_color_classes " << num_ints_in_color_classes << std::endl;
        }

        ii_builder.finalize(num_ints_in_color_classes);
        out_map.close();

        /* close tmp files and remove them */
        for (uint64_t i = 0; i != m_num_sorted_runs; ++i) mm_files[i].close();
        remove_tmp_files();
    }
};

}  // namespace fulgor
#pragma once

#include <fstream>
#include <limits>
#include <string>
#include <iostream>
#include <chrono>
#include <filesystem>

#ifdef __linux__
#include <parallel/algorithm>
#endif

#include "util.hpp"
#include "../external/sshash/external/pthash/external/essentials/include/essentials.hpp"
#include "../external/sshash/external/pthash/external/mm_file/include/mm_file/mm_file.hpp"
#include "inverted_index.hpp"

namespace fulgor {

#pragma pack(push, 4)
struct pair_t {
    seg_id_t seg_id;
    uint32_t ref_id; /* we always assume we have less than 2^32 references */
    bool operator==(const pair_t rhs) const {
        return ((*this).seg_id == rhs.seg_id) and ((*this).ref_id == rhs.ref_id);
    }
};
#pragma pack(pop)

struct inverted_index_builder {
    inverted_index_builder(build_configuration const& build_config)
        : m_build_config(build_config)
        , m_num_pairs(0)
        , m_num_sorted_runs(0)
        , m_num_unitigs(0)
        , m_curr_ref_id(0)
        , m_prev_ref_id(-1)
        , m_num_refs(0)
        , m_num_docs_in_curr_run(0)
        , m_run_identifier(std::chrono::high_resolution_clock::now().time_since_epoch().count()) {
        uint64_t available_ram = build_config.ram_limit_in_GiB * essentials::GiB;
        if (build_config.num_refs != 0) {  // if the user specified num_refs
            m_num_refs = build_config.num_refs;
            if (available_ram <= m_num_refs * 4 * sizeof(uint32_t)) {
                throw std::runtime_error("please use more RAM");
            }
            /* reserve the memory for a uncompressed list */
            available_ram -= m_num_refs * 4 * sizeof(uint32_t);
        }
        m_num_pairs = available_ram / sizeof(pair_t);
#ifdef __linux__
        m_num_pairs /= 2; /* since __gnu::parallel_sort is not in-place */
#endif
        if (m_build_config.verbose) std::cout << "m_num_pairs = " << m_num_pairs << std::endl;
        m_pairs.reserve(m_num_pairs);
        m_buffer.reserve(256);
    }

    void build() {
        try {
            essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
            timer.start();
            create_runs();
            merge_runs();
            timer.stop();
            std::cout << "** inversion took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes / "
                      << (timer.elapsed() * 1000000) / m_num_unitigs << " musec/unitig"
                      << std::endl;
        } catch (std::runtime_error const& e) {
            std::cout << e.what() << std::endl;
            remove_tmp_files();
        }
    }

private:
    build_configuration const& m_build_config;
    uint64_t m_num_pairs;
    std::vector<pair_t> m_pairs;

    /* small buffer to contain a single list */
    std::vector<uint32_t> m_buffer;

    uint64_t m_num_sorted_runs;
    uint64_t m_num_unitigs;
    uint32_t m_curr_ref_id;
    uint32_t m_prev_ref_id;
    uint32_t m_num_refs;
    uint32_t m_num_docs_in_curr_run;
    uint64_t m_run_identifier;

    void remove_tmp_files() {
        for (uint64_t i = 0; i != m_num_sorted_runs; ++i) {
            std::string index_filename =
                util::get_tmp_filename(m_build_config.tmp_dirname, m_run_identifier, i);
            std::remove(index_filename.c_str());
        }
    }

    void create_runs() {
        mm::file_source<uint64_t> in(m_build_config.file_base_name + ".cf_inv_col",
                                     mm::advice::sequential);
        uint64_t const* begin = in.data();
        uint64_t const* end = begin + in.size();
        if (!in.is_open()) throw std::runtime_error("error in opening file");

        constexpr uint64_t max_doc_id = std::numeric_limits<decltype(m_curr_ref_id)>::max();

        /*
            Format is:
            8-byte uint reference id (1-based)
            8-byte uint for number of unitig ids, say n
            n 8-byte uint for the unitig ids
        */
        while (begin != end) {
            uint64_t b = *begin;

            // Make sure we don't see an ID that exceeds what we can store.
            if (b > max_doc_id) {
                std::string msg = "found ref id [" + std::to_string(b) + "] > max value [" +
                                  std::to_string(max_doc_id) + "]";
                throw std::runtime_error(msg);
            }

            m_curr_ref_id = b;
            assert(m_curr_ref_id > 0);
            m_curr_ref_id -= 1;  // docids assigned from 0
            if (m_curr_ref_id == m_num_refs - 1) break;
            if (m_curr_ref_id != m_prev_ref_id) {
                m_num_docs_in_curr_run += 1;
                if (m_curr_ref_id and m_curr_ref_id % 100 == 0) {
                    std::cout << "processed " << m_curr_ref_id << " references" << std::endl;
                }
            }
            m_prev_ref_id = m_curr_ref_id;
            begin += 1;
            uint64_t n = *begin;
            begin += 1;
            for (uint64_t i = 0; i != n; ++i, ++begin) {
                m_pairs.push_back({static_cast<seg_id_t>(*begin), m_curr_ref_id});
            }
            m_num_unitigs += n;
            if (m_pairs.size() >= m_num_pairs) write_run_to_disk();
        }

        in.close();

        write_run_to_disk();  // last run

        if (m_build_config.verbose) {
            std::cout << "processed " << m_num_unitigs << " unitigs" << std::endl;
            std::cout << "processed " << m_curr_ref_id + 1 << " references" << std::endl;
            std::cout << "num_sorted_runs " << m_num_sorted_runs << std::endl;
        }
    }

    void write_run_to_disk() {
        if (m_pairs.empty()) return;

        if (m_build_config.verbose) std::cout << " == sorting buffer..." << std::endl;
#ifdef __linux__
        __gnu_parallel::sort
#else
        std::sort
#endif
            (m_pairs.begin(), m_pairs.end(), [](auto const& x, auto const& y) {
                if (x.seg_id != y.seg_id) return x.seg_id < y.seg_id;
                return x.ref_id < y.ref_id;
            });

        std::string output_filename =
            util::get_tmp_filename(m_build_config.tmp_dirname, m_run_identifier, m_num_sorted_runs);

        if (m_build_config.verbose) {
            std::cout << " == writing to file '" << output_filename << "'..." << std::endl;
        }

        inverted_index::builder ii_builder(m_num_docs_in_curr_run, output_filename);

        uint64_t num_ints = 0;
        uint64_t num_lists = 0;
        for (uint64_t i = 0; i != m_pairs.size();) {
            seg_id_t seg_id = m_pairs[i].seg_id;
            uint64_t end = i;
            while (end != m_pairs.size() and m_pairs[end].seg_id == seg_id) ++end;

            uint32_t prev_ref_id = -1;
            for (uint64_t begin = i; begin != end; ++begin) {
                uint32_t ref_id = m_pairs[begin].ref_id;
                assert(m_pairs[begin].seg_id == seg_id);
                if (ref_id != prev_ref_id) m_buffer.push_back(ref_id);  // skip equal ids
                prev_ref_id = ref_id;
            }

            assert(m_buffer.size() >= 1);
            ii_builder.write_list(seg_id, m_buffer.begin(), m_buffer.size());
            num_ints += m_buffer.size();
            num_lists += 1;

            m_buffer.clear();
            i = end;
        }
        ii_builder.finalize(num_ints);

        if (m_build_config.verbose) {
            std::cout << " == DONE" << std::endl;
            std::cout << "    (written " << num_lists << " lists)" << std::endl;
            std::cout << "    (written " << num_ints << " ints)" << std::endl;
        }

        m_num_docs_in_curr_run = 1;
        m_num_sorted_runs += 1;
        m_pairs.clear();
    }

    void merge_runs() {
        if (m_num_sorted_runs == 0) return; /* all empty files */

        std::string output_filename = m_build_config.file_base_name + ".inv_idx";

        if (m_num_sorted_runs == 1) {
            try {
                std::filesystem::rename(
                    util::get_tmp_filename(m_build_config.tmp_dirname, m_run_identifier, 0).c_str(),
                    output_filename.c_str());
            } catch (std::filesystem::filesystem_error const& e) { throw e; }
            return;
        }

        assert(m_num_sorted_runs > 1);

        uint32_t num_refs = m_curr_ref_id + 1;
        inverted_index::builder ii_builder(num_refs, output_filename);

        /* input iterators and heap */
        constexpr bool with_cache = false;
        std::vector<inverted_index::iterator<with_cache>> iterators;
        std::vector<uint32_t> idx_heap;
        iterators.reserve(m_num_sorted_runs);
        idx_heap.reserve(m_num_sorted_runs);

        /* heap comparator */
        auto stdheap_idx_comparator = [&](uint32_t i, uint32_t j) {
            return iterators[i].list().seg_id() > iterators[j].list().seg_id();
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

        /* percolate down the head */
        auto advance_heap_head = [&]() {
            auto idx = idx_heap.front();
            iterators[idx].advance_to_next_list();
            if (iterators[idx].has_next()) {
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
                /* file has been read completely, so close it and remove it */
                mm_files[idx].close();
                std::string index_filename =
                    util::get_tmp_filename(m_build_config.tmp_dirname, m_run_identifier, idx);
                std::remove(index_filename.c_str());
                std::pop_heap(idx_heap.begin(), idx_heap.end(), stdheap_idx_comparator);
                idx_heap.pop_back();
            }
        };

        std::cout << "writing lists..." << std::endl;

        seg_id_t seg_id = iterators[idx_heap.front()].list().seg_id();

        m_buffer.clear();
        m_buffer.reserve(2 * num_refs);
        uint64_t num_lists = 0;
        uint64_t num_integers = 0;

        auto write = [&]() {
            std::sort(m_buffer.begin(), m_buffer.end());
            auto it = std::unique(m_buffer.begin(), m_buffer.end());
            uint32_t list_size = std::distance(m_buffer.begin(), it);
            ii_builder.write_list(seg_id, m_buffer.begin(), list_size);
            m_buffer.clear();
            ++num_lists;
            num_integers += list_size;
            if (num_lists % 1000000 == 0 and m_build_config.verbose) {
                std::cout << " == written " << num_lists << " lists (" << num_integers
                          << " integers)" << std::endl;
            }
        };

        while (!idx_heap.empty()) {
            auto& it = iterators[idx_heap.front()];
            seg_id_t new_seg_id = it.list().seg_id();
            assert(new_seg_id >= seg_id);
            if (new_seg_id > seg_id) {
                write();
                seg_id = new_seg_id;
            }
            std::copy(it.list().begin(), it.list().end(), std::back_inserter(m_buffer));
            advance_heap_head();
        }

        assert(m_buffer.size());
        write();

        ii_builder.finalize(num_integers);

        if (m_build_config.verbose) {
            std::cout << " == written " << num_lists << " lists (" << num_integers << " integers)"
                      << std::endl;
        }

        /* All files should have been already closed and removed, but double-check anyway */
        for (uint64_t i = 0; i != m_num_sorted_runs; ++i) mm_files[i].close();
        remove_tmp_files();
    }
};

}  // namespace fulgor

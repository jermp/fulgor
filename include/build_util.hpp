#pragma once

#include "external/sketch/include/sketch/hll.h"
#include "external/kmeans/include/kmeans.hpp"
#include "color_sets/meta.hpp"

namespace fulgor {

void build_reference_sketches(index_type const& index,
                              uint64_t p,                  // use 2^p bytes per HLL sketch
                              uint64_t num_threads,        // num. threads for construction
                              std::string output_filename  // where the sketches will be serialized
) {
    assert(num_threads > 0);

    const uint64_t num_colors = index.num_colors();
    typename sketch::hll_t::HashType hasher;
    auto const& u2c = index.get_u2c();
    auto const& ccs = index.get_color_sets();
    const uint64_t num_color_sets = ccs.num_color_sets();
    const uint64_t num_ones = u2c.num_ones();
    assert(num_color_sets == num_ones);

    if (num_ones < num_threads) {
        throw std::runtime_error("there are only " + std::to_string(num_color_sets) +
                                 ": reduce the number of threads.");
    }

    std::vector<std::vector<sketch::hll_t>> thread_sketches(
        num_threads, std::vector<sketch::hll_t>(num_colors, sketch::hll_t(p)));

    struct slice {
        uint64_t begin;                         // start position in u2c
        uint64_t color_id_begin, color_id_end;  // [..)
    };
    std::vector<slice> thread_slices;

    /* compute load */
    uint64_t load = 0;
    {
        uint64_t pop_count = 0;
        uint64_t prev_pos = 0;
        pthash::bit_vector::unary_iterator unary_it(u2c);
        for (uint64_t color_id = 0; color_id != num_color_sets; ++color_id) {
            uint64_t curr_pos = pop_count != num_ones ? unary_it.next() : (u2c.size() - 1);
            uint64_t num_unitigs = curr_pos - prev_pos + 1;
            auto it = ccs.color_set(color_id);
            uint64_t size = it.size();
            load += size * num_unitigs;
            pop_count += 1;
            prev_pos = curr_pos + 1;
        }
    }

    const uint64_t load_per_thread = load / num_threads;
    if (load_per_thread == 0) {
        throw std::runtime_error("load is too small: reduce the number of threads");
    }

    {
        uint64_t prev_pos = 0;
        pthash::bit_vector::unary_iterator unary_it(u2c);
        slice s;
        s.begin = 0;
        s.color_id_begin = 0;
        uint64_t cur_load = 0;

        for (uint64_t color_id = 0; color_id != num_color_sets; ++color_id) {
            uint64_t curr_pos = color_id != num_color_sets - 1 ? unary_it.next() : (u2c.size() - 1);

            uint64_t num_unitigs = curr_pos - prev_pos + 1;
            auto it = ccs.color_set(color_id);
            uint64_t size = it.size();
            cur_load += size * num_unitigs;
            prev_pos = curr_pos + 1;

            if (cur_load >= load_per_thread or color_id == num_color_sets - 1) {
                s.color_id_end = color_id + 1;
                thread_slices.push_back(s);
                s.begin = prev_pos;
                s.color_id_begin = color_id + 1;
                cur_load = 0;
            }
        }

        num_threads = thread_slices.size();
    }

    auto exe = [&](uint64_t thread_id) {
        assert(thread_id < thread_slices.size());
        auto& sketches = thread_sketches[thread_id];
        auto s = thread_slices[thread_id];
        uint64_t prev_pos = s.begin;
        std::vector<uint64_t> hashes;
        pthash::bit_vector::unary_iterator unary_it(u2c, s.begin);
        for (uint64_t color_id = s.color_id_begin; color_id != s.color_id_end; ++color_id) {
            uint64_t curr_pos = color_id != num_color_sets - 1 ? unary_it.next() : (u2c.size() - 1);
            auto it = ccs.color_set(color_id);
            const uint64_t size = it.size();
            hashes.reserve(curr_pos - prev_pos + 1);
            for (uint64_t unitig_id = prev_pos; unitig_id <= curr_pos; ++unitig_id) {
                assert(unitig_id < u2c.size());
                assert(index.u2c(unitig_id) == color_id);
                hashes.push_back(hasher.hash(unitig_id));
            }
            for (uint64_t i = 0; i != size; ++i, ++it) {
                uint32_t ref_id = *it;
                assert(ref_id < num_colors);
                for (auto hash : hashes) sketches[ref_id].add(hash);
            }
            prev_pos = curr_pos + 1;
            hashes.clear();
        }
    };

    std::vector<std::thread> threads(num_threads);
    for (uint64_t thread_id = 0; thread_id != num_threads; ++thread_id) {
        threads[thread_id] = std::thread(exe, thread_id);
    }
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    /* merge sketches into thread_sketches[0] */
    for (uint64_t i = 0; i != num_colors; ++i) {
        auto& sketch = thread_sketches[0][i];
        for (uint64_t thread_id = 1; thread_id != num_threads; ++thread_id) {
            sketch += thread_sketches[thread_id][i];
        }
    }

    std::ofstream out(output_filename, std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("cannot open file");
    const uint64_t num_bytes = 1ULL << p;
    out.write(reinterpret_cast<char const*>(&num_bytes), 8);
    out.write(reinterpret_cast<char const*>(&num_colors), 8);
    for (auto const& x : thread_sketches[0]) {
        assert(x.m() == num_bytes);
        assert(x.m() == x.core().size());
        uint8_t const* data = x.data();
        out.write(reinterpret_cast<char const*>(data), num_bytes);
    }
    out.close();
}

template <typename Iterator>
void build_colors_sketches_sliced(
    uint64_t num_colors, uint64_t num_color_sets, function<Iterator(uint64_t)> colors,
    uint64_t p,                   // use 2^p bytes per HLL sketch
    uint64_t num_threads,         // num. threads for construction
    std::string output_filename,  // where the sketches will be serialized
    double left, double right)    //
{
    assert(num_threads > 0);

    const double min_size = left * num_colors;
    const double max_size = right * num_colors;
    assert(min_size >= 0);
    assert(max_size <= num_colors);

    if (num_color_sets < num_threads) { num_threads = num_color_sets; }

    std::vector<Iterator> filtered_colors;
    std::vector<uint64_t> filtered_colors_ids;
    for (uint64_t color_id = 0; color_id != num_color_sets; ++color_id) {
        auto it = colors(color_id);
        uint64_t size = it.size();
        if (size > min_size && size <= max_size) {
            filtered_colors.push_back(it);
            filtered_colors_ids.push_back(color_id);
        }
    }
    const uint64_t partition_size = filtered_colors.size();

    struct slice {
        uint64_t begin, end;  // [..)
    };
    std::vector<slice> thread_slices;

    uint64_t load = 0;
    {
        for (auto it : filtered_colors) { load += it.size(); }
    }

    uint64_t load_per_thread = load / num_threads;
    {
        slice s;
        s.begin = 0;
        uint64_t curr_load = 0;

        for (uint64_t i = 0; i != partition_size; ++i) {
            auto it = filtered_colors[i];
            curr_load += it.size();
            if (curr_load >= load_per_thread || i == partition_size - 1) {
                s.end = i + 1;
                thread_slices.push_back(s);
                s.begin = i + 1;
                curr_load = 0;
            }
        }
        assert(thread_slices.size() <= num_threads);
    }
    num_threads = thread_slices.size();
    std::vector<std::vector<sketch::hll_t>> thread_sketches(num_threads);

    auto exe = [&](uint64_t thread_id) {
        assert(thread_id < thread_slices.size());
        auto& sketches = thread_sketches[thread_id];
        auto s = thread_slices[thread_id];
        sketches = std::vector<sketch::hll_t>(s.end - s.begin, sketch::hll_t(p));

        for (uint64_t color_id = s.begin; color_id != s.end; ++color_id) {
            auto it = filtered_colors[color_id];
            const uint64_t size = it.size();
            assert(size > 0);
            for (uint64_t i = 0; i < size; ++i, ++it) {
                uint64_t ref_id = *it;
                assert(ref_id < num_colors);
                sketches[color_id - s.begin].addh(ref_id);
            }
        }
    };

    std::vector<std::thread> threads(num_threads);
    for (uint64_t thread_id = 0; thread_id != num_threads; ++thread_id) {
        threads[thread_id] = std::thread(exe, thread_id);
    }
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    std::ofstream out(output_filename, std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("cannot open file");
    const uint64_t num_bytes = 1ULL << p;
    out.write(reinterpret_cast<char const*>(&num_bytes), 8);
    out.write(reinterpret_cast<char const*>(&num_colors), 8);
    out.write(reinterpret_cast<char const*>(&partition_size), 8);
    for (auto const color_id : filtered_colors_ids) {
        out.write(reinterpret_cast<char const*>(&color_id), 8);
    }
    for (auto const& sketch : thread_sketches) {
        for (auto const& x : sketch) {
            assert(x.m() == num_bytes);
            assert(x.m() == x.core().size());
            uint8_t const* data = x.data();
            out.write(reinterpret_cast<char const*>(data), num_bytes);
        }
    }
    out.close();
}

}  // namespace fulgor

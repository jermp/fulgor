#pragma once

#include "sketch/hll.h"
#include "kmeans.hpp"
#include "index_types.hpp"

namespace fulgor {

template <typename ColorClasses>
void build_partial_color_sketches(
    index<ColorClasses>& index,
    uint64_t p,                  // use 2^p bytes per HLL sketch
    uint64_t num_threads,        // num. threads for construction
    std::string output_filename  // where the sketches will be serialized
) {
    assert(num_threads > 0);

    typename sketch::hll_t::HashType hasher;
    const uint64_t num_docs = index.num_docs();
    const uint64_t num_color_classes = index.num_color_classes();

    cout << "num_color_classes: " << num_color_classes << '\n';

    if (num_color_classes < num_threads) {
        throw std::runtime_error("there are only " + std::to_string(num_color_classes) +
                                 ": reduce the number of threads.");
    }

    std::vector<std::vector<sketch::hll_t>> thread_sketches(
        num_threads, std::vector<sketch::hll_t>(num_color_classes, sketch::hll_t(p)));

    struct slice {
        uint64_t begin, end;  // [..)
    };
    std::vector<slice> thread_slices;

    uint64_t load = 0;
    {
        for (uint64_t color_id = 0; color_id != num_color_classes; ++color_id) {
            auto it = index.colors(color_id);
            load += it.size();
        }
    }

    uint64_t load_per_thread = load / num_threads;

    {
        slice s;
        s.begin = 0;
        uint64_t cur_load = 0;

        for (uint64_t color_id = 0; color_id != num_color_classes; ++color_id) {
            auto it = index.colors(color_id);
            cur_load += it.size();
            if (cur_load >= load_per_thread || color_id == num_color_classes - 1) {
                s.end = color_id + 1;
                thread_slices.push_back(s);
                s.begin = color_id + 1;
                cur_load = 0;
            }
        }
        assert(thread_slices.size() == num_threads);
    }

    auto exe = [&](uint64_t thread_id) {
        assert(thread_id < thread_slices.size());
        auto& sketches = thread_sketches[thread_id];
        auto s = thread_slices[thread_id];

        for (uint64_t color_id = s.begin; color_id != s.end; ++color_id) {
            auto it = index.colors(color_id);
            const uint64_t size = it.size();
            for (uint64_t i = 0; i < size; ++i, ++it) {
                uint64_t ref_id = *it;
                assert(ref_id < num_docs);
                sketches[color_id].add(ref_id);
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

    /* merge sketches into thread_sketches[0] */
    for (uint64_t i = 0; i != num_color_classes; ++i) {
        auto& sketch = thread_sketches[0][i];
        for (uint64_t thread_id = 1; thread_id != num_threads; ++thread_id) {
            sketch += thread_sketches[thread_id][i];
        }
    }

    uint64_t i = 0;
    std::ofstream out(output_filename, std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("cannot open file");
    const uint64_t num_bytes = 1ULL << p;
    out.write(reinterpret_cast<char const*>(&num_bytes), 8);
    out.write(reinterpret_cast<char const*>(&num_docs), 8);
    for (auto const& x : thread_sketches[0]) {
        assert(x.m() == num_bytes);
        assert(x.m() == x.core().size());
        uint8_t const* data = x.data();
        std::cout << "color_sketch_" << i++ << ' ' << unsigned(*data) << '\n';
        out.write(reinterpret_cast<char const*>(data), num_bytes);
    }
    std::cout << '\n';
    out.close();

    cout << '\n';
}

}  // namespace fulgor
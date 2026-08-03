#pragma once

#include "external/sketch/include/sketch/hll.h"

namespace fulgor {

inline void merge_sketches(sketch::hll_t& dest, const sketch::hll_t& src) {
    std::transform(dest.mutable_core().begin(), dest.mutable_core().end(), src.core().begin(),
                   dest.mutable_core().begin(),
                   [](uint8_t x, uint8_t y) { return std::max(x, y); });  //[cite: 1]
    dest.not_ready();                                                     //[cite: 1]
}

inline void build_reference_sketches(
    const uint64_t num_colors,
    const uint64_t p,                   // use 2^p bytes per HLL sketch
    const uint64_t num_threads,         // num. threads for construction
    std::string const& input_basename,  // where the sketches will be serialized
    std::string const& output_filename) {
    assert(num_threads > 0);
    const uint64_t max_queue_size = num_threads * 2;

    std::vector sketches(num_colors, sketch::hll_t(p));  // 10^6 colors & p == 10 -> 1GB
    std::vector<std::mutex> mutexes(num_colors);

    cdbg::unitigs_color_set_stream stream(input_basename, max_queue_size);
    stream.start();

    auto process = [&] {
        for (auto opt = stream.get(); opt != std::nullopt; opt = stream.get()) {
            sketch::hll_t sketch(p);
            auto& [unitig_start, num_unitigs, cs_id, color_set] = opt.value();
            for (uint64_t unitig_id = unitig_start; unitig_id < unitig_start + num_unitigs;
                 ++unitig_id) {
                sketch.addh(unitig_id);
            }
            for (const auto color : color_set) {
                std::lock_guard lock(mutexes[color]);
                merge_sketches(sketches[color], sketch);
            }
        }
    };

    std::vector<std::thread> threads(num_threads - 1);
    for (uint64_t thread_id = 0; thread_id != num_threads - 1; ++thread_id) {
        threads[thread_id] = std::thread(process);
    }
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    std::ofstream out(output_filename, std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("cannot open file");
    const uint64_t num_bytes = 1ULL << p;
    out.write(reinterpret_cast<char const*>(&num_bytes), 8);
    out.write(reinterpret_cast<char const*>(&num_colors), 8);
    for (auto const& x : sketches) {
        assert(x.m() == num_bytes);
        assert(x.m() == x.core().size());
        uint8_t const* data = x.data();
        out.write(reinterpret_cast<char const*>(data), num_bytes);
    }
    out.close();
}

template <typename Index>
void build_colors_sketches_sliced(
    Index const& index,
    uint64_t p,                   // use 2^p bytes per HLL sketch
    uint64_t num_threads,         // num. threads for construction
    std::string output_filename,  // where the sketches will be serialized
    double left, double right)    //
{
    assert(num_threads > 0);

    const uint64_t num_colors = index.num_colors();
    const uint64_t num_color_sets = index.num_color_sets();

    const double min_size = left * num_colors;
    const double max_size = right * num_colors;
    assert(min_size >= 0);
    assert(max_size <= num_colors);

    if (num_color_sets < num_threads) {
        num_threads = num_color_sets;
    }

    uint64_t load = 0;
    std::vector<uint64_t> filtered_colors_ids;
    filtered_colors_ids.reserve(num_color_sets);
    for (uint64_t color_id = 0; color_id != num_color_sets; ++color_id) {
        auto it = index.color_set(color_id);
        uint64_t size = it.size();
        if (size > min_size && size <= max_size) {
            load += size;
            filtered_colors_ids.push_back(color_id);
        }
    }
    const uint64_t partition_size = filtered_colors_ids.size();

    struct slice {
        uint64_t begin, end;  // [..)
    };
    std::vector<slice> thread_slices;

    uint64_t load_per_thread = load / num_threads;
    {
        slice s;
        s.begin = 0;
        uint64_t curr_load = 0;

        for (uint64_t i = 0; i != partition_size; ++i) {
            auto color_id = filtered_colors_ids[i];
            auto it = index.color_set(color_id);
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

        for (uint64_t i = s.begin; i != s.end; ++i) {
            auto color_id = filtered_colors_ids[i];
            auto it = index.color_set(color_id);
            const uint64_t size = it.size();
            assert(size > 0);
            for (uint64_t j = 0; j < size; ++j, ++it) {
                uint64_t ref_id = *it;
                assert(ref_id < num_colors);
                sketches[i - s.begin].addh(ref_id);
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

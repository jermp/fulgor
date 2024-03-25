#pragma once

#include "sketch/hll.h"
#include "external/kmeans/include/kmeans.hpp"
#include "index_types.hpp"

namespace fulgor {

template <typename ColorClass>
void sketch_color_lists(index<ColorClass>& index,
                        uint64_t p,                   // use 2^p bytes per HLL sketch
                        uint64_t num_threads,         // num. threads for construction
                        double left, double right) {
    assert(num_threads > 0);

    const uint64_t num_docs = index.num_docs();
    const uint64_t num_color_classes = index.num_color_classes();
    const double min_colors = left * num_docs;
    const double max_colors = right * num_docs;

    if (num_color_classes < num_threads) {
        throw std::runtime_error("there are only " + std::to_string(num_color_classes) +
                                 ": reduce the number of threads.");
    }

    std::vector<typename ColorClass::forward_iterator> filtered_colors;
    std::vector<uint64_t> filtered_colors_ids;
    for (uint64_t color_id = 0; color_id != num_color_classes; ++color_id) {
        auto it = index.colors(color_id);
        uint64_t size = it.size();
        if (size > min_colors && size <= max_colors) {
            filtered_colors.push_back(it);
            filtered_colors_ids.push_back(color_id);
        }
    }
    const uint64_t num_filtered = filtered_colors.size();

    std::vector<std::vector<sketch::hll_t>> thread_sketches(
        num_threads, std::vector<sketch::hll_t>(num_filtered, sketch::hll_t(p)));

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
        uint64_t cur_load = 0;

        for (uint64_t i = 0; i != num_filtered; ++i) {
            auto it = filtered_colors[i];
            cur_load += it.size();
            if (cur_load >= load_per_thread || i == num_filtered - 1) {
                s.end = i + 1;
                thread_slices.push_back(s);
                s.begin = i + 1;
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
            auto it = filtered_colors[color_id];
            const uint64_t size = it.size();
            for (uint64_t i = 0; i < size; ++i, ++it) {
                uint64_t ref_id = *it;
                assert(ref_id < num_docs);
                sketches[color_id].addh(ref_id);
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
    for (uint64_t i = 0; i != num_filtered; ++i) {
        auto& sketch = thread_sketches[0][i];
        for (uint64_t thread_id = 1; thread_id != num_threads; ++thread_id) {
            sketch += thread_sketches[thread_id][i];
        }
    }

    std::ofstream out("./sketch.tmp", std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("cannot open file");
    const uint64_t num_bytes = 1ULL << p;
    out.write(reinterpret_cast<char const*>(&num_bytes), 8);
    out.write(reinterpret_cast<char const*>(&num_docs), 8);
    out.write(reinterpret_cast<char const*>(&num_filtered), 8);
    for(auto const color_id: filtered_colors_ids){
        out.write(reinterpret_cast<char const*>(&color_id), 8);
    }
    for (auto const& x : thread_sketches[0]) {
        assert(x.m() == num_bytes);
        assert(x.m() == x.core().size());
        uint8_t const* data = x.data();
        out.write(reinterpret_cast<char const*>(data), num_bytes);
    }
    out.close();
}

template <typename ColorClass>
void cluster_sketches(index<ColorClass>& index, std::string output_filename) {

    std::ifstream in("./sketch.tmp", std::ios::binary);
    if (!in.is_open()) throw std::runtime_error("error in opening file");

    std::vector<kmeans::point> points;
    std::vector<uint64_t> color_ids;
    uint64_t num_bytes_per_point;
    uint64_t num_points;
    uint64_t num_docs;
    in.read(reinterpret_cast<char*>(&num_bytes_per_point), sizeof(uint64_t));
    in.read(reinterpret_cast<char*>(&num_docs), sizeof(uint64_t));
    in.read(reinterpret_cast<char*>(&num_points), sizeof(uint64_t));
    points.resize(num_points, kmeans::point(num_bytes_per_point));
    color_ids.resize(num_points);
    for (uint64_t i = 0; i != num_points; ++i){
        in.read(reinterpret_cast<char*>(&color_ids[i]), sizeof(uint64_t));
    }
    for (auto& point : points) {
        in.read(reinterpret_cast<char*>(point.data()), num_bytes_per_point);
    }
    in.close();
    std::remove("./sketch.tmp");

    kmeans::clustering_parameters params;
    constexpr float min_delta = 0.0001;
    constexpr float max_iteration = 10;
    constexpr uint64_t min_cluster_size = 50;
    constexpr uint64_t seed = 42;
    params.set_min_delta(min_delta);
    params.set_max_iteration(max_iteration);
    params.set_min_cluster_size(min_cluster_size);
    params.set_random_seed(seed);
    auto clustering_data = kmeans::kmeans_divisive(points.begin(), points.end(), params);

    std::cout << "** clustering completed" << std::endl;

    const uint64_t m_num_partitions = clustering_data.num_clusters;
    std::vector<uint64_t> m_partition_size(m_num_partitions + 1, 0);
    for (auto c : clustering_data.clusters) m_partition_size[c] += 1;

    /* prefix sum */
    {
        uint64_t val = 0;
        for (auto& size : m_partition_size) {
            uint64_t tmp = size;
            size = val;
            val += tmp;
        }
    }

    const uint64_t num_color_classes = num_points;

    auto clusters_pos = m_partition_size;
    std::vector<uint64_t> m_permutation(num_color_classes);
    assert(clustering_data.clusters.size() == num_color_classes);
    for (uint64_t i = 0; i != num_color_classes; ++i) {
        uint64_t cluster_id = clustering_data.clusters[i];
        m_permutation[i] = clusters_pos[cluster_id];
        clusters_pos[cluster_id] += 1;
    }

    std::vector<uint64_t> m_color_classes_ids(num_color_classes);
    for (uint64_t i = 0; i != num_color_classes; ++i) { m_color_classes_ids[m_permutation[i]] = color_ids[i]; }

    std::cout << "Computed " << m_num_partitions << " partitions\n";

    std::ofstream cluster_dump(output_filename + ".clst", std::ios::binary);
    cluster_dump.write(reinterpret_cast<char const*>(&num_docs), 8);
    cluster_dump.write(reinterpret_cast<char const*>(&num_color_classes), 8);
    cluster_dump.write(reinterpret_cast<char const*>(&m_num_partitions), 8);

    for(uint64_t i = 1; i < m_partition_size.size(); ++i){
        cluster_dump.write(reinterpret_cast<char const*>(&m_partition_size[i]), 8);
    }

    for (uint64_t i = 0; i != num_color_classes; i++) {
        // cluster_dump << "cluster=" << cluster << "; color=" << m_color_classes_ids[i] << ": ";
        auto it = index.colors(m_color_classes_ids[i]);
        const uint64_t size = it.size();
        cluster_dump.write(reinterpret_cast<char const*>(&size), 8);
        for (uint64_t j = 0; j < size; ++j, ++it) {
            uint64_t val = *it;
            cluster_dump.write(reinterpret_cast<char const*>(&val), 8);
        }
    }
    cluster_dump.close();
}

}  // namespace fulgor
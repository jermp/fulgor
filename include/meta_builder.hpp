#pragma once

#include "index.hpp"
#include "build_util.hpp"

namespace fulgor {

template <typename ColorClasses>
struct index<ColorClasses>::meta_builder {
    meta_builder() {}

    struct partition_endpoint {
        uint64_t begin, end;
    };

    meta_builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        uint64_t num_partitions = 0;
        uint64_t max_partition_size = 0;
        uint64_t num_integers_in_metacolors = 0;
        uint64_t num_lists = 0;

        std::vector<uint32_t> permutation;
        std::vector<uint32_t> permuted_list;

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        essentials::logger("step 1. loading index to be partitioned...");
        essentials::load(m_index, m_build_config.index_filename_to_partition.c_str());
        essentials::logger("DONE");

        const uint64_t num_docs = m_index.num_docs();
        const uint64_t num_color_classes = m_index.num_color_classes();

        {
            essentials::logger("step 2.1. build sketches");
            timer.start();
            constexpr uint64_t p = 10;  // use 2^p bytes per HLL sketch
            build_reference_sketches(m_index, p, m_build_config.tmp_dirname + "/sketches.bin");
            timer.stop();
            std::cout << "** building sketches took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 2.2. clustering sketches");
            timer.start();

            std::ifstream in(m_build_config.tmp_dirname + "/sketches.bin", std::ios::binary);
            if (!in.is_open()) throw std::runtime_error("error in opening file");

            std::vector<kmeans::point> points;
            uint64_t num_bytes_per_point = 0;
            uint64_t num_points = 0;
            in.read(reinterpret_cast<char*>(&num_bytes_per_point), sizeof(uint64_t));
            in.read(reinterpret_cast<char*>(&num_points), sizeof(uint64_t));
            points.resize(num_points, kmeans::point(num_bytes_per_point));
            for (auto& point : points) {
                in.read(reinterpret_cast<char*>(point.data()), num_bytes_per_point);
            }
            in.close();

            std::remove((m_build_config.tmp_dirname + "/sketches.bin").c_str());

            kmeans::clustering_parameters params;

            /* kmeans_divisive */
            constexpr float min_delta = 0.0001;
            constexpr float max_iteration = 10;
            constexpr uint64_t min_cluster_size = 50;
            constexpr uint64_t seed = 13;
            params.set_min_delta(min_delta);
            params.set_max_iteration(max_iteration);
            params.set_min_cluster_size(min_cluster_size);
            params.set_random_seed(seed);
            auto clustering_data = kmeans::kmeans_divisive(points, params);

            timer.stop();
            std::cout << "** clustering sketches took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();

            num_partitions = clustering_data.num_clusters;
            std::cout << "num_partitions = " << num_partitions << std::endl;

            m_partition_size.resize(num_partitions + 1, 0);
            for (auto c : clustering_data.clusters) { m_partition_size[c] += 1; }

            /* take prefix sums */
            uint64_t val = 0;
            for (auto& size : m_partition_size) {
                if (size > max_partition_size) max_partition_size = size;
                uint64_t tmp = size;
                size = val;
                val += tmp;
            }

            m_hashes.resize(num_partitions);

            /* build permutation */
            auto counts = m_partition_size;  // copy
            permutation.resize(num_docs);
            assert(clustering_data.clusters.size() == num_docs);
            for (uint64_t i = 0; i != num_docs; ++i) {
                uint32_t cluster_id = clustering_data.clusters[i];
                permutation[i] = counts[cluster_id];
                counts[cluster_id] += 1;
            }
        }

        {
            essentials::logger("step 2.1. build colors");
            timer.start();

            std::vector<uint32_t> partial_color;
            permuted_list.reserve(num_docs);
            partial_color.reserve(max_partition_size);

            typename ColorClasses::builder colors_builder;

            colors_builder.init_colors_builder(num_docs, num_partitions);
            for (uint64_t partition_id = 0; partition_id != num_partitions; ++partition_id) {
                auto endpoints = partition_endpoints(partition_id);
                uint64_t num_docs_in_partition = endpoints.end - endpoints.begin;
                colors_builder.init_color_partition(partition_id, num_docs_in_partition);
            }

            uint64_t partition_id = 0;

            auto hash_and_compress = [&]() {
                auto hash = util::hash128(reinterpret_cast<char const*>(partial_color.data()),
                                          partial_color.size() * sizeof(uint32_t));
                if (auto it = m_hashes[partition_id].find(hash);
                    it == m_hashes[partition_id].cend()) {
                    uint32_t id = m_hashes[partition_id].size();
                    m_hashes[partition_id].insert({hash, id});
                    colors_builder.process_colors(partition_id, partial_color.data(),
                                                  partial_color.size());
                }
                partial_color.clear();
                num_integers_in_metacolors += 1;
            };

            for (uint64_t color_class_id = 0; color_class_id != num_color_classes;
                 ++color_class_id) {
                /* permute list */
                permuted_list.clear();
                auto it = m_index.colors(color_class_id);
                uint64_t list_size = it.size();
                for (uint64_t i = 0; i != list_size; ++i, ++it) {
                    uint32_t ref_id = *it;
                    permuted_list.push_back(permutation[ref_id]);
                }
                std::sort(permuted_list.begin(), permuted_list.end());

                /* partition list */
                partition_id = 0;
                partition_endpoint curr_partition = partition_endpoints(0);
                assert(partial_color.empty());

                for (uint64_t i = 0; i != list_size; ++i) {
                    uint32_t ref_id = permuted_list[i];
                    while (ref_id >= curr_partition.end) {
                        if (!partial_color.empty()) hash_and_compress();
                        partition_id += 1;
                        curr_partition = partition_endpoints(partition_id);
                    }
                    assert(ref_id >= curr_partition.begin);
                    partial_color.push_back(ref_id - curr_partition.begin);
                }
                if (!partial_color.empty()) hash_and_compress();
            }

            uint64_t prev_id = 0;
            m_num_lists_in_partition.reserve(num_partitions);
            for (partition_id = 0; partition_id != num_partitions; ++partition_id) {
                uint64_t num_lists_in_partition = m_hashes[partition_id].size();
                num_lists += num_lists_in_partition;
                m_num_lists_in_partition.push_back(num_lists_in_partition);
                std::cout << "num_lists in partition-" << partition_id << ": "
                          << num_lists_in_partition << std::endl;
                for (auto& p : m_hashes[partition_id]) { p.second += prev_id; }
                prev_id += num_lists_in_partition;
            }

            std::cout << "total num_colors = " << num_lists << std::endl;

            essentials::logger("step 2.2. build metacolors");

            colors_builder.init_meta_colors_builder(
                num_integers_in_metacolors + m_index.num_color_classes(), num_lists,
                m_partition_size, m_num_lists_in_partition);

            std::vector<uint32_t> metacolors;
            metacolors.reserve(num_partitions);  // at most

            for (uint64_t color_class_id = 0; color_class_id != num_color_classes;
                 ++color_class_id) {
                metacolors.clear();
                partial_color.clear();

                partition_id = 0;

                auto hash_and_write_metacolor = [&]() {
                    auto hash = util::hash128(reinterpret_cast<char const*>(partial_color.data()),
                                              partial_color.size() * sizeof(uint32_t));
                    auto it = m_hashes[partition_id].find(hash);
                    assert(it != m_hashes[partition_id].cend());  // must be found
                    uint32_t metacolor = (*it).second;
                    metacolors.push_back(metacolor);
                    partial_color.clear();
                };

                auto it = m_index.colors(color_class_id);
                uint64_t list_size = it.size();
                partition_endpoint curr_partition = partition_endpoints(0);

                permuted_list.clear();
                for (uint64_t i = 0; i != list_size; ++i, ++it) {
                    uint32_t ref_id = *it;
                    permuted_list.push_back(permutation[ref_id]);
                }
                std::sort(permuted_list.begin(), permuted_list.end());

                for (uint64_t i = 0; i != list_size; ++i) {
                    uint32_t ref_id = permuted_list[i];
                    while (ref_id >= curr_partition.end) {
                        if (!partial_color.empty()) hash_and_write_metacolor();
                        partition_id += 1;
                        curr_partition = partition_endpoints(partition_id);
                    }
                    assert(ref_id >= curr_partition.begin);
                    partial_color.push_back(ref_id - curr_partition.begin);
                }
                if (!partial_color.empty()) hash_and_write_metacolor();

                colors_builder.process_metacolors(metacolors.data(), metacolors.size());
            }

            colors_builder.build(idx.m_ccs);

            timer.stop();
            std::cout << "** building and compressing colors/meta-colors took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 3. copy m_u2c, m_k2u");
            timer.start();
            idx.m_u2c = m_index.get_u2c();
            idx.m_k2u = m_index.get_dict();
            timer.stop();
            std::cout << "** copying m_u2c and m_k2u took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 4. permuting m_filenames");
            timer.start();
            std::vector<std::string> filenames;
            filenames.resize(num_docs);
            for (uint64_t i = 0; i != num_docs; ++i) {
                filenames[permutation[i]] = m_index.filename(i);
            }
            idx.m_filenames.build(filenames);
            timer.stop();
            std::cout << "** permuting m_filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        if (m_build_config.check) {
            essentials::logger("step 5. check correctness...");
            for (uint64_t color_class_id = 0; color_class_id != num_color_classes;
                 ++color_class_id) {
                auto it_exp = m_index.colors(color_class_id);
                auto it_got = idx.colors(color_class_id);
                uint64_t exp_size = it_exp.size();
                uint64_t got_size = it_got.size();

                if (exp_size != got_size) {
                    std::cout << "got colors list of size " << got_size << " but expected "
                              << exp_size << std::endl;
                    return;
                }

                permuted_list.clear();
                for (uint64_t i = 0; i != exp_size; ++i, ++it_exp) {
                    uint32_t ref_id = *it_exp;
                    permuted_list.push_back(permutation[ref_id]);
                }
                std::sort(permuted_list.begin(), permuted_list.end());

                for (uint64_t i = 0; i != got_size; ++i, ++it_got) {
                    if (permuted_list[i] != *it_got) {
                        std::cout << "got ref " << *it_got << " BUT expected " << permuted_list[i]
                                  << std::endl;
                        return;
                    }
                }
            }
            essentials::logger("DONE!");
        }
    }

private:
    build_configuration m_build_config;
    index_type m_index;
    std::vector<uint32_t> m_partition_size;
    std::vector<uint32_t> m_num_lists_in_partition;
    std::vector<std::unordered_map<__uint128_t, uint32_t>> m_hashes;  // (hash,id)

    partition_endpoint partition_endpoints(uint64_t partition_id) const {
        assert(partition_id + 1 < m_partition_size.size());
        return {m_partition_size[partition_id], m_partition_size[partition_id + 1]};
    }
};

}  // namespace fulgor

#pragma once

#include "index.hpp"
#include "build_util.hpp"

namespace fulgor {

struct partition_endpoint {
    uint64_t begin, end;
};

struct permuter {
    permuter(build_configuration const& build_config)
        : m_build_config(build_config), m_num_partitions(0), m_max_partition_size(0) {}

    void permute(index_type const& index) {
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            essentials::logger("step 2. build sketches");
            timer.start();
            constexpr uint64_t p = 10;  // use 2^p bytes per HLL sketch
            build_reference_sketches(index, p, m_build_config.tmp_dirname + "/sketches.bin");
            timer.stop();
            std::cout << "** building sketches took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 3. clustering sketches");
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
            constexpr uint64_t seed = 0;
            params.set_min_delta(min_delta);
            params.set_max_iteration(max_iteration);
            params.set_min_cluster_size(min_cluster_size);
            params.set_random_seed(seed);
            auto clustering_data = kmeans::kmeans_divisive(points.begin(), points.end(), params);

            timer.stop();
            std::cout << "** clustering sketches took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();

            m_num_partitions = clustering_data.num_clusters;

            m_partition_size.resize(m_num_partitions + 1, 0);
            for (auto c : clustering_data.clusters) m_partition_size[c] += 1;

            /* take prefix sums */
            uint64_t val = 0;
            for (auto& size : m_partition_size) {
                if (size > m_max_partition_size) m_max_partition_size = size;
                uint64_t tmp = size;
                size = val;
                val += tmp;
            }

            const uint64_t num_docs = index.num_docs();

            /* build permutation */
            auto counts = m_partition_size;  // copy
            m_permutation.resize(num_docs);
            assert(clustering_data.clusters.size() == num_docs);
            for (uint64_t i = 0; i != num_docs; ++i) {
                uint32_t cluster_id = clustering_data.clusters[i];
                m_permutation[i] = counts[cluster_id];
                counts[cluster_id] += 1;
            }

            /* permute filenames */
            m_filenames.resize(num_docs);
            for (uint64_t i = 0; i != num_docs; ++i) {
                m_filenames[m_permutation[i]] = index.filename(i);
            }
        }
    }

    partition_endpoint partition_endpoints(uint64_t partition_id) const {
        assert(partition_id + 1 < m_partition_size.size());
        return {m_partition_size[partition_id], m_partition_size[partition_id + 1]};
    }

    uint64_t num_partitions() const { return m_num_partitions; }
    uint64_t max_partition_size() const { return m_max_partition_size; }
    std::vector<uint32_t> permutation() const { return m_permutation; }
    std::vector<uint32_t> partition_size() const { return m_partition_size; }
    std::vector<std::string> filenames() const { return m_filenames; }

private:
    build_configuration m_build_config;
    uint64_t m_num_partitions;
    uint64_t m_max_partition_size;
    std::vector<uint32_t> m_permutation;
    std::vector<uint32_t> m_partition_size;
    std::vector<std::string> m_filenames;
};

template <typename ColorClasses>
struct index<ColorClasses>::meta_builder {
    meta_builder() {}

    meta_builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        index_type index;
        essentials::logger("step 1. loading index to be partitioned...");
        essentials::load(index, m_build_config.index_filename_to_partition.c_str());
        essentials::logger("DONE");

        const uint64_t num_docs = index.num_docs();
        const uint64_t num_color_classes = index.num_color_classes();

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        permuter p(m_build_config);
        p.permute(index);
        auto const& permutation = p.permutation();

        const uint64_t num_partitions = p.num_partitions();
        const uint64_t max_partition_size = p.max_partition_size();
        std::cout << "num_partitions = " << num_partitions << std::endl;
        std::cout << "max_partition_size = " << max_partition_size << std::endl;

        {
            essentials::logger("step 4. building partial/meta colors");
            timer.start();

            std::ofstream metacolors_out(m_build_config.tmp_dirname + "/metacolors.bin",
                                         std::ios::binary);
            if (!metacolors_out.is_open()) throw std::runtime_error("error in opening file");

            uint64_t num_integers_in_metacolors = 0;
            uint64_t num_partial_colors = 0;

            std::vector<uint32_t> partial_color;
            std::vector<uint32_t> permuted_list;
            partial_color.reserve(max_partition_size);
            permuted_list.reserve(num_docs);

            typename ColorClasses::builder colors_builder;

            colors_builder.init_colors_builder(num_docs, num_partitions);
            for (uint64_t partition_id = 0; partition_id != num_partitions; ++partition_id) {
                auto endpoints = p.partition_endpoints(partition_id);
                uint64_t num_docs_in_partition = endpoints.end - endpoints.begin;
                colors_builder.init_color_partition(partition_id, num_docs_in_partition);
            }

            uint64_t partition_id = 0;
            uint32_t meta_color_list_size = 0;

            std::vector<std::unordered_map<__uint128_t, uint32_t>> hashes;  // (hash, id)
            hashes.resize(num_partitions);

            auto hash_and_compress = [&]() {
                assert(!partial_color.empty());
                auto hash = util::hash128(reinterpret_cast<char const*>(partial_color.data()),
                                          partial_color.size() * sizeof(uint32_t));
                uint32_t partial_color_id = 0;
                auto it = hashes[partition_id].find(hash);
                if (it == hashes[partition_id].cend()) {  // new partial color
                    partial_color_id = hashes[partition_id].size();
                    hashes[partition_id].insert({hash, partial_color_id});
                    colors_builder.process_colors(partition_id, partial_color.data(),
                                                  partial_color.size());
                } else {
                    partial_color_id = (*it).second;
                }

                /*  write meta color: (partition_id, partial_color_id)
                    Note: at this stage, partial_color_id is relative
                          to its partition (is not global yet).
                */
                metacolors_out.write(reinterpret_cast<char const*>(&partition_id),
                                     sizeof(uint32_t));
                metacolors_out.write(reinterpret_cast<char const*>(&partial_color_id),
                                     sizeof(uint32_t));

                partial_color.clear();
                meta_color_list_size += 1;
            };

            for (uint64_t color_class_id = 0; color_class_id != num_color_classes;
                 ++color_class_id) {
                /* permute list */
                permuted_list.clear();
                auto it = index.colors(color_class_id);
                uint64_t list_size = it.size();
                for (uint64_t i = 0; i != list_size; ++i, ++it) {
                    uint32_t ref_id = *it;
                    permuted_list.push_back(permutation[ref_id]);
                }
                std::sort(permuted_list.begin(), permuted_list.end());

                /* partition list */
                meta_color_list_size = 0;
                partition_id = 0;
                partition_endpoint curr_partition = p.partition_endpoints(0);
                assert(partial_color.empty());

                /* reserve space to hold the size of the meta color list */
                metacolors_out.write(reinterpret_cast<char const*>(&meta_color_list_size),
                                     sizeof(uint32_t));

                for (uint64_t i = 0; i != list_size; ++i) {
                    uint32_t ref_id = permuted_list[i];
                    while (ref_id >= curr_partition.end) {
                        if (!partial_color.empty()) hash_and_compress();
                        partition_id += 1;
                        curr_partition = p.partition_endpoints(partition_id);
                    }
                    assert(ref_id >= curr_partition.begin);
                    partial_color.push_back(ref_id - curr_partition.begin);
                }
                if (!partial_color.empty()) hash_and_compress();

                num_integers_in_metacolors += meta_color_list_size;

                /* write size of meta color list */
                uint64_t current_pos = metacolors_out.tellp();
                uint64_t num_bytes_in_meta_color_list =
                    2 * meta_color_list_size * sizeof(uint32_t) + sizeof(uint32_t);
                assert(current_pos >= num_bytes_in_meta_color_list);
                uint64_t pos = current_pos - num_bytes_in_meta_color_list;
                metacolors_out.seekp(pos);
                metacolors_out.write(reinterpret_cast<char const*>(&meta_color_list_size),
                                     sizeof(uint32_t));
                metacolors_out.seekp(current_pos);
            }

            metacolors_out.close();

            std::vector<uint64_t> num_partial_colors_before;
            std::vector<uint32_t> num_lists_in_partition;
            num_partial_colors_before.reserve(num_partitions);
            num_lists_in_partition.reserve(num_partitions);
            num_partial_colors = 0;
            for (partition_id = 0; partition_id != num_partitions; ++partition_id) {
                num_partial_colors_before.push_back(num_partial_colors);
                uint64_t num_partial_colors_in_partition = hashes[partition_id].size();
                num_partial_colors += num_partial_colors_in_partition;
                num_lists_in_partition.push_back(num_partial_colors_in_partition);
                std::cout << "num_partial_colors_in_partition-" << partition_id << ": "
                          << num_partial_colors_in_partition << std::endl;
            }

            std::cout << "total num. partial colors = " << num_partial_colors << std::endl;

            colors_builder.init_meta_colors_builder(num_integers_in_metacolors + num_color_classes,
                                                    num_partial_colors, p.partition_size(),
                                                    num_lists_in_partition);

            std::vector<uint32_t> metacolors;
            metacolors.reserve(num_partitions);  // at most

            std::ifstream metacolors_in(m_build_config.tmp_dirname + "/metacolors.bin",
                                        std::ios::binary);
            if (!metacolors_in.is_open()) throw std::runtime_error("error in opening file");

            for (uint64_t color_class_id = 0; color_class_id != num_color_classes;
                 ++color_class_id) {
                assert(metacolors.empty());
                uint32_t meta_color_list_size = 0;
                metacolors_in.read(reinterpret_cast<char*>(&meta_color_list_size),
                                   sizeof(uint32_t));
                for (uint32_t i = 0; i != meta_color_list_size; ++i) {
                    uint32_t partition_id = 0;
                    uint32_t partial_color_id = 0;
                    metacolors_in.read(reinterpret_cast<char*>(&partition_id), sizeof(uint32_t));
                    metacolors_in.read(reinterpret_cast<char*>(&partial_color_id),
                                       sizeof(uint32_t));
                    /* transform the partial_color_id into a global id */
                    metacolors.push_back(partial_color_id +
                                         num_partial_colors_before[partition_id]);
                }
                colors_builder.process_metacolors(metacolors.data(), metacolors.size());
                metacolors.clear();
            }

            metacolors_in.close();
            std::remove((m_build_config.tmp_dirname + "/metacolors.bin").c_str());
            colors_builder.build(idx.m_ccs);

            timer.stop();
            std::cout << "** building partial/meta colors took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 5. copy u2c and k2u");
            timer.start();
            idx.m_u2c = index.get_u2c();
            idx.m_k2u = index.get_dict();
            timer.stop();
            std::cout << "** copying u2c and k2u took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 6. building filenames");
            timer.start();
            idx.m_filenames.build(p.filenames());
            timer.stop();
            std::cout << "** building filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        if (m_build_config.check) {
            essentials::logger("step 7. check correctness...");

            std::vector<uint32_t> permuted_list;
            permuted_list.reserve(num_docs);

            for (uint64_t color_class_id = 0; color_class_id != num_color_classes;
                 ++color_class_id) {
                auto it_exp = index.colors(color_class_id);
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
};

}  // namespace fulgor

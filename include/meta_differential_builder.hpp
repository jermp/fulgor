#pragma once

#include <map>
#include "index.hpp"
#include "build_util.hpp"
#include "color_classes/meta.hpp"

namespace fulgor {

struct md_diff_permuter {
    md_diff_permuter(build_configuration const& build_config)
        : m_build_config(build_config), m_num_partitions(0) {}

    void permute_meta(meta<hybrid> const& m, uint64_t num_docs) {
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            essentials::logger("step 2. build sketches");
            timer.start();

            constexpr uint64_t p = 10;
            build_differential_sketches_from_meta(
                m, num_docs, m.num_color_classes(), p, m_build_config.num_threads,
                m_build_config.tmp_dirname + "/sketches.bin");

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
            std::vector<uint64_t> color_ids;
            uint64_t num_bytes_per_point = 0;
            uint64_t num_points = 0;
            uint64_t num_docs = 0;
            in.read(reinterpret_cast<char*>(&num_bytes_per_point), sizeof(uint64_t));
            in.read(reinterpret_cast<char*>(&num_docs), sizeof(uint64_t));
            in.read(reinterpret_cast<char*>(&num_points), sizeof(uint64_t));
            m_num_docs = num_docs;
            points.resize(num_points, kmeans::point(num_bytes_per_point));
            color_ids.resize(num_points);
            for (uint64_t i = 0; i != num_points; ++i) {
                in.read(reinterpret_cast<char*>(&color_ids[i]), sizeof(uint64_t));
            }
            for (auto& point : points) {
                in.read(reinterpret_cast<char*>(point.data()), num_bytes_per_point);
            }
            in.close();

            std::remove((m_build_config.tmp_dirname + "/sketches.bin").c_str());

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

            timer.stop();
            std::cout << "** clustering sketches took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();

            m_num_partitions = clustering_data.num_clusters;
            m_partition_size.resize(m_num_partitions + 1, 0);
            for (auto c : clustering_data.clusters) { m_partition_size[c] += 1; }

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
            std::vector<uint32_t> permutation(num_color_classes);
            assert(clustering_data.clusters.size() == num_color_classes);
            for (uint64_t i = 0; i != num_color_classes; ++i) {
                uint64_t cluster_id = clustering_data.clusters[i];
                permutation[i] = clusters_pos[cluster_id];
                clusters_pos[cluster_id] += 1;
            }

            m_color_classes_ids.resize(num_color_classes);
            for (uint64_t i = 0; i != num_color_classes; ++i) {
                m_color_classes_ids[permutation[i]] = color_ids[i];
            }

            std::cout << "Computed " << m_num_partitions << " partitions\n";
            std::vector<uint32_t> distribution_stats(num_docs, 0);

            m_permutation.resize(num_color_classes);
            m_references.resize(m_num_partitions);
            std::vector<uint32_t> distribution(num_docs, 0);
            uint64_t cluster_size = 0;
            for (uint64_t color_id = 0, cluster_id = 0; color_id != num_color_classes + 1;
                 ++color_id, ++cluster_size) {
                if (color_id == m_partition_size[cluster_id + 1]) {
                    auto& reference = m_references[cluster_id];
                    for (uint32_t i = 0; i != num_docs; ++i) {
                        if (distribution[i] >= ceil(1. * cluster_size / 2.))
                            reference.emplace_back(i);
                    }
                    fill(distribution.begin(), distribution.end(), 0);
                    cluster_id++;
                    cluster_size = 0;
                    if (color_id == num_color_classes) break;
                }
                meta<hybrid>::forward_iterator it = m.colors(m_color_classes_ids[color_id]);
                uint64_t size = it.meta_color_list_size();
                while (size-- > 0) {
                    uint64_t val = it.partition_id();
                    distribution[val]++;
                    distribution_stats[val]++;
                    it.next_partition_id();
                }
                m_permutation[color_id] = {cluster_id, m_color_classes_ids[color_id]};
            }
            /*
            for (uint64_t i = 0; i < num_docs; i++) {
                cout << i << ": " << distribution_stats[i] << ", ";
                if ((i + 1) % 10 == 0) cout << endl;
            }
            */
            cout << "FINISHED PERMUTING;" << endl;
        }
    }

    void permute(hybrid const& h) {
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            essentials::logger("step 2. build sketches");
            timer.start();

            constexpr uint64_t p = 10;
            build_differential_sketches_from_hybrid(h, h.num_docs(), h.num_color_classes(), p,
                                                    m_build_config.num_threads,
                                                    m_build_config.tmp_dirname + "/sketches.bin");

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
            std::vector<uint64_t> color_ids;
            uint64_t num_bytes_per_point = 0;
            uint64_t num_points = 0;
            uint64_t num_docs = 0;
            in.read(reinterpret_cast<char*>(&num_bytes_per_point), sizeof(uint64_t));
            in.read(reinterpret_cast<char*>(&num_docs), sizeof(uint64_t));
            in.read(reinterpret_cast<char*>(&num_points), sizeof(uint64_t));
            m_num_docs = num_docs;
            points.resize(num_points, kmeans::point(num_bytes_per_point));
            color_ids.resize(num_points);
            for (uint64_t i = 0; i != num_points; ++i) {
                in.read(reinterpret_cast<char*>(&color_ids[i]), sizeof(uint64_t));
            }
            for (auto& point : points) {
                in.read(reinterpret_cast<char*>(point.data()), num_bytes_per_point);
            }
            in.close();

            std::remove((m_build_config.tmp_dirname + "/sketches.bin").c_str());

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

            timer.stop();
            std::cout << "** clustering sketches took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();

            m_num_partitions = clustering_data.num_clusters;
            m_partition_size.resize(m_num_partitions + 1, 0);
            for (auto c : clustering_data.clusters) { m_partition_size[c] += 1; }

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
            std::vector<uint32_t> permutation(num_color_classes);
            assert(clustering_data.clusters.size() == num_color_classes);
            for (uint64_t i = 0; i != num_color_classes; ++i) {
                uint64_t cluster_id = clustering_data.clusters[i];
                permutation[i] = clusters_pos[cluster_id];
                clusters_pos[cluster_id] += 1;
            }

            m_color_classes_ids.resize(num_color_classes);
            for (uint64_t i = 0; i != num_color_classes; ++i) {
                m_color_classes_ids[permutation[i]] = color_ids[i];
            }

            std::cout << "Computed " << m_num_partitions << " partitions\n";

            m_permutation.resize(num_color_classes);
            m_references.resize(m_num_partitions);
            std::vector<uint32_t> distribution(num_docs, 0);
            uint64_t cluster_size = 0;
            for (uint64_t color_id = 0, cluster_id = 0; color_id != num_color_classes + 1;
                 ++color_id, ++cluster_size) {
                if (color_id == m_partition_size[cluster_id + 1]) {
                    auto& reference = m_references[cluster_id];
                    for (uint32_t i = 0; i != num_docs; ++i) {
                        if (distribution[i] >= ceil(1. * cluster_size / 2.))
                            reference.emplace_back(i);
                    }
                    fill(distribution.begin(), distribution.end(), 0);
                    cluster_id++;
                    cluster_size = 0;
                    if (color_id == num_color_classes) break;
                }
                auto it = h.colors(m_color_classes_ids[color_id]);
                for (uint32_t i = 0; i != it.size(); ++i, ++it) { distribution[*it]++; }
                m_permutation[color_id] = {cluster_id, m_color_classes_ids[color_id]};
            }
        }
    }

    uint64_t num_partitions() const { return m_num_partitions; }
    uint64_t num_docs() const { return m_num_docs; }
    std::vector<std::pair<uint32_t, uint32_t>> permutation() const { return m_permutation; }
    std::vector<uint32_t> color_classes_ids() const { return m_color_classes_ids; }
    std::vector<std::vector<uint32_t>> references() const { return m_references; }

private:
    build_configuration m_build_config;
    uint64_t m_num_partitions;
    uint64_t m_num_docs;
    std::vector<std::pair<uint32_t, uint32_t>> m_permutation;
    std::vector<std::vector<uint32_t>> m_references;
    std::vector<uint32_t> m_partition_size;
    std::vector<uint32_t> m_color_classes_ids;
};

template <typename ColorClasses>
struct index<ColorClasses>::meta_differential_builder {
    meta_differential_builder() {}

    meta_differential_builder(build_configuration const& build_config)
        : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        meta_index_type meta_index;
        essentials::logger("step 1. loading index to be mega-partitioned");
        essentials::load(meta_index, m_build_config.index_filename_to_partition.c_str());
        essentials::logger("DONE");

        const uint64_t num_docs = meta_index.num_docs();
        const uint64_t num_color_classes = meta_index.num_color_classes();

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        uint64_t num_partitions = meta_index.get_color_classes().num_partitions();

        {
            essentials::logger("step 4. building partial/meta colors");

            std::vector<hybrid> pc = meta_index.get_color_classes().partial_colors();
            assert(pc.size() == num_partitions);

            for (uint64_t i = 0; i < num_partitions; i++) {
                md_diff_permuter dp(m_build_config);
                dp.permute(pc[i]);

                differential::builder diff_builder;
                diff_builder.init_colors_builder(dp.num_docs());

                auto const& permutation = dp.permutation();
                auto const& references = dp.references();

                for (auto& reference : references) { diff_builder.encode_reference(reference); }
                for (auto& [cluster_id, color_id] : permutation) {
                    auto it = pc[i].colors(color_id);
                    diff_builder.encode_list(
                        cluster_id, 
                        references[cluster_id],
                        it.size(), 
                        [&it]() -> void {++it;},
                        [&it]() -> uint64_t {return *it;}
                    );
                }
                differential d;
                diff_builder.build(d);
                // TODO: add to index
                d.print_stats();
            }

            {
                essentials::logger("step infty. build differential-meta colors");

                std::map<std::vector<uint64_t>, uint64_t> meta_partitions;
                std::vector<uint64_t> partition_map(num_color_classes);
                std::vector<std::vector<uint64_t>> partition_bases;
                uint64_t num_integers_in_metacolors = 0;
                uint64_t num_meta_partitions = 0;
                for (uint64_t color_id = 0; color_id < num_color_classes; color_id++) {
                    auto it = meta_index.get_color_classes().colors(color_id);
                    uint64_t size = it.meta_color_list_size();
                    num_integers_in_metacolors += size;

                    std::vector<uint64_t> partition_list(size);
                    for (uint64_t i = 0; i < size; ++i, it.next_partition_id()) {
                        partition_list[i] = it.partition_id();
                    }
                    if (meta_partitions.count(partition_list) == 0){
                        meta_partitions[partition_list] = num_meta_partitions++;
                        partition_bases.push_back(partition_list);
                    }
                    partition_map[color_id] = meta_partitions[partition_list];
                }


                pthash::compact_vector::builder partial_colors_ids_builder(
                    num_integers_in_metacolors + num_color_classes,
                    std::ceil(std::log2(meta_index.get_color_classes().num_max_lists_in_partition()))
                );
                for (uint64_t color_id = 0; color_id < num_color_classes; color_id++) {
                    auto it = meta_index.get_color_classes().colors(color_id);
                    partial_colors_ids_builder.push_back(it.meta_color_list_size());
                    for(uint64_t i = 0; i < it.meta_color_list_size(); i++, it.next_partition_id()){
                        it.update_partition();
                        partial_colors_ids_builder.push_back(it.meta_color() - it.num_lists_before());
                    }
                }
                pthash::compact_vector partial_colors_ids;
                partial_colors_ids_builder.build(partial_colors_ids);
                cout << "  PARTIAL COLORS IDS SIZE: " << partial_colors_ids.bytes() << endl; 

                // TODO: add to index
            }

            {
/*
                essentials::logger("step 7. check correctness...");

                for (uint64_t color_id = 0; color_id < num_color_classes; color_id++) {
                    uint64_t meta_list_start =
                        idx.m_ccs.m_meta_colors_offsets.access(permutation[color_id].second);
                    pthash::compact_vector::iterator exp_it =
                        idx.m_ccs.m_meta_colors.at(meta_list_start);
                    uint64_t exp_it_size = *exp_it;
                    auto res_it = d.colors(color_id);
                    if (res_it.size() != exp_it_size) {
                        cout << "Error while checking color " << color_id
                             << ", different sizes: expected " << exp_it_size << " but got "
                             << res_it.size() << ")\n";
                        continue;
                    }
                    for (uint64_t j = 0; j < exp_it_size; ++j, ++res_it) {
                        auto exp = *exp_it;
                        auto got = *res_it;
                        if (exp != got) {
                            cout << "Error while checking color " << color_id
                                 << ", mismatch at position " << j << ": expected " << exp
                                 << " but got " << got << std::endl;
                        }
                    }
                }
                
                cout << " META-COLORS DONE.\n";
*/
            }

            timer.stop();
            std::cout << "** building partial/meta colors took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 5. copy u2c and k2u");
            timer.start();
            idx.m_u2c = meta_index.get_u2c();
            idx.m_k2u = meta_index.get_k2u();
            timer.stop();
            std::cout << "** copying u2c and k2u took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 6. building filenames");
            timer.start();
            idx.m_filenames = meta_index.get_filenames();
            timer.stop();
            std::cout << "** building filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        if (m_build_config.check) {
            essentials::logger("step 7. check correctness...");

            essentials::logger("DONE!");
        }
    }

private:
    build_configuration m_build_config;
};

}  // namespace fulgor
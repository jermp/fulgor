#pragma once

#include "index.hpp"

namespace fulgor {
struct differential_permuter {
    differential_permuter(build_configuration const& build_config)
        : m_build_config(build_config), m_num_partitions(0) {}

    void permute(index_type const& index) {
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            essentials::logger("step 2. build sketches");
            timer.start();

            constexpr uint64_t p = 10;
            build_reference_sketches_partitioned(index, p, m_build_config.num_threads,
                                                 m_build_config.tmp_dirname + "/sketches.bin", 0.,
                                                 1.);

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
                auto it = index.colors(m_color_classes_ids[color_id]);
                for (uint32_t i = 0; i != it.size(); ++i, ++it) { distribution[*it]++; }
                m_permutation[color_id] = {cluster_id, m_color_classes_ids[color_id]};
            }
        }
    }

    uint64_t num_partitions() const { return m_num_partitions; }
    std::vector<std::pair<uint32_t, uint32_t>> permutation() const { return m_permutation; }
    std::vector<uint32_t> color_classes_ids() const { return m_color_classes_ids; }
    std::vector<std::vector<uint32_t>> references() const { return m_references; }

private:
    build_configuration m_build_config;
    uint64_t m_num_partitions;
    std::vector<std::pair<uint32_t, uint32_t>> m_permutation;
    std::vector<std::vector<uint32_t>> m_references;
    std::vector<uint32_t> m_partition_size;
    std::vector<uint32_t> m_color_classes_ids;
};

template <typename ColorClasses>
struct index<ColorClasses>::differential_builder {
    differential_builder() {}

    differential_builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        index_type index;
        essentials::logger("step 1. loading index to be differentiated...");
        essentials::load(index, m_build_config.index_filename_to_differentiate.c_str());
        essentials::logger("DONE");

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        differential_permuter p(m_build_config);
        p.permute(index);
        auto const& permutation = p.permutation();
        auto const& references = p.references();

        const uint64_t num_partitions = p.num_partitions();
        const uint64_t num_color_classes = index.num_color_classes();
        std::cout << "num_partitions = " << num_partitions << std::endl;

        {
            essentials::logger("step 4. building differential colors");
            timer.start();

            typename ColorClasses::builder colors_builder;
            colors_builder.init_colors_builder(index.num_docs());

            for (auto& reference : references) { colors_builder.encode_reference(reference); }
            for (auto& [cluster_id, color_id] : permutation) {
                colors_builder.encode_list(cluster_id, references[cluster_id],
                                           index.colors(color_id));
            }
            colors_builder.build(idx.m_ccs);
        }

        {
            essentials::logger("step 5. permute unitigs and rebuild sshash");

            const std::string permuted_unitigs_filename =
                m_build_config.tmp_dirname + "/permuted_unitigs.fa";
            std::ofstream out(permuted_unitigs_filename.c_str());
            if (!out.is_open()) throw std::runtime_error("cannot open output file");

            pthash::darray1 d;  // for select_1 on index.u2c
            d.build(index.get_u2c());

            const uint64_t num_unitigs = index.get_u2c().size();
            pthash::bit_vector_builder u2c_builder(num_unitigs+1, 0);

            auto const& dict = index.get_k2u();
            const uint64_t k = dict.k();

            uint64_t pos = 0;
            for (uint64_t new_color_id = 0; new_color_id != num_color_classes; ++new_color_id) {
                auto [_, old_color_id] = permutation[new_color_id];
                uint64_t old_unitig_id_end = num_unitigs;
                if (old_color_id < num_color_classes -1 ){
                    old_unitig_id_end = d.select(index.get_u2c(), old_color_id) + 1;
                }
                uint64_t old_unitig_id_begin = 0;
                if (old_color_id > 0) {
                    old_unitig_id_begin = d.select(index.get_u2c(), old_color_id - 1) + 1;
                }

                // num. unitigs that have the same color
                pos += old_unitig_id_end - old_unitig_id_begin;
                // cout << "[" << new_color_id << "] " << pos << "\n";
                assert(pos-1 < u2c_builder.size());

                u2c_builder.set(pos-1, 1);

                for (uint64_t i = old_unitig_id_begin; i != old_unitig_id_end; ++i) {
                    auto it = dict.at_contig_id(i);
                    out << ">\n";
                    auto [_, kmer] = it.next();
                    out << kmer;
                    while (it.has_next()) {
                        auto [_, kmer] = it.next();
                        out << kmer[k - 1];  // overlaps!
                    }
                    out << '\n';
                }
            }

            assert(pos == num_unitigs);
            out.close();
            idx.m_u2c.build(&u2c_builder);

            /* build a new sshash::dictionary on the permuted unitigs */
            sshash::build_configuration sshash_config;
            sshash_config.k = m_build_config.k;
            sshash_config.m = m_build_config.m;
            sshash_config.canonical_parsing = m_build_config.canonical_parsing;
            sshash_config.verbose = m_build_config.verbose;
            sshash_config.tmp_dirname = m_build_config.tmp_dirname;
            sshash_config.print();
            idx.m_k2u.build(permuted_unitigs_filename, sshash_config);
            assert(idx.get_k2u().size() == dict.size());
            try {  // remove unitig file
                std::remove(permuted_unitigs_filename.c_str());
            } catch (std::exception const& e) { std::cerr << e.what() << std::endl; }
        }

        {
            essentials::logger("step 6. building filenames");
            timer.start();
            idx.m_filenames = index.get_filenames();
            timer.stop();
            std::cout << "** building filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        if (m_build_config.check) {
            essentials::logger("step 7. check correctness...");

            for (uint64_t color_id = 0; color_id < num_color_classes; color_id++) {
                auto exp_it = index.colors(permutation[color_id].second);
                auto res_it = idx.colors(color_id);
                if (res_it.size() != exp_it.size()) {
                    cout << "Error while checking color " << color_id
                         << ", different sizes: expected " << exp_it.size() << " but got "
                         << res_it.size() << ")\n";
                    continue;
                }
                for (uint64_t j = 0; j < exp_it.size(); ++j, ++exp_it, ++res_it) {
                    auto exp = *exp_it;
                    auto got = *res_it;
                    if (exp != got) {
                        cout << "Error while checking color " << color_id
                             << ", mismatch at position " << j << ": expected " << exp
                             << " but got " << got << std::endl;
                    }
                }
            }
            cout << " COLORS DONE.\n";

            for (uint64_t unitig_id = 0; unitig_id < idx.m_k2u.num_contigs(); ++unitig_id) {
                auto it = idx.get_k2u().at_contig_id(unitig_id);
                while (it.has_next()) {
                    auto [_, kmer] = it.next();
                    uint64_t new_contig_id = idx.get_k2u().lookup_advanced(kmer.c_str()).contig_id;
                    if (new_contig_id != unitig_id) {
                        std::cout << "expected " << unitig_id << " but found " << new_contig_id
                                  << std::endl;
                        continue;
                    }
                    uint64_t old_contig_id =
                        index.get_k2u().lookup_advanced(kmer.c_str()).contig_id;

                    uint64_t new_color_id = idx.u2c(new_contig_id);
                    uint64_t old_color_id = index.u2c(old_contig_id);
                    // cout << "[" << new_contig_id << "] " << new_color_id << " -> " << old_color_id << " (got: "<< permutation[new_color_id].second <<")\n";

                    auto exp_it = index.colors(old_color_id);
                    auto res_it = idx.colors(new_color_id);
                    if (res_it.size() != exp_it.size()) {
                        std::cout << "Error while checking color " << new_color_id
                                  << ", different sizes: expected " << exp_it.size() << " but got "
                                  << res_it.size() << std::endl;
                        continue;
                    }
                    for (uint64_t j = 0; j < exp_it.size(); ++j, ++exp_it, ++res_it) {
                        auto exp = *exp_it;
                        auto got = *res_it;
                        if (exp != got) {
                            std::cout << "Error while checking color " << new_color_id
                                      << ", mismatch at position " << j << ": expected " << exp
                                      << " but got " << got << std::endl;
                        }
                    }
                }
            }
        }
    }

private:
    build_configuration m_build_config;
};
}  // namespace fulgor
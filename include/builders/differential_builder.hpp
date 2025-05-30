#pragma once

#include "include/index.hpp"

namespace fulgor {
struct differential_permuter {
    differential_permuter(build_configuration const& build_config)
        : m_build_config(build_config), m_num_partitions(0) {}

    template <typename Index>
    void permute(Index const& index) {
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        const std::vector<float> slices = {0, 0.25, 0.5, 0.75, 1};
        const uint64_t num_slices = slices.size() - 1;

        {
            essentials::logger("step 2. build sketches");

            constexpr uint64_t p = 10;
            for (uint64_t slice_id = 0; slice_id != num_slices; slice_id++) {
                timer.start();
                build_colors_sketches_sliced<hybrid::forward_iterator>(
                    index.num_colors(), index.num_color_sets(),
                    [&](uint64_t color_set_id) -> hybrid::forward_iterator {
                        return index.color_set(color_set_id);
                    },
                    p, m_build_config.num_threads,
                    m_build_config.tmp_dirname + "/sketches" + std::to_string(slice_id) + ".bin",
                    slices[slice_id], slices[slice_id + 1]);
                timer.stop();
                std::cout << "** building sketches took " << timer.elapsed() << " seconds / "
                          << timer.elapsed() / 60 << " minutes" << std::endl;
                timer.reset();
            }
        }

        {
            essentials::logger("step 3. clustering sketches");

            std::vector<uint64_t> color_set_ids;
            std::vector<kmeans::cluster_data> clustering_data(num_slices);
            std::vector<uint64_t> num_points(num_slices);

            for (uint64_t slice_id = 0; slice_id < num_slices; slice_id++) {
                num_points[slice_id] = cluster("/sketches" + std::to_string(slice_id) + ".bin",
                                               clustering_data[slice_id], color_set_ids);
            }

            timer.start();

            m_num_partitions = 0;
            for (uint64_t slice_id = 0; slice_id < num_slices; ++slice_id) {
                if (num_points[slice_id] == 0) continue;
                m_num_partitions += clustering_data[slice_id].num_clusters;
            }

            m_partition_size.resize(m_num_partitions + 1, 0);

            uint64_t prev_num_clusters = 0;
            for (uint64_t slice_id = 0; slice_id < num_slices; slice_id++) {
                if (num_points[slice_id] == 0) continue;
                for (auto c : clustering_data[slice_id].clusters) {
                    m_partition_size[c + prev_num_clusters] += 1;
                }
                prev_num_clusters += clustering_data[slice_id].num_clusters;
            }

            /* prefix sum */
            {
                uint64_t val = 0;
                for (auto& size : m_partition_size) {
                    uint64_t tmp = size;
                    size = val;
                    val += tmp;
                }
            }

            const uint64_t num_color_sets = index.num_color_sets();
            m_num_colors = index.num_colors();
            std::vector<uint32_t> color_sets_ids;
            color_sets_ids.resize(num_color_sets);

            auto clusters_pos = m_partition_size;
            uint64_t clusters_size = 0;
            for (uint64_t slice_id = 0; slice_id < num_slices; ++slice_id) {
                if (num_points[slice_id] == 0) continue;
                clusters_size += clustering_data[slice_id].clusters.size();
            }
            assert(clusters_size == num_color_sets);
            (void)clusters_size;

            std::vector<uint32_t> permutation(num_color_sets);
            prev_num_clusters = 0;
            uint64_t prev_num_color_sets = 0;
            for (uint64_t slice_id = 0; slice_id < num_slices; slice_id++) {
                uint64_t num_cc = num_points[slice_id];
                for (uint64_t color_set_id = 0; color_set_id != num_cc; ++color_set_id) {
                    uint64_t cluster_id =
                        clustering_data[slice_id].clusters[color_set_id] + prev_num_clusters;
                    permutation[color_set_id + prev_num_color_sets] = clusters_pos[cluster_id];
                    clusters_pos[cluster_id] += 1;
                }
                prev_num_clusters += clustering_data[slice_id].num_clusters;
                prev_num_color_sets += num_cc;
            }

            for (uint64_t i = 0; i != num_color_sets; ++i) {
                color_sets_ids[permutation[i]] = color_set_ids[i];
            }

            std::cout << "Computed " << m_num_partitions << " partitions\n";

            m_permutation.resize(num_color_sets);
            for (uint64_t i = 0, cluster_id = 0; i != num_color_sets + 1; ++i) {
                if (i == m_partition_size[cluster_id + 1]) {
                    cluster_id++;
                    if (i == num_color_sets) break;
                }
                m_permutation[i] = {cluster_id, color_sets_ids[i]};
            }

            timer.stop();
            std::cout << "  ** OTHER operations took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
        }
    }

    uint64_t num_partitions() const { return m_num_partitions; }
    uint64_t num_colors() const { return m_num_colors; }
    std::vector<std::pair<uint32_t, uint32_t>> permutation() const { return m_permutation; }

private:
    build_configuration m_build_config;
    uint64_t m_num_partitions, m_num_colors;
    std::vector<std::pair<uint32_t, uint32_t>> m_permutation;
    std::vector<uint32_t> m_partition_size;

    uint64_t cluster(std::string filename, kmeans::cluster_data& clustering_data,
                     std::vector<uint64_t>& color_set_ids) {
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        timer.start();

        std::ifstream in(m_build_config.tmp_dirname + filename, std::ios::binary);
        if (!in.is_open()) throw std::runtime_error("error in opening file");

        std::vector<kmeans::point> points;
        std::vector<uint64_t> group_color_set_ids;
        uint64_t num_bytes_per_point = 0;
        uint64_t num_points = 0;
        uint64_t num_colors = 0;
        in.read(reinterpret_cast<char*>(&num_bytes_per_point), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&num_colors), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&num_points), sizeof(uint64_t));
        points.resize(num_points, kmeans::point(num_bytes_per_point));
        group_color_set_ids.resize(num_points);
        for (uint64_t i = 0; i != num_points; ++i) {
            in.read(reinterpret_cast<char*>(&group_color_set_ids[i]), sizeof(uint64_t));
        }
        for (auto& point : points) {
            in.read(reinterpret_cast<char*>(point.data()), num_bytes_per_point);
        }
        in.close();
        std::remove((m_build_config.tmp_dirname + filename).c_str());

        {
            essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> timer;
            timer.start();
            if (num_points == 0) {
                std::cout << "Found empty partition" << endl;
                clustering_data.num_clusters = 0;
                clustering_data.clusters = {};
                return 0;
            } else {
                kmeans::clustering_parameters params;
                float min_delta = 0.0001;
                float max_iteration = 10;
                uint64_t min_cluster_size = 0;
                uint64_t seed = 0;
                params.set_min_delta(min_delta);
                params.set_max_iteration(max_iteration);
                params.set_min_cluster_size(min_cluster_size);
                params.set_random_seed(seed);
                params.set_num_threads(m_build_config.num_threads);
                clustering_data = kmeans::kmeans_divisive(points.begin(), points.end(), params);
            }
            timer.stop();
            std::cout << "  ** ONLY clustering sketches took " << timer.elapsed() / 1000
                      << " seconds / " << timer.elapsed() / (60 * 1000) << " minutes" << std::endl;
        }

        color_set_ids.insert(color_set_ids.end(), group_color_set_ids.begin(),
                             group_color_set_ids.end());

        timer.stop();
        std::cout << "** clustering sketches took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;

        return num_points;
    }
};

template <typename ColorSets>
struct index<ColorSets>::differential_builder {
    differential_builder() {}

    differential_builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        const uint32_t num_threads = m_build_config.num_threads;

        index_type index;
        essentials::logger("step 1. loading index to be differentiated...");
        essentials::load(index, m_build_config.index_filename_to_partition.c_str());
        essentials::logger("DONE");

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        differential_permuter p(m_build_config);
        p.permute(index);
        auto const& permutation = p.permutation();
        const uint64_t num_partitions = p.num_partitions();
        const uint64_t num_color_sets = index.num_color_sets();
        std::cout << "num_partitions = " << num_partitions << std::endl;

        {
            essentials::logger("step 4. building differential color sets");
            timer.start();

            struct slice {
                uint64_t begin, end;
            };
            std::vector<slice> thread_slices;

            uint64_t load = 0;
            for(uint32_t color_set_id = 0; color_set_id < num_color_sets; ++color_set_id){
                load += index.color_set(color_set_id).size();
            }
            const uint64_t load_per_thread = load/num_threads;

            slice s = {0, 0};
            uint64_t curr_load = 0;
            uint64_t prev_cluster = permutation[0].first;
            for (uint64_t i = 0; i < num_color_sets; ++i){
                auto& [cluster_id, color_set_id] = permutation[i];
                if (cluster_id != prev_cluster && curr_load >= load_per_thread){
                    s.end = i;
                    thread_slices.push_back(s);
                    s.begin = i;
                    curr_load = 0;
                }
                curr_load += index.color_set(color_set_id).size();
            }
            s.end = num_color_sets;
            thread_slices.push_back(s);
            
            std::vector<typename ColorSets::builder> thread_builders(thread_slices.size(), m_build_config.num_colors);
            std::vector<std::thread> threads(thread_slices.size());
            
            auto encode_color_sets = [&](uint64_t thread_id) {
                auto& color_sets_builder = thread_builders[thread_id];
                auto& [begin, end] = thread_slices[thread_id];
                color_sets_builder.reserve_num_bits(16 * essentials::GB * 8);

                std::vector<uint32_t> group_endpoints;
                uint64_t curr_group = permutation[begin].first + 1; // different from first group
                for (uint64_t i = begin; i < end; i++){
                    auto& [group_id, color_set_id] = permutation[i];
                    if (group_id != curr_group) {
                        group_endpoints.push_back(i);
                    }
                }
                group_endpoints.push_back(end);
                
                std::vector<uint32_t> distribution(num_color_sets, 0);
                for (uint64_t group = 0; group < group_endpoints.size()-1; ++group) {
                    uint64_t g_begin = group_endpoints[group];
                    uint64_t g_end = group_endpoints[group+1];
                    std::vector<uint32_t> representative;

                    for (uint64_t i = g_begin; i < g_end; ++i) {
                        auto& [group_id, color_set_id] = permutation[i];
                        auto it = index.color_set(color_set_id);
                        for(; *it != num_color_sets; ++it) distribution[*it]++;
                    }
                    uint64_t g_size = g_end - g_begin;
                    for (uint64_t color = 0; color < num_color_sets; ++color){
                        if (distribution[color] >= ceil(1. * g_size / 2.)) representative.push_back(color);
                    }
                    color_sets_builder.process_partition(representative);

                    for (uint64_t i = g_begin; i < g_end; ++i) {
                        auto& [group_id, color_set_id] = permutation[i];
                        auto it = index.color_set(color_set_id);
                        color_sets_builder.process_color_set(it);
                    }
                    std::fill(distribution.begin(), distribution.end(), 0);
                }
            };

            for (uint64_t thread_id = 0; thread_id < thread_slices.size(); thread_id++){
                threads[thread_id] = std::thread(encode_color_sets, thread_id);
            }
            for (auto& t : threads){
                if (t.joinable()) t.join();
            }

            for (uint64_t thread_id = 1; thread_id < thread_builders.size(); thread_id++){
                thread_builders[0].append(thread_builders[thread_id]);
            }
            thread_builders[0].build(idx.m_color_sets);

            timer.stop();
            std::cout << "** building color sets took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 5. permute unitigs and rebuild k2u");
            timer.start();

            const std::string permuted_unitigs_filename =
                m_build_config.tmp_dirname + "/permuted_unitigs.fa";
            std::ofstream out(permuted_unitigs_filename.c_str());
            if (!out.is_open()) throw std::runtime_error("cannot open output file");

            auto const& u2c = index.get_u2c();
            bits::darray1 d;  // for select_1 on u2c
            d.build(u2c);

            const uint64_t num_unitigs = u2c.num_bits();
            bits::bit_vector::builder u2c_builder(num_unitigs + 1, 0);

            auto const& dict = index.get_k2u();
            const uint64_t k = dict.k();

            uint64_t pos = 0;
            for (uint64_t new_color_set_id = 0; new_color_set_id != num_color_sets;
                 ++new_color_set_id) {
                auto [_, old_color_set_id] = permutation[new_color_set_id];
                uint64_t old_unitig_id_end = num_unitigs;
                if (old_color_set_id < num_color_sets - 1) {
                    old_unitig_id_end = d.select(u2c, old_color_set_id) + 1;
                }
                uint64_t old_unitig_id_begin = 0;
                if (old_color_set_id > 0)
                    old_unitig_id_begin = d.select(u2c, old_color_set_id - 1) + 1;

                // num. unitigs that have the same color
                pos += old_unitig_id_end - old_unitig_id_begin;
                assert(pos - 1 < u2c_builder.num_bits());

                u2c_builder.set(pos - 1, 1);

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
            u2c_builder.build(idx.m_u2c);
            idx.m_u2c_rank1_index.build(idx.m_u2c);

            /* build a new sshash::dictionary on the permuted unitigs */
            sshash::build_configuration sshash_config;
            sshash_config.k = dict.k();
            sshash_config.m = dict.m();
            sshash_config.canonical_parsing = dict.canonicalized();
            sshash_config.verbose = m_build_config.verbose;
            sshash_config.tmp_dirname = m_build_config.tmp_dirname;
            sshash_config.num_threads = m_build_config.num_threads;
            sshash_config.print();
            idx.m_k2u.build(permuted_unitigs_filename, sshash_config);
            assert(idx.get_k2u().size() == dict.size());
            try {  // remove unitig file
                std::remove(permuted_unitigs_filename.c_str());
            } catch (std::exception const& e) { std::cerr << e.what() << std::endl; }

            timer.stop();
            std::cout << "** permuting unitigs and rebuilding k2u took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
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

            for (uint64_t color_set_id = 0; color_set_id < num_color_sets; color_set_id++) {
                auto exp_it = index.color_set(permutation[color_set_id].second);
                auto res_it = idx.color_set(color_set_id);
                if (res_it.size() != exp_it.size()) {
                    std::cout << "Error while checking color " << color_set_id
                              << ", different sizes: expected " << exp_it.size() << " but got "
                              << res_it.size() << ")" << std::endl;
                    continue;
                }

                for (uint64_t j = 0; j < exp_it.size(); ++j, ++exp_it, ++res_it) {
                    auto exp = *exp_it;
                    auto got = *res_it;
                    if (exp != got) {
                        std::cout << "Error while checking color " << color_set_id
                                  << ", mismatch at position " << j << ": expected " << exp
                                  << " but got " << got << std::endl;
                    }
                }
            }

            std::cout << " COLORS DONE." << std::endl;

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

                    uint64_t new_color_set_id = idx.u2c(new_contig_id);
                    uint64_t old_color_set_id = index.u2c(old_contig_id);

                    auto exp_it = index.color_set(old_color_set_id);
                    auto res_it = idx.color_set(new_color_set_id);
                    if (res_it.size() != exp_it.size()) {
                        std::cout << "Error while checking color " << new_color_set_id
                                  << ", different sizes: expected " << exp_it.size() << " but got "
                                  << res_it.size() << std::endl;
                        continue;
                    }
                    for (uint64_t j = 0; j < exp_it.size(); ++j, ++exp_it, ++res_it) {
                        auto exp = *exp_it;
                        auto got = *res_it;
                        if (exp != got) {
                            std::cout << "Error while checking color " << new_color_set_id
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

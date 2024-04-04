#pragma once

#include "index.hpp"
#include "GGCAT.hpp"
#include "color_list_sketcher.hpp"

namespace fulgor {

template <typename ColorClasses>
struct index<ColorClasses>::builder {
    builder() {}

    builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            essentials::logger("step 1. build colored compacted dBG");
            timer.start();
            m_ccdbg.build(m_build_config);
            m_build_config.num_docs = m_ccdbg.num_docs();
            timer.stop();
            std::cout << "** building the ccdBG took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 2. build m_u2c and m_ccs");
            timer.start();

            uint64_t num_unitigs = 0;
            uint64_t num_distinct_colors = 0;

            pthash::bit_vector_builder u2c_builder;

            /* write unitigs to fasta file for SSHash */
            std::ofstream out((m_build_config.file_base_name + ".fa").c_str());
            if (!out.is_open()) throw std::runtime_error("cannot open output file");

            typename ColorClasses::builder colors_builder(m_build_config.num_docs);

            m_ccdbg.loop_through_unitigs([&](ggcat::Slice<char> const unitig,
                                             ggcat::Slice<uint32_t> const colors, bool same_color) {
                try {
                    if (!same_color) {
                        num_distinct_colors += 1;
                        if (num_unitigs > 0) u2c_builder.set(num_unitigs - 1, 1);

                        /* compress colors */
                        colors_builder.process(colors.data, colors.size);
                    }
                    u2c_builder.push_back(0);

                    out << ">\n";
                    out.write(unitig.data, unitig.size);
                    out << '\n';

                    num_unitigs += 1;

                } catch (std::exception const& e) {
                    std::cerr << e.what() << std::endl;
                    exit(1);
                }
            });

            out.close();

            assert(num_unitigs > 0);
            assert(num_unitigs < (uint64_t(1) << 32));

            std::cout << "num_unitigs " << num_unitigs << std::endl;
            std::cout << "num_distinct_colors " << num_distinct_colors << std::endl;

            u2c_builder.set(num_unitigs - 1, 1);
            idx.m_u2c.build(&u2c_builder);
            assert(idx.m_u2c.size() == num_unitigs);
            assert(idx.m_u2c.num_ones() == num_distinct_colors);

            std::cout << "m_u2c.size() " << idx.m_u2c.size() << std::endl;
            std::cout << "m_u2c.num_ones() " << idx.m_u2c.num_ones() << std::endl;
            std::cout << "m_u2c.num_zeros() " << idx.m_u2c.num_zeros() << std::endl;

            colors_builder.build(idx.m_ccs);

            timer.stop();
            std::cout << "** building m_u2c and m_ccs took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 3. build m_k2u");
            timer.start();

            sshash::build_configuration sshash_config;
            sshash_config.k = m_build_config.k;
            sshash_config.m = m_build_config.m;
            sshash_config.canonical_parsing = m_build_config.canonical_parsing;
            sshash_config.verbose = m_build_config.verbose;
            sshash_config.tmp_dirname = m_build_config.tmp_dirname;
            sshash_config.print();
            idx.m_k2u.build(m_build_config.file_base_name + ".fa", sshash_config);
            try {  // remove unitig file
                std::remove((m_build_config.file_base_name + ".fa").c_str());
            } catch (std::exception const& e) { std::cerr << e.what() << std::endl; }

            timer.stop();
            std::cout << "** building m_k2u took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 3 + 1/2. cluster m_ccs");
            timer.start();

            std::string filename = "/m_ccs_sketches.bin";

            constexpr uint64_t p = 10;
            build_color_list_sketches(idx, p, m_build_config.num_threads,
                                         m_build_config.tmp_dirname + filename);

            std::ifstream in(m_build_config.tmp_dirname + filename, std::ios::binary);
            if (!in.is_open()) throw std::runtime_error("error in opening file");

            std::vector<kmeans::point> points;
            uint64_t num_bytes_per_point;
            uint64_t num_points;
            in.read(reinterpret_cast<char*>(&num_bytes_per_point), sizeof(uint64_t));
            in.read(reinterpret_cast<char*>(&num_points), sizeof(uint64_t));
            points.resize(num_points, kmeans::point(num_bytes_per_point));
            for(auto & point: points){
                in.read(reinterpret_cast<char*>(point.data()), num_bytes_per_point);
            }
            in.close();
            std::remove((m_build_config.tmp_dirname + filename).c_str());

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
            std::cout << "** clustering took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();

            const uint64_t m_num_partitions = clustering_data.num_clusters;
            std::vector<uint64_t> m_partition_size(m_num_partitions+1, 0);
            for (auto c : clustering_data.clusters) m_partition_size[c] += 1;

            /* prefix sum */
            {
                uint64_t val = 0;
                for (auto& size: m_partition_size){
                    uint64_t tmp = size;
                    size = val;
                    val += tmp;
                }
            }

            const uint64_t num_color_classes = idx.num_color_classes();

            auto clusters_pos = m_partition_size;
            std::vector<uint32_t> m_permutation(num_color_classes);
            assert(clustering_data.clusters.size() == num_color_classes);
            for(uint64_t i = 0; i != num_color_classes; ++i){
                uint32_t cluster_id = clustering_data.clusters[i];
                m_permutation[i] = clusters_pos[cluster_id];
                clusters_pos[cluster_id] += 1;
            }

            std::vector<uint64_t> m_color_classes_ids(num_color_classes);
            for(uint64_t i = 0; i != num_color_classes; ++i){
                m_color_classes_ids[m_permutation[i]] = i;
            }

            std::cout << "Computed " << m_num_partitions << " partitions\n";

            std::ofstream cluster_dump(m_build_config.file_base_name + ".clst", std::ios::binary);
            uint64_t num_docs = idx.num_docs();
            cluster_dump.write(reinterpret_cast<char const*>(&num_docs), 8);
            cluster_dump.write(reinterpret_cast<char const*>(&num_color_classes), 8);
            cluster_dump.write(reinterpret_cast<char const*>(&m_num_partitions), 8);

            for(uint64_t i = 1; i < m_partition_size.size(); ++i){
                cluster_dump.write(reinterpret_cast<char const*>(&m_partition_size[i]), 8);
            }

            for (uint64_t i = 0; i != num_color_classes; i++) {
                auto it = idx.colors(m_color_classes_ids[i]);
                const uint64_t size = it.size();
                cluster_dump.write(reinterpret_cast<char const*>(&size), 8);
                for (uint64_t j = 0; j < size; ++j, ++it) {
                    uint64_t val = *it;
                    cluster_dump.write(reinterpret_cast<char const*>(&val), 8);
                }
            }
            cluster_dump.close();
        }

        {
            essentials::logger("step 4. write filenames");
            timer.start();
            idx.m_filenames.build(m_ccdbg.filenames());
            timer.stop();
            std::cout << "** writing filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        if (m_build_config.check)  //
        {
            essentials::logger("step 5. check correctness...");
            m_ccdbg.loop_through_unitigs(
                [&](ggcat::Slice<char> const unitig, ggcat::Slice<uint32_t> const colors,
                    bool /* same_color */)  //
                {
                    auto lookup_result = idx.m_k2u.lookup_advanced(unitig.data);
                    const uint64_t unitig_id = lookup_result.contig_id;
                    const uint64_t color_id = idx.u2c(unitig_id);
                    for (uint64_t i = 1; i != unitig.size - idx.m_k2u.k() + 1; ++i) {
                        uint64_t got = idx.m_k2u.lookup_advanced(unitig.data + i).contig_id;
                        if (got != unitig_id) {
                            std::cout << "got unitig_id " << got << " but expected " << unitig_id
                                      << std::endl;
                            return;
                        }
                    }
                    auto fwd_it = idx.m_ccs.colors(color_id);
                    const uint64_t size = fwd_it.size();
                    if (size != colors.size) {
                        std::cout << "got colors list of size " << size << " but expected "
                                  << colors.size << std::endl;
                        return;
                    }
                    for (uint64_t i = 0; i != size; ++i, ++fwd_it) {
                        uint32_t ref = *fwd_it;
                        if (ref != colors.data[i]) {
                            std::cout << "got ref " << ref << " but expected " << colors.data[i]
                                      << std::endl;
                            return;
                        }
                    }
                },
                m_build_config.num_threads  //
            );
            essentials::logger("DONE!");
        }
    }

private:
    build_configuration m_build_config;
    GGCAT m_ccdbg;
};

}  // namespace fulgor
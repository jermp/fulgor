#pragma once

#include <map>
#include "index.hpp"
#include "build_util.hpp"

namespace fulgor {

template <typename ColorSets>
struct index<ColorSets>::meta_differential_builder {
    meta_differential_builder() {}

    meta_differential_builder(build_configuration const& build_config)
        : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        meta_index_type meta_index;
        essentials::logger("step 1. loading index to be partitioned");
        essentials::load(meta_index, m_build_config.index_filename_to_partition.c_str());
        essentials::logger("DONE");

        const uint64_t num_color_sets = meta_index.num_color_sets();

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        uint64_t num_partitions = meta_index.get_color_sets().num_partitions();

        meta_differential::builder builder;
        builder.init(meta_index.num_colors(), num_partitions);

        std::vector<std::vector<uint64_t>> partial_permutations(num_partitions);

        {
            essentials::logger("step 2. building differential partial/meta colors");
            timer.start();

            std::vector<hybrid> pc = meta_index.get_color_sets().partial_colors();
            assert(pc.size() == num_partitions);

            for (uint64_t i = 0; i < num_partitions; i++) {
                std::cout << " Partition " << i << " / " << num_partitions << std::endl;
                differential_permuter dp(m_build_config);
                dp.permute(pc[i]);

                differential::builder diff_builder;
                diff_builder.init_colors_builder(dp.num_colors());

                auto const& permutation = dp.permutation();
                auto const& references = dp.references();

                partial_permutations[i].resize(permutation.size());
                uint64_t original_id = 0;

                for (auto& reference : references) {
                    diff_builder.encode_representative(reference);
                }
                for (auto& [cluster_id, color_id] : permutation) {
                    auto it = pc[i].color_set(color_id);
                    diff_builder.encode_list(
                        cluster_id, references[cluster_id], it.size(), [&it]() -> void { ++it; },
                        [&it]() -> uint64_t { return *it; });
                    partial_permutations[i][color_id] = original_id++;
                }
                differential d;
                diff_builder.build(d);
                builder.process_partition(d);
                d.print_stats();
            }

            timer.stop();
            std::cout << "** building partial/meta colors took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        std::vector<uint64_t> permutation(num_color_sets);

        {
            essentials::logger("step 5. build differential-meta colors");
            timer.start();

            std::vector<uint32_t> counts;
            std::map<std::vector<uint64_t>, uint64_t> meta_partitions;
            std::vector<std::vector<uint64_t>> partition_sets;
            std::vector<uint64_t> color_set_to_partition_set(num_color_sets);
            uint64_t num_partition_sets = 0;
            for (uint64_t color_id = 0; color_id < num_color_sets; color_id++) {
                auto it = meta_index.get_color_sets().color_set(color_id);
                uint64_t size = it.meta_color_list_size();

                std::vector<uint64_t> partition_list(size);
                for (uint64_t i = 0; i < size; ++i, it.next_partition_id()) {
                    partition_list[i] = it.partition_id();
                }
                if (meta_partitions.count(partition_list) == 0) {
                    meta_partitions[partition_list] = num_partition_sets++;
                    partition_sets.push_back(partition_list);
                    counts.push_back(0);
                }
                color_set_to_partition_set[color_id] = meta_partitions[partition_list];
                counts[meta_partitions[partition_list]]++;
            }

            builder.init_meta_color_partition_sets(num_partition_sets);
            for (uint64_t partition_set_id = 0; partition_set_id < num_partition_sets;
                 partition_set_id++) {
                builder.process_meta_color_partition_set(partition_sets[partition_set_id]);
            }

            std::vector<uint64_t> cum_sum = {0};
            uint64_t prev_val = 0;
            for (uint64_t count : counts) {
                prev_val += count;
                cum_sum.push_back(prev_val);
            }

            for (uint64_t color_id = 0; color_id < num_color_sets; color_id++) {
                permutation[cum_sum[color_set_to_partition_set[color_id]]++] = color_id;
            }
            for (uint64_t permuted_id = 0; permuted_id < num_color_sets; permuted_id++) {
                uint64_t original_color_id = permutation[permuted_id];
                uint64_t partition_set_id = color_set_to_partition_set[original_color_id];
                auto it = meta_index.get_color_sets().color_set(original_color_id);
                uint64_t size = it.meta_color_list_size();
                std::vector<uint64_t> relative_colors;
                relative_colors.reserve(size);

                for (uint64_t i = 0; i < size; i++, it.next_partition_id()) {
                    it.update_partition();
                    uint64_t partition_id = partition_sets[partition_set_id][i];
                    relative_colors.push_back(
                        partial_permutations[partition_id]
                                            [it.meta_color() - it.num_lists_before()]);
                }
                builder.process_metacolors(partition_set_id, partition_sets[partition_set_id],
                                           relative_colors);
            }

            builder.build(idx.m_color_sets);

            timer.stop();
            std::cout << "** building differential-meta colors took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 6. build u2c and k2u");
            timer.start();

            const std::string permuted_unitigs_filename =
                m_build_config.tmp_dirname + "/permuted_unitigs.fa";
            std::ofstream out(permuted_unitigs_filename.c_str());
            if (!out.is_open()) throw std::runtime_error("cannot open output file");

            pthash::darray1 d;  // for select_1 on index.u2c
            d.build(meta_index.get_u2c());

            const uint64_t num_unitigs = meta_index.get_u2c().size();
            pthash::bit_vector_builder u2c_builder(num_unitigs + 1, 0);

            auto const& dict = meta_index.get_k2u();
            const uint64_t k = dict.k();

            uint64_t pos = 0;
            for (uint64_t new_color_id = 0; new_color_id != num_color_sets; ++new_color_id) {
                uint64_t old_color_id = permutation[new_color_id];
                uint64_t old_unitig_id_end = num_unitigs;
                if (old_color_id < num_color_sets - 1) {
                    old_unitig_id_end = d.select(meta_index.get_u2c(), old_color_id) + 1;
                }
                uint64_t old_unitig_id_begin = 0;
                if (old_color_id > 0) {
                    old_unitig_id_begin = d.select(meta_index.get_u2c(), old_color_id - 1) + 1;
                }

                // num. unitigs that have the same color
                pos += old_unitig_id_end - old_unitig_id_begin;
                // cout << "[" << new_color_id << "] " << pos << "\n";
                assert(pos - 1 < u2c_builder.size());

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
            idx.m_u2c.build(&u2c_builder);

            /* build a new sshash::dictionary on the permuted unitigs */
            sshash::build_configuration sshash_config;
            sshash_config.k = dict.k();
            sshash_config.m = dict.m();
            sshash_config.canonical_parsing = dict.canonicalized();
            sshash_config.verbose = m_build_config.verbose;
            sshash_config.tmp_dirname = m_build_config.tmp_dirname;
            sshash_config.print();
            idx.m_k2u.build(permuted_unitigs_filename, sshash_config);
            assert(idx.get_k2u().size() == dict.size());
            try {  // remove unitig file
                std::remove(permuted_unitigs_filename.c_str());
            } catch (std::exception const& e) { std::cerr << e.what() << std::endl; }

            timer.stop();
            std::cout << "** building u2c and k2u took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 7. copying filenames");
            timer.start();
            idx.m_filenames = meta_index.get_filenames();
            timer.stop();
            std::cout << "** copying filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        if (m_build_config.check) {
            essentials::logger("step 8. check correctness...");
            timer.start();

            uint64_t slice_size = ceil(idx.m_k2u.num_contigs() / m_build_config.num_threads);

            auto exe = [&](uint64_t thread_id) {
                uint64_t l = slice_size * thread_id;
                uint64_t r = min(slice_size * (thread_id + 1), idx.m_k2u.num_contigs());

                for (uint64_t unitig_id = l; unitig_id < r; ++unitig_id) {
                    auto it = idx.get_k2u().at_contig_id(unitig_id);
                    while (it.has_next()) {
                        auto [_, kmer] = it.next();
                        uint64_t new_contig_id =
                            idx.get_k2u().lookup_advanced(kmer.c_str()).contig_id;
                        if (new_contig_id != unitig_id) {
                            std::cout << "expected " << unitig_id << " but found " << new_contig_id
                                      << std::endl;
                            continue;
                        }
                        uint64_t old_contig_id =
                            meta_index.get_k2u().lookup_advanced(kmer.c_str()).contig_id;

                        uint64_t new_color_id = idx.u2c(new_contig_id);
                        uint64_t old_color_id = meta_index.u2c(old_contig_id);

                        auto exp_it = meta_index.color_set(old_color_id);
                        auto res_it = idx.color_set(new_color_id);
                        if (res_it.size() != exp_it.size()) {
                            std::cout << "Error while checking color " << new_color_id
                                      << ", different sizes: expected " << exp_it.size()
                                      << " but got " << res_it.size() << std::endl;
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
            };

            std::vector<std::thread> threads(m_build_config.num_threads);
            for (uint64_t thread_id = 0; thread_id != m_build_config.num_threads; ++thread_id) {
                threads[thread_id] = std::thread(exe, thread_id);
            }
            for (auto& t : threads) {
                if (t.joinable()) t.join();
            }

            timer.stop();
            std::cout << "** checking correctness took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            essentials::logger("DONE!");
        }
    }

private:
    build_configuration m_build_config;
};

}  // namespace fulgor

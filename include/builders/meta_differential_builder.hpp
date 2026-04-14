#pragma once

#include "include/index.hpp"
#include "include/build_util.hpp"

namespace fulgor {

template <typename ColorSets>
struct index<ColorSets>::meta_differential_builder {
    meta_differential_builder() {}

    meta_differential_builder(build_configuration const& build_config)
        : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        essentials::logger("step 1. loading index to be partitioned");
        essentials::load(meta_index, m_build_config.index_filename_to_partition.c_str());
        essentials::logger("DONE");
        const uint32_t num_threads = m_build_config.num_threads;

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        const uint64_t num_partitions = meta_index.get_color_sets().num_partitions();
        const uint64_t num_color_sets = meta_index.num_color_sets();

        meta_differential::builder builder;
        builder.init(meta_index.num_colors(), num_partitions);

        std::vector<std::vector<uint32_t>> partial_permutations(num_partitions);

        {
            essentials::logger("step 2. building differential partial/meta color sets");
            timer.start();

            std::vector<hybrid> const& pc = meta_index.get_color_sets().partial_colors();
            assert(pc.size() == num_partitions);

            for (uint64_t meta_partition_id = 0; meta_partition_id < num_partitions;
                 meta_partition_id++) {
                std::cout << " Partition " << meta_partition_id << " / " << num_partitions - 1
                          << std::endl;
                auto& meta_partition = pc[meta_partition_id];
                const uint64_t num_partition_color_sets = meta_partition.num_color_sets();
                const uint64_t num_partition_colors = meta_partition.num_colors();

                differential_permuter dp(m_build_config);
                dp.permute(meta_partition);
                auto const& permutation = dp.permutation();
                partial_permutations[meta_partition_id].resize(num_partition_color_sets);

                essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
                timer.start();

                struct slice {
                    uint64_t begin, end;
                };
                std::vector<slice> thread_slices;

                uint64_t load = 0;
                for (uint32_t color_set_id = 0; color_set_id < num_partition_color_sets;
                     ++color_set_id) {
                    load += meta_partition.color_set(color_set_id).size();
                }
                const uint64_t load_per_thread = load / num_threads;

                slice s = {0, 0};
                uint64_t curr_load = 0;
                uint64_t prev_cluster = permutation[0].first;
                for (uint64_t i = 0; i < num_partition_color_sets; ++i) {
                    auto& [cluster_id, color_set_id] = permutation[i];
                    if (cluster_id != prev_cluster && curr_load >= load_per_thread) {
                        s.end = i;
                        thread_slices.push_back(s);
                        s.begin = i;
                        curr_load = 0;
                    }
                    curr_load += meta_partition.color_set(color_set_id).size();
                }
                s.end = num_partition_color_sets;
                thread_slices.push_back(s);

                std::vector<differential::builder> thread_builders(thread_slices.size(),
                                                                   num_partition_colors);
                std::vector<std::thread> threads(thread_slices.size());

                auto encode_color_sets = [&thread_builders, &thread_slices, &permutation,
                                          &meta_partition,
                                          num_partition_colors](uint64_t thread_id) {
                    auto& color_sets_builder = thread_builders[thread_id];
                    auto& [begin, end] = thread_slices[thread_id];
                    color_sets_builder.reserve_num_bits(16 * essentials::GB * 8);

                    std::vector<uint64_t> group_endpoints;
                    uint64_t curr_group =
                        permutation[begin].first + 1;  // different from first group
                    for (uint64_t i = begin; i < end; i++) {
                        auto& [group_id, color_set_id] = permutation[i];
                        if (group_id != curr_group) {
                            group_endpoints.push_back(i);
                            curr_group = group_id;
                        }
                    }
                    group_endpoints.push_back(end);

                    std::vector<uint32_t> distribution(num_partition_colors, 0);
                    for (uint64_t group = 0; group < group_endpoints.size() - 1; ++group) {
                        uint64_t g_begin = group_endpoints[group];
                        uint64_t g_end = group_endpoints[group + 1];
                        std::vector<uint32_t> representative;
                        representative.reserve(num_partition_colors);

                        for (uint64_t i = g_begin; i < g_end; ++i) {
                            auto& [group_id, color_set_id] = permutation[i];
                            auto it = meta_partition.color_set(color_set_id);
                            uint64_t it_size = it.size();
                            for (uint64_t pos = 0; pos < it_size; ++pos, ++it) {
                                distribution[*it]++;
                            }
                        }
                        uint64_t g_size = g_end - g_begin;
                        for (uint64_t color = 0; color < num_partition_colors; ++color) {
                            if (distribution[color] >= ceil(1. * g_size / 2.))
                                representative.push_back(color);
                        }
                        color_sets_builder.process_partition(representative);

                        for (uint64_t i = g_begin; i < g_end; ++i) {
                            auto& [group_id, color_set_id] = permutation[i];
                            auto it = meta_partition.color_set(color_set_id);
                            color_sets_builder.process_color_set(it);
                        }
                        std::fill(distribution.begin(), distribution.end(), 0);
                    }
                };

                for (uint64_t i = 0; i < num_partition_color_sets; i++) {
                    auto& [group_id, color_set_id] = permutation[i];
                    partial_permutations[meta_partition_id][color_set_id] = i;
                }

                for (uint64_t thread_id = 0; thread_id < thread_slices.size(); thread_id++) {
                    threads[thread_id] = std::thread(encode_color_sets, thread_id);
                }
                for (auto& t : threads) {
                    if (t.joinable()) t.join();
                }

                for (uint64_t thread_id = 1; thread_id < thread_builders.size(); thread_id++) {
                    thread_builders[0].append(thread_builders[thread_id]);
                }
                differential d;
                thread_builders[0].build(d);
                builder.process_partition(d);

                timer.stop();
                std::cout << "  ** building the color sets for partition " << meta_partition_id
                          << " took " << timer.elapsed() << " seconds / " << timer.elapsed() / 60
                          << " minutes" << std::endl;
            }

            timer.stop();
            std::cout << "** building partial/meta color sets took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        std::vector<uint32_t> permutation(num_color_sets);

        {
            essentials::logger("step 5. build differential-meta color sets");
            timer.start();

            std::vector<uint32_t> endpoints = {0};
            std::vector<meta<hybrid>::iterator_type> iterators;
            std::queue<std::pair<uint32_t, uint32_t>> slices;
            slices.emplace(0, num_color_sets);

            for (uint64_t color_set_id = 0; color_set_id < num_color_sets; color_set_id++) {
                permutation[color_set_id] = color_set_id;
                iterators.push_back(meta_index.color_set(color_set_id));
            }
            while (!slices.empty()) {
                auto [start, end] = slices.front();
                slices.pop();

                std::sort(permutation.begin() + start, permutation.begin() + end,
                          [&](const int a, const int b) {
                              auto& it_a = iterators[a];
                              auto& it_b = iterators[b];
                              uint32_t val_a = it_a.partition_id();
                              uint32_t val_b = it_b.partition_id();

                              return val_a < val_b;
                          });

                uint32_t prev = iterators[permutation[start]].partition_id();
                bool broken = false;

                if (prev == num_partitions) continue;

                iterators[permutation[start]].next_partition_id();
                for (uint64_t pos = start + 1; pos < end; pos++) {
                    uint32_t color_set_id = permutation[pos];
                    auto& it = iterators[color_set_id];
                    uint32_t partition_id = it.partition_id();

                    if (partition_id != prev) {
                        endpoints.push_back(pos);
                        slices.emplace(start, pos);
                        prev = partition_id;
                        start = pos;
                    }

                    if (partition_id == num_partitions) {
                        broken = true;
                        break;
                    }

                    it.next_partition_id();
                }

                if (!broken) { slices.emplace(start, end); }
            }

            std::sort(endpoints.begin(), endpoints.end());
            builder.init_meta_color_partition_sets(endpoints.size());

            uint64_t partition_set_id = 0;
            for (uint64_t permuted_id = 0; permuted_id < num_color_sets; permuted_id++) {
                uint64_t color_set_id = permutation[permuted_id];

                if (permuted_id == endpoints[partition_set_id]) {
                    auto it = meta_index.color_set(color_set_id);
                    uint32_t size = it.meta_color_set_size();
                    std::vector<uint32_t> partition_set(size);
                    for (uint32_t i = 0; i < size; ++i, it.next_partition_id()) {
                        partition_set[i] = it.partition_id();
                    }
                    builder.process_meta_color_partition_set(partition_set);
                    partition_set_id++;
                }

                auto it = meta_index.color_set(color_set_id);
                const uint64_t size = it.meta_color_set_size();
                std::vector<uint32_t> relative_colors;
                relative_colors.reserve(size);

                for (uint64_t i = 0; i < size; i++, it.next_partition_id()) {
                    it.update_partition();
                    uint64_t partition_id = it.partition_id();
                    relative_colors.push_back(
                        partial_permutations[partition_id]
                                            [it.meta_color() - it.num_color_sets_before()]);
                }
                builder.process_metacolor_set(relative_colors);
            }

            builder.build(idx.m_color_sets);

            timer.stop();
            std::cout << "** building differential-meta color sets took " << timer.elapsed()
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

            auto const& u2c = meta_index.get_u2c();
            bits::darray1 d;  // for select_1 on u2c
            d.build(u2c);

            const uint64_t num_unitigs = u2c.num_bits();
            bits::bit_vector::builder u2c_builder(num_unitigs + 1, 0);

            auto const& dict = meta_index.get_k2u();
            const uint64_t k = dict.k();

            uint64_t pos = 0;
            for (uint64_t new_color_id = 0; new_color_id != num_color_sets; ++new_color_id) {
                uint64_t old_color_id = permutation[new_color_id];
                uint64_t old_unitig_id_end = num_unitigs;
                if (old_color_id < num_color_sets - 1) {
                    old_unitig_id_end = d.select(u2c, old_color_id) + 1;
                }
                uint64_t old_unitig_id_begin = 0;
                if (old_color_id > 0) old_unitig_id_begin = d.select(u2c, old_color_id - 1) + 1;

                // num. unitigs that have the same color
                pos += old_unitig_id_end - old_unitig_id_begin;
                // cout << "[" << new_color_id << "] " << pos << "\n";
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
            assert(dict.canonical() == true);
            sshash_config.canonical = dict.canonical();
            sshash_config.verbose = m_build_config.verbose;
            sshash_config.tmp_dirname = m_build_config.tmp_dirname;
            sshash_config.num_threads = util::largest_power_of_2(m_build_config.num_threads);
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
    }

    void check(index& idx) {
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        essentials::logger("step 8. check correctness...");
        timer.start();
        std::atomic<uint64_t> num_checked_unitigs(0);

        auto exe = [this, &idx, &num_checked_unitigs](uint64_t start, uint64_t end) {
            assert(start < end);
            assert(end <= idx.m_k2u.num_contigs());

            uint64_t prev_new_color_set_id = idx.num_color_sets();
            uint64_t prev_old_color_set_id = idx.num_color_sets();

            for (uint64_t unitig_id = start; unitig_id < end; ++unitig_id) {
                auto it = idx.get_k2u().at_contig_id(unitig_id);
                while (it.has_next()) {
                    auto [_, kmer] = it.next();
                    uint64_t new_contig_id = idx.get_k2u().lookup_advanced(kmer.c_str()).contig_id;
                    if (new_contig_id != unitig_id) {
                        std::cout << "\033[1;31m"
                                  << "expected " << unitig_id << " but found " << new_contig_id
                                  << ")\033[0m" << std::endl;
                        continue;
                    }
                    uint64_t old_contig_id =
                        meta_index.get_k2u().lookup_advanced(kmer.c_str()).contig_id;

                    uint64_t new_color_set_id = idx.u2c(new_contig_id);
                    uint64_t old_color_set_id = meta_index.u2c(old_contig_id);

                    if (new_color_set_id == prev_new_color_set_id &&
                        old_color_set_id != prev_old_color_set_id) {
                        std::cout << "\033[1;31m"
                                  << "Error while checking unitig " << unitig_id
                                  << ", expected color set id " << prev_old_color_set_id
                                  << " but got " << old_color_set_id << ")\033[0m" << std::endl;
                        continue;
                    }

                    if (new_color_set_id == prev_new_color_set_id) continue;
                    prev_new_color_set_id = new_color_set_id;
                    prev_old_color_set_id = old_color_set_id;

                    auto exp_it = meta_index.color_set(old_color_set_id);
                    auto res_it = idx.color_set(new_color_set_id);
                    if (res_it.size() != exp_it.size()) {
                        std::cout << "\033[1;31m"
                                  << "Error while checking color " << new_color_set_id
                                  << ", different sizes: expected " << exp_it.size() << " but got "
                                  << res_it.size() << ")\033[0m" << std::endl;
                        continue;
                    }
                    for (uint64_t j = 0; j < exp_it.size(); ++j, ++exp_it, ++res_it) {
                        auto exp = *exp_it;
                        auto got = *res_it;
                        if (exp != got) {
                            std::cout << "\033[1;31m"
                                      << "Error while checking color " << new_color_set_id
                                      << ", mismatch at position " << j << ": expected " << exp
                                      << " but got " << got << ")\033[0m" << std::endl;
                        }
                    }
                }

                if (++num_checked_unitigs % 1000 == 0) {
                    std::cout << "\rChecked " << num_checked_unitigs << "/"
                              << idx.m_k2u.num_contigs() << " unitigs " << std::flush;
                }
            }
        };

        kmeans::thread_pool threads(m_build_config.num_threads);
        const uint64_t num_unitigs = idx.m_k2u.num_contigs();
        const uint64_t load_per_thread = num_unitigs / (m_build_config.num_threads << 10);
        uint64_t start = 0, end = load_per_thread;
        while (end < num_unitigs) {
            threads.enqueue([&, start, end] { exe(start, std::min(end, num_unitigs)); });
            start = end;
            end = std::min(end + load_per_thread, num_unitigs);
        }
        threads.enqueue([&, start, num_unitigs] { exe(start, num_unitigs); });  // last one
        threads.wait();

        std::cout << "\rChecked " << num_checked_unitigs << "/" << idx.m_k2u.num_contigs()
                  << " unitigs (" << std::endl;

        timer.stop();
        std::cout << "** checking correctness took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        essentials::logger("CHECKS DONE!");
    }

private:
    build_configuration m_build_config;
    mfur_index_t meta_index;
};

}  // namespace fulgor

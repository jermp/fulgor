#pragma once

#include "include/index.hpp"
#include "include/build_util.hpp"

#include <shared_mutex>

namespace fulgor {

struct partition_endpoint {
    uint64_t begin, end;  // [..)
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
            build_reference_sketches(index, p, m_build_config.num_threads,
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
            params.set_num_threads(m_build_config.num_threads);
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

            const uint64_t num_colors = index.num_colors();

            /* build permutation */
            auto counts = m_partition_size;  // copy
            m_permutation.resize(num_colors);
            assert(clustering_data.clusters.size() == num_colors);
            for (uint64_t i = 0; i != num_colors; ++i) {
                uint32_t cluster_id = clustering_data.clusters[i];
                m_permutation[i] = counts[cluster_id];
                counts[cluster_id] += 1;
            }

            /* permute filenames */
            m_filenames.resize(num_colors);
            for (uint64_t i = 0; i != num_colors; ++i) {
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

template <typename ColorSets>
struct index<ColorSets>::meta_builder {
    meta_builder() {}

    meta_builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        index_type index;
        essentials::logger("step 1. loading index to be partitioned...");
        essentials::load(index, m_build_config.index_filename_to_partition.c_str());
        essentials::logger("DONE");

        const uint64_t num_threads = m_build_config.num_threads;
        const uint64_t num_colors = index.num_colors();
        const uint64_t num_color_sets = index.num_color_sets();

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        permuter p(m_build_config);
        p.permute(index);
        auto const& permutation = p.permutation();

        const uint64_t num_partitions = p.num_partitions();
        const uint64_t max_partition_size = p.max_partition_size();
        std::cout << "num_partitions = " << num_partitions << std::endl;
        std::cout << "max_partition_size = " << max_partition_size << std::endl;

        {
            essentials::logger("step 4. building partial/meta color sets");
            timer.start();

            atomic_uint64_t num_integers_in_metacolor_sets = 0;
            uint64_t num_partial_color_sets = 0;

            typename ColorSets::builder color_sets_builder;

            color_sets_builder.init_color_sets_builder(num_colors, num_partitions);
            for (uint64_t partition_id = 0; partition_id != num_partitions; ++partition_id) {
                auto endpoints = p.partition_endpoints(partition_id);
                uint64_t num_colors_in_partition = endpoints.end - endpoints.begin;
                color_sets_builder.init_partition(partition_id, num_colors_in_partition);
                color_sets_builder.reserve_num_bits(partition_id, 8 * essentials::GB * 8);
            }

            std::vector<std::unordered_map<__uint128_t,            // key
                                           uint32_t,               // value
                                           util::hasher_uint128_t  // key's hasher
                                           >>
                hashes;  // (hash, id)
            hashes.resize(num_partitions);

            std::vector<std::thread> threads(num_threads);
            std::vector<uint32_t> thread_slices(num_threads + 1);
            std::vector<std::shared_mutex> partitions_mutex(num_partitions);

            for (uint64_t i = 0; i < num_threads; ++i) {
                thread_slices[i] = index.num_color_sets() / num_threads * i;
            }
            thread_slices[num_threads] = index.num_color_sets();

            auto exe = [&](uint64_t thread_id) {
                string tmp_filename = metacolor_set_file_name(thread_id);
                uint64_t partition_id = 0;
                uint32_t meta_color_set_size = 0;
                std::vector<uint32_t> partial_color_set;
                std::vector<uint32_t> permuted_set;
                std::ofstream metacolor_sets_ofstream(tmp_filename, std::ios::binary);
                if (!metacolor_sets_ofstream.is_open()) {
                    throw std::runtime_error("error in opening file");
                }

                partial_color_set.reserve(max_partition_size);
                permuted_set.reserve(num_colors);

                auto hash_and_compress = [&]() {
                    std::lock_guard<std::shared_mutex> lock(partitions_mutex[partition_id]);
                    assert(!partial_color_set.empty());
                    uint32_t partial_color_set_id = 0;
                    auto hash =
                        util::hash128(reinterpret_cast<char const*>(partial_color_set.data()),
                                      partial_color_set.size() * sizeof(uint32_t));
                    auto it = hashes[partition_id].find(hash);

                    if (it == hashes[partition_id].cend()) {  // new partial color
                        partial_color_set_id = hashes[partition_id].size();
                        hashes[partition_id].insert({hash, partial_color_set_id});
                        color_sets_builder.encode_color_set(partition_id, partial_color_set.data(),
                                                            partial_color_set.size());
                    } else {
                        partial_color_set_id = (*it).second;
                    }

                    /*  write meta color: (partition_id, partial_color_set_id)
                        Note: at this stage, partial_color_set_id is relative
                              to its partition (is not global yet).
                    */
                    metacolor_sets_ofstream.write(reinterpret_cast<char const*>(&partition_id),
                                                  sizeof(uint32_t));
                    metacolor_sets_ofstream.write(
                        reinterpret_cast<char const*>(&partial_color_set_id), sizeof(uint32_t));

                    partial_color_set.clear();
                    meta_color_set_size += 1;
                };

                for (uint64_t color_set_id = thread_slices[thread_id];
                     color_set_id != thread_slices[thread_id + 1]; ++color_set_id) {
                    /* permute set */
                    permuted_set.clear();
                    auto it = index.color_set(color_set_id);
                    uint64_t set_size = it.size();
                    for (uint64_t i = 0; i != set_size; ++i, ++it) {
                        uint32_t ref_id = *it;
                        permuted_set.push_back(permutation[ref_id]);
                    }
                    std::sort(permuted_set.begin(), permuted_set.end());

                    /* partition set */
                    meta_color_set_size = 0;
                    partition_id = 0;
                    partition_endpoint curr_partition = p.partition_endpoints(0);
                    assert(partial_color_set.empty());

                    /* reserve space to hold the size of the meta color set */
                    metacolor_sets_ofstream.write(
                        reinterpret_cast<char const*>(&meta_color_set_size), sizeof(uint32_t));

                    for (uint64_t i = 0; i != set_size; ++i) {
                        uint32_t ref_id = permuted_set[i];
                        while (ref_id >= curr_partition.end) {
                            if (!partial_color_set.empty()) hash_and_compress();
                            partition_id += 1;
                            curr_partition = p.partition_endpoints(partition_id);
                        }
                        assert(ref_id >= curr_partition.begin);
                        partial_color_set.push_back(ref_id - curr_partition.begin);
                    }
                    if (!partial_color_set.empty()) hash_and_compress();

                    num_integers_in_metacolor_sets += meta_color_set_size;

                    /* write size of meta color set */
                    uint64_t current_pos = metacolor_sets_ofstream.tellp();
                    uint64_t num_bytes_in_meta_color_set =
                        2 * meta_color_set_size * sizeof(uint32_t) + sizeof(uint32_t);
                    assert(current_pos >= num_bytes_in_meta_color_set);
                    uint64_t pos = current_pos - num_bytes_in_meta_color_set;
                    metacolor_sets_ofstream.seekp(pos);
                    metacolor_sets_ofstream.write(
                        reinterpret_cast<char const*>(&meta_color_set_size), sizeof(uint32_t));
                    metacolor_sets_ofstream.seekp(current_pos);
                }

                metacolor_sets_ofstream.close();
            };

            for (uint64_t i = 0; i != num_threads; ++i) threads[i] = std::thread(exe, i);

            for (auto& t : threads) {
                if (t.joinable()) t.join();
            }

            std::vector<uint64_t> num_partial_color_sets_before;
            std::vector<uint32_t> num_sets_in_partition;
            num_partial_color_sets_before.reserve(num_partitions);
            num_sets_in_partition.reserve(num_partitions);
            num_partial_color_sets = 0;
            for (uint64_t partition_id = 0; partition_id != num_partitions; ++partition_id) {
                num_partial_color_sets_before.push_back(num_partial_color_sets);
                uint64_t num_partial_color_sets_in_partition = hashes[partition_id].size();
                num_partial_color_sets += num_partial_color_sets_in_partition;
                num_sets_in_partition.push_back(num_partial_color_sets_in_partition);
                std::cout << "num_partial_color_sets_in_partition-" << partition_id << ": "
                          << num_partial_color_sets_in_partition << std::endl;
            }

            std::cout << "total num. partial color sets = " << num_partial_color_sets << std::endl;

            color_sets_builder.init_meta_color_sets_builder(
                num_integers_in_metacolor_sets + num_color_sets, num_partial_color_sets,
                p.partition_size(), num_sets_in_partition);

            std::vector<uint32_t> metacolor_set;
            metacolor_set.reserve(num_partitions);  // at most

            std::ifstream metacolor_set_in(metacolor_set_file_name(0), std::ios::binary);
            if (!metacolor_set_in.is_open()) throw std::runtime_error("error in opening file");

            uint64_t thread_id = 0;
            for (uint64_t color_set_id = 0; color_set_id != num_color_sets; ++color_set_id) {
                if (color_set_id >= thread_slices[thread_id + 1]) {
                    metacolor_set_in.close();
                    std::remove(metacolor_set_file_name(thread_id).c_str());

                    thread_id++;
                    string tmp_filename = metacolor_set_file_name(thread_id);
                    metacolor_set_in = std::ifstream(tmp_filename, std::ios::binary);
                    if (!metacolor_set_in.is_open())
                        throw std::runtime_error("error in opening file: " + tmp_filename);
                }

                assert(metacolor_set.empty());
                uint32_t meta_color_set_size = 0;
                metacolor_set_in.read(reinterpret_cast<char*>(&meta_color_set_size),
                                      sizeof(uint32_t));
                for (uint32_t i = 0; i != meta_color_set_size; ++i) {
                    uint32_t partition_id = 0;
                    uint32_t partial_color_set_id = 0;
                    metacolor_set_in.read(reinterpret_cast<char*>(&partition_id), sizeof(uint32_t));
                    metacolor_set_in.read(reinterpret_cast<char*>(&partial_color_set_id),
                                          sizeof(uint32_t));
                    /* transform the partial_color_set_id into a global id */
                    metacolor_set.push_back(partial_color_set_id +
                                            num_partial_color_sets_before[partition_id]);
                }
                color_sets_builder.encode_metacolor_set(metacolor_set.data(), metacolor_set.size());
                metacolor_set.clear();
            }

            metacolor_set_in.close();
            std::remove(metacolor_set_file_name(thread_id).c_str());
            color_sets_builder.build(idx.m_color_sets);

            timer.stop();
            std::cout << "** building partial/meta color sets took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 5. copy u2c + rank1_index and k2u");
            timer.start();
            idx.m_u2c = index.get_u2c();
            idx.m_u2c_rank1_index = index.get_u2c_rank1_index();
            idx.m_k2u = index.get_k2u();
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

            std::vector<uint32_t> permuted_set;
            permuted_set.reserve(num_colors);

            uint8_t progress_bar_size = 20;
            uint8_t curr_progress = 0;
            std::string progress_bar(progress_bar_size, ' ');
            uint64_t color_set_id = 0;
            for (; color_set_id != num_color_sets; ++color_set_id)  //
            {
                if (color_set_id >= 1.0 * curr_progress * num_color_sets / progress_bar_size) {
                    progress_bar[curr_progress++] = '#';
                }
                if (color_set_id % 1000 == 0) {
                    std::cout << "\r Progress: [" << progress_bar << "] " << color_set_id << "/"
                              << num_color_sets << std::flush;
                }
                auto it_exp = index.color_set(color_set_id);
                auto it_got = idx.color_set(color_set_id);
                const uint64_t exp_size = it_exp.size();
                const uint64_t got_size = it_got.size();

                if (exp_size != got_size) {
                    std::cout << "\033[1;31m"
                              << "got colors set of size " << got_size << " but expected "
                              << exp_size << " (color_set: " << color_set_id << ")\033[0m"
                              << std::endl;
                    return;
                }

                permuted_set.clear();
                for (uint64_t i = 0; i != exp_size; ++i, ++it_exp) {
                    uint32_t ref_id = *it_exp;
                    permuted_set.push_back(permutation[ref_id]);
                }
                std::sort(permuted_set.begin(), permuted_set.end());

                for (uint64_t i = 0; i != got_size; ++i, ++it_got) {
                    if (permuted_set[i] != *it_got) {
                        std::cout << "\033[1;31m"
                                  << "got ref " << *it_got << " BUT expected " << permuted_set[i]
                                  << "(color_set: " << color_set_id << ")"
                                  << "\033[0m" << std::endl;
                        return;
                    }
                }
            }
            std::cout << "\r Progress: [" << progress_bar << "] " << color_set_id << "/"
                      << num_color_sets << std::endl;
            essentials::logger("DONE!");
        }
    }

private:
    build_configuration m_build_config;

    std::string metacolor_set_file_name(uint32_t id) {
        return m_build_config.tmp_dirname + "/metacolor_set_" + std::to_string(id) + ".bin";
    }
};

}  // namespace fulgor

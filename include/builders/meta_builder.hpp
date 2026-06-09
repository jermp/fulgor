#pragma once

#include "include/index.hpp"
#include "include/build_util.hpp"
#include <span>

#include <shared_mutex>

namespace fulgor {

struct partition_endpoint {
    uint64_t begin, end;  // [..)
};

struct permuter {
    permuter(build_configuration const& build_config)
        : m_build_config(build_config), m_num_partitions(0), m_max_partition_size(0) {}

    void permute(hfur_index_t const& index) {
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

            m_partition_starts.resize(m_num_partitions + 1, 0);
            for (auto c : clustering_data.clusters) m_partition_starts[c] += 1;

            /* take prefix sums */
            uint64_t val = 0;
            for (auto& size : m_partition_starts) {
                if (size > m_max_partition_size) m_max_partition_size = size;
                uint64_t tmp = size;
                size = val;
                val += tmp;
            }

            const uint64_t num_colors = index.num_colors();

            /* build permutation */
            auto counts = m_partition_starts;  // copy
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

    partition_endpoint partition_endpoints(const uint64_t partition_id) const {
        assert(partition_id + 1 < m_partition_starts.size());
        return {m_partition_starts[partition_id], m_partition_starts[partition_id + 1]};
    }

    uint64_t num_partitions() const { return m_num_partitions; }
    uint64_t max_partition_size() const { return m_max_partition_size; }
    std::vector<uint32_t>& permutation() { return m_permutation; }
    std::vector<uint32_t> partition_starts() const { return m_partition_starts; }
    std::vector<std::string> filenames() const { return m_filenames; }

private:
    build_configuration m_build_config;
    uint64_t m_num_partitions;
    uint64_t m_max_partition_size;
    std::vector<uint32_t> m_permutation;
    std::vector<uint32_t> m_partition_starts;
    std::vector<std::string> m_filenames;
};

template <typename ColorSets>
struct index<ColorSets>::meta_builder {
    meta_builder() {}

    meta_builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.num_kmers() != 0) throw std::runtime_error("index already built");

        essentials::logger("step 1. loading index to be partitioned...");
        essentials::load(
            m_base_index,
            m_build_config.index_filename_to_partition.c_str());  // TODO: requires custom loader
        essentials::logger("DONE");

        const uint64_t num_colors = m_base_index.num_colors();
        const uint64_t num_color_sets = m_base_index.num_color_sets();

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        permuter p(m_build_config);
        p.permute(m_base_index);
        std::swap(m_permutation, p.permutation());

        const uint64_t num_partitions = p.num_partitions();
        const uint64_t max_partition_size = p.max_partition_size();
        if (m_build_config.verbose) {
            std::cout << "num_partitions = " << num_partitions << std::endl;
            std::cout << "max_partition_size = " << max_partition_size << std::endl;
        }

        {
            essentials::logger("step 4. building partial/meta color sets");
            timer.start();

            std::atomic<uint64_t> num_integers_in_metacolor_sets = 0;
            typename ColorSets::builder color_sets_builder(
                num_colors, p.partition_starts(), m_build_config.tmp_dirname,
                m_build_config.ram_limit_in_GiB << 30, m_build_config.verbose);

            std::string metacolor_sets_filename =
                m_build_config.tmp_dirname + "/metacolor_sets.bin";
            std::ofstream metacolor_sets_ofstream(metacolor_sets_filename,
                                                  std::ios::binary | std::ios::trunc);
            if (!metacolor_sets_ofstream.is_open()) {
                throw std::runtime_error("error in opening file");
            }
            using QueueElem = std::pair<uint64_t, std::vector<uint32_t>>;
            // TODO: capacity based on max RAM?
            util::bounded_priority_queue<QueueElem, util::compare_first> q(1 << 10);
            std::mutex write_mutex;
            std::atomic<uint64_t> num_written_meta_sets = 0;

            auto flush_meta_queue = [&num_written_meta_sets, &q, &metacolor_sets_ofstream] {
                QueueElem item;
                while (q.try_pop_if(item, [&num_written_meta_sets](const QueueElem& top_item) {
                    return top_item.first == num_written_meta_sets;
                })) {
                    const auto& meta_set = item.second;
                    const uint32_t size = meta_set.size() / 2;

                    metacolor_sets_ofstream.write(reinterpret_cast<char const*>(&size),
                                                  sizeof(uint32_t));
                    metacolor_sets_ofstream.write(reinterpret_cast<char const*>(meta_set.data()),
                                                  meta_set.size() * sizeof(uint32_t));
                    ++num_written_meta_sets;
                }
            };

            auto process_color_set = [this, num_colors, &color_sets_builder, &num_written_meta_sets,
                                      &num_integers_in_metacolor_sets, &write_mutex, &q,
                                      &flush_meta_queue](const uint64_t color_set_id) {
                std::vector<uint32_t> permuted_set;
                permuted_set.reserve(num_colors);

                const auto color_set = util::range_view(m_base_index.color_set(color_set_id));
                for (const auto& color : color_set) {
                    permuted_set.push_back(m_permutation[color]);
                }
                std::ranges::sort(permuted_set);

                std::vector<uint32_t> metacolor_set = color_sets_builder.encode(permuted_set);
                const uint32_t metacolor_set_size = metacolor_set.size() / 2;
                num_integers_in_metacolor_sets += metacolor_set_size;

                q.push(std::make_pair(color_set_id, std::move(metacolor_set)));

                const bool is_next_in_sequence = color_set_id == num_written_meta_sets;
                std::unique_lock write_lock(write_mutex, std::defer_lock);

                if (is_next_in_sequence) {
                    write_lock.lock();
                } else {
                    write_lock.try_lock();
                }
                if (write_lock.owns_lock()) {
                    flush_meta_queue();
                }
            };

            kmeans::thread_pool threads(m_build_config.num_threads);
            for (uint64_t color_set_id = 0; color_set_id < num_color_sets; ++color_set_id) {
                threads.enqueue([&, color_set_id] { process_color_set(color_set_id); });
            }
            threads.wait();

            color_sets_builder.flush();
            flush_meta_queue();
            metacolor_sets_ofstream.close();

            color_sets_builder.init_meta_color_sets_builder(num_integers_in_metacolor_sets +
                                                            num_color_sets);

            std::vector<std::pair<uint32_t, uint32_t>> metacolor_set;
            metacolor_set.reserve(num_partitions);  // at most

            std::ifstream metacolor_set_in(metacolor_sets_filename, std::ios::binary);
            if (!metacolor_set_in.is_open()) throw std::runtime_error("error in opening file");

            for (uint64_t color_set_id = 0; color_set_id != num_color_sets; ++color_set_id) {
                assert(metacolor_set.empty());
                uint32_t size = 0;
                metacolor_set_in.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
                for (uint32_t i = 0; i != size; ++i) {
                    std::array<uint32_t, 2> metacolor;
                    metacolor_set_in.read(reinterpret_cast<char*>(metacolor.data()),
                                          sizeof(metacolor));
                    metacolor_set.emplace_back(metacolor);
                }
                color_sets_builder.encode_metacolor_set(metacolor_set);
                metacolor_set.clear();
            }

            metacolor_set_in.close();
            std::remove(metacolor_sets_filename.c_str());
            color_sets_builder.build(idx.m_color_sets);

            timer.stop();
            std::cout << "** building partial/meta color sets took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 5. copy u2c + rank1_index and k2u");
            timer.start();
            idx.m_u2c = m_base_index.get_u2c();
            idx.m_u2c_rank1_index = m_base_index.get_u2c_rank1_index();
            idx.m_k2u = m_base_index.get_k2u();
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
    }

    void check(index const& idx) {
        const uint64_t num_color_sets = idx.num_color_sets();
        const uint64_t num_colors = idx.num_colors();
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        essentials::logger("checking correctness...");
        timer.start();

        std::atomic<uint64_t> num_checked_color_sets(0);

        uint64_t load = 0;
        for (uint64_t color_set_id = 0; color_set_id != num_color_sets; ++color_set_id) {
            load += idx.color_set(color_set_id).size();
        }
        const uint64_t load_per_thread = load / m_build_config.num_threads + 1;

        auto exe = [this, &idx, &num_checked_color_sets, num_colors, num_color_sets](
                       const uint64_t start, const uint64_t end) {
            assert(end > start);
            std::vector<uint32_t> permuted_set;
            permuted_set.reserve(num_colors);

            for (uint64_t color_set_id = start; color_set_id != end; ++color_set_id) {
                auto it_exp = m_base_index.color_set(color_set_id);
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
                    const uint32_t ref_id = *it_exp;
                    permuted_set.push_back(m_permutation[ref_id]);
                }
                std::ranges::sort(permuted_set);

                for (uint64_t i = 0; i != got_size; ++i, ++it_got) {
                    if (permuted_set[i] != *it_got) {
                        std::cout << "\033[1;31m"
                                  << "got ref " << *it_got << " but expected " << permuted_set[i]
                                  << "(color_set: " << color_set_id << ")"
                                  << "\033[0m" << std::endl;
                        return;
                    }
                }

                if (++num_checked_color_sets % 1000 == 0) {
                    std::cout << "\rChecked " << num_checked_color_sets << "/" << num_color_sets
                              << " color sets" << std::flush;
                }
            }
        };

        std::vector<std::thread> threads(m_build_config.num_threads);
        uint64_t start = 0, curr_set_id = 0, curr_load = 0;
        for (uint64_t thread_id = 0; thread_id != m_build_config.num_threads; ++thread_id) {
            while (curr_load < load_per_thread && curr_set_id < num_color_sets) {
                curr_load += idx.color_set(curr_set_id).size();
                ++curr_set_id;
            }
            threads[thread_id] = std::thread(exe, start, curr_set_id);
            start = curr_set_id;
            curr_load = 0;
        }
        for (auto& t : threads) {
            if (t.joinable()) t.join();
        }

        std::cout << "\rChecked " << num_checked_color_sets << "/" << num_color_sets
                  << " color sets" << std::endl;

        timer.stop();
        std::cout << "** checking correctness took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        essentials::logger("DONE!");
    }

private:
    build_configuration m_build_config;
    hfur_index_t m_base_index;
    std::vector<uint32_t> m_permutation;

    std::string metacolor_set_file_name(const uint32_t id) const {
        return m_build_config.tmp_dirname + "/metacolor_set_" + std::to_string(id) + ".bin";
    }
};

}  // namespace fulgor

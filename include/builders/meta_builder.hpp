#pragma once

#include "stream.hpp"
#include "include/index.hpp"
#include "include/build_util.hpp"
#include <span>

#include <shared_mutex>
#include <utility>

namespace fulgor {

struct partition_endpoint {
    uint64_t begin, end;  // [..)
};

struct permuter {
    explicit permuter(build_configuration build_config)
        : m_build_config(std::move(build_config)), m_num_partitions(0), m_max_partition_size(0) {}

    void compute_permutation() {
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            essentials::logger("step 2.1. build sketches");
            timer.start();
            constexpr uint64_t p = 10;  // use 2^p bytes per HLL sketch
            build_reference_sketches(
                m_build_config.num_colors, p, m_build_config.num_threads,
                m_build_config.tmp_filename(m_build_config.output_filename.filename()),
                m_build_config.tmp_filename("sketches.bin"));
            timer.stop();
            std::cout << "** building sketches took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 2.2. clustering sketches");
            timer.start();

            std::ifstream in(m_build_config.tmp_filename("sketches.bin"), std::ios::binary);
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

            std::remove(m_build_config.tmp_filename("sketches.bin").c_str());

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

            /* build permutation */
            auto counts = m_partition_starts;  // copy
            m_permutation.resize(m_build_config.num_colors);
            assert(clustering_data.clusters.size() == m_build_config.num_colors);
            for (uint64_t i = 0; i != m_build_config.num_colors; ++i) {
                uint32_t cluster_id = clustering_data.clusters[i];
                m_permutation[i] = counts[cluster_id];
                counts[cluster_id] += 1;
            }
        }
    }

    partition_endpoint partition_endpoints(const uint64_t partition_id) const {
        assert(partition_id + 1 < m_partition_starts.size());
        return {m_partition_starts[partition_id], m_partition_starts[partition_id + 1]};
    }

    uint64_t num_partitions() const { return m_num_partitions; }
    uint64_t max_partition_size() const { return m_max_partition_size; }
    std::vector<uint32_t> permutation() { return m_permutation; }
    std::vector<uint32_t> partition_starts() const { return m_partition_starts; }

    filenames permuted_filenames() const {
        std::ifstream file(m_build_config.filenames_list, std::ios::binary);

        std::vector<std::string> permuted_filenames;
        permuted_filenames.resize(m_build_config.num_colors);
        for (uint64_t i = 0; i != m_build_config.num_colors; ++i) {
            std::getline(file, permuted_filenames[m_permutation[i]]);
        }

        filenames fn;
        fn.build(permuted_filenames);
        return fn;
    }

    void apply(std::vector<uint32_t>& color_set) const {
        for (auto& color : color_set) {
            color = m_permutation[color];
        }
        std::ranges::sort(color_set);
    }

private:
    build_configuration m_build_config;
    uint64_t m_num_partitions;
    uint64_t m_max_partition_size;
    std::vector<uint32_t> m_permutation;
    std::vector<uint32_t> m_partition_starts;
};

template <typename ColorSets>
struct index<ColorSets>::meta_builder {
    meta_builder(build_configuration const& build_config)
        : m_build_config(build_config), m_saver(build_config.output_filename) {}

    void build(index& idx) {
        if (idx.m_k2u.num_kmers() != 0) throw std::runtime_error("index already built");
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        cdbg::build_config cdbg_build_config;
        cdbg_build_config.filenames_list = m_build_config.filenames_list;
        cdbg_build_config.out_basename =
            m_build_config.tmp_filename(m_build_config.output_filename.filename());
        cdbg_build_config.k = m_build_config.k;
        cdbg_build_config.m = m_build_config.m;
        cdbg_build_config.num_threads = m_build_config.num_threads;
        cdbg_build_config.max_ram_gb = m_build_config.ram_limit_in_GiB;
        uint64_t curr_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                               std::chrono::system_clock::now().time_since_epoch())
                               .count();
        cdbg_build_config.tmp_dir =
            std::format("{}/cdbg_build_{}", m_build_config.tmp_dirname.string(), curr_ms);
        cdbg::builder cdbg_builder(cdbg_build_config);

        {
            essentials::logger("step 1. build colored compacted dBG...");
            timer.start();

            cdbg_builder.build();
            m_build_config.num_colors = cdbg_builder.num_colors();

            timer.stop();
            std::cout << "** building the ccdBG took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        const uint64_t num_colors = cdbg_builder.num_colors();
        const uint64_t num_color_sets = cdbg_builder.num_color_sets();

        permuter p(m_build_config);
        {
            essentials::logger("step 2. compute partition and permutation...");
            timer.start();

            p.compute_permutation();
            timer.stop();
            std::cout << "** computing the partition took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        const uint64_t num_partitions = p.num_partitions();
        const uint64_t max_partition_size = p.max_partition_size();
        if (m_build_config.verbose) {
            std::cout << "num_partitions = " << num_partitions << std::endl;
            std::cout << "max_partition_size = " << max_partition_size << std::endl;
        }

        {
            essentials::logger("step 3. build partial/meta color sets");
            timer.start();

            std::atomic<uint64_t> num_integers_in_metacolor_sets = 0;
            typename ColorSets::builder color_sets_builder(
                num_colors, m_saver, p.partition_starts(), m_build_config.tmp_dirname,
                m_build_config.ram_limit_in_GiB << 30, m_build_config.verbose);

            std::string metacolor_sets_filename = m_build_config.tmp_filename("metacolor_sets.bin");
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

            const uint64_t max_queue_size = m_build_config.num_threads * 2;
            cdbg::unitigs_color_set_stream stream(cdbg_build_config.out_basename, max_queue_size);
            stream.start();

            auto process_color_set = [this, &color_sets_builder, &num_written_meta_sets,
                                      &num_integers_in_metacolor_sets, &write_mutex, &q,
                                      &flush_meta_queue, &stream, &p] {
                for (auto opt = stream.get(); opt != std::nullopt; opt = stream.get()) {
                    auto& [unitig_start, num_unitigs, cs_id, color_set] = opt.value();
                    p.apply(color_set);

                    std::vector<uint32_t> metacolor_set = color_sets_builder.encode(color_set);
                    const uint32_t metacolor_set_size = metacolor_set.size() / 2;
                    num_integers_in_metacolor_sets += metacolor_set_size;

                    q.push(std::make_pair(cs_id, std::move(metacolor_set)));

                    const bool is_next_in_sequence = cs_id == num_written_meta_sets;
                    std::unique_lock write_lock(write_mutex, std::defer_lock);

                    if (is_next_in_sequence) {
                        write_lock.lock();
                    } else {
                        write_lock.try_lock();
                    }
                    if (write_lock.owns_lock()) {
                        flush_meta_queue();
                    }
                }
            };

            std::vector<std::thread> threads(m_build_config.num_threads - 1);
            for (uint64_t thread_id = 0; thread_id < m_build_config.num_threads - 1; ++thread_id) {
                threads[thread_id] = std::thread(process_color_set);
            }
            for (auto& thread : threads) {
                if (thread.joinable()) thread.join();
            }

            color_sets_builder.flush();
            flush_meta_queue();
            metacolor_sets_ofstream.close();

            color_sets_builder.init_meta_color_sets_builder(num_integers_in_metacolor_sets +
                                                            num_color_sets);

            std::vector<std::pair<uint32_t, uint32_t>> metacolor_set;
            metacolor_set.reserve(num_partitions);  // at most

            std::ifstream metacolor_set_in(metacolor_sets_filename, std::ios::binary);
            if (!metacolor_set_in.is_open()) throw std::runtime_error("error in opening file");

            for (uint64_t color_set_id = 0; color_set_id < num_color_sets; ++color_set_id) {
                assert(metacolor_set.empty());
                uint32_t size = 0;
                metacolor_set_in.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
                for (uint32_t i = 0; i != size; ++i) {
                    uint32_t metacolor[2];
                    metacolor_set_in.read(reinterpret_cast<char*>(metacolor), sizeof(metacolor));
                    metacolor_set.emplace_back(metacolor[0], metacolor[1]);
                }
                color_sets_builder.encode_metacolor_set(metacolor_set);
                metacolor_set.clear();
            }

            metacolor_set_in.close();
            std::remove(metacolor_sets_filename.c_str());
            color_sets_builder.build();

            timer.stop();
            std::cout << "** building partial/meta color sets took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 4. copy u2c and build rank1_index");
            timer.start();

            bits::bit_vector u2c;
            bits::rank9 u2c_rank1_index;
            essentials::load(u2c, cdbg_build_config.u2c_filename().c_str());
            u2c_rank1_index.build(u2c);
            m_saver.visit(u2c);
            m_saver.visit(u2c_rank1_index);

            std::cout << "m_u2c.num_bits() " << u2c.num_bits() << std::endl;
            std::cout << "m_u2c_rank1_index.num_ones() " << u2c_rank1_index.num_ones() << std::endl;

            timer.stop();
            std::cout << "** copying u2c and building rank1_index took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 5. building SSHash...");
            timer.start();

            sshash::build_configuration sshash_config;
            sshash_config.k = m_build_config.k;
            sshash_config.m = m_build_config.m;
            sshash_config.canonical = true;
            sshash_config.verbose = m_build_config.verbose;
            sshash_config.tmp_dirname = m_build_config.tmp_dirname;
            sshash_config.num_threads = m_build_config.num_threads;
            sshash_config.print();

            sshash::dictionary_type k2u;
            k2u.build(cdbg_build_config.fa_filename(), sshash_config);
            m_saver.visit(k2u);

            try {  // remove unitig file
                std::remove(cdbg_build_config.fa_filename().c_str());
            } catch (std::exception const& e) {
                std::cerr << e.what() << std::endl;
            }

            timer.stop();
            std::cout << "** building SSHash took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 6. building filenames");
            timer.start();

            m_saver.visit(p.permuted_filenames());

            timer.stop();
            std::cout << "** building filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            m_saver.write(constants::current_version_number::major);
            m_saver.write(constants::current_version_number::minor);
            m_saver.write(constants::current_version_number::patch);
        }
    }

    void check(index const& idx) {
        // const uint64_t num_color_sets = idx.num_color_sets();
        // const uint64_t num_colors = idx.num_colors();
        // essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        // essentials::logger("checking correctness...");
        // timer.start();
        //
        // std::atomic<uint64_t> num_checked_color_sets(0);
        //
        // uint64_t load = 0;
        // for (uint64_t color_set_id = 0; color_set_id != num_color_sets; ++color_set_id) {
        //     load += idx.color_set(color_set_id).size();
        // }
        // const uint64_t load_per_thread = load / m_build_config.num_threads + 1;
        //
        // auto exe = [this, &idx, &num_checked_color_sets, num_colors, num_color_sets](
        //                const uint64_t start, const uint64_t end) {
        //     assert(end > start);
        //     std::vector<uint32_t> permuted_set;
        //     permuted_set.reserve(num_colors);
        //
        //     for (uint64_t color_set_id = start; color_set_id != end; ++color_set_id) {
        //         auto it_exp = m_base_index.color_set(color_set_id);
        //         auto it_got = idx.color_set(color_set_id);
        //         const uint64_t exp_size = it_exp.size();
        //         const uint64_t got_size = it_got.size();
        //
        //         if (exp_size != got_size) {
        //             std::cout << "\033[1;31m"
        //                       << "got colors set of size " << got_size << " but expected "
        //                       << exp_size << " (color_set: " << color_set_id << ")\033[0m"
        //                       << std::endl;
        //             return;
        //         }
        //
        //         permuted_set.clear();
        //         for (uint64_t i = 0; i != exp_size; ++i, ++it_exp) {
        //             const uint32_t ref_id = *it_exp;
        //             permuted_set.push_back(m_permutation[ref_id]);
        //         }
        //         std::ranges::sort(permuted_set);
        //
        //         for (uint64_t i = 0; i != got_size; ++i, ++it_got) {
        //             if (permuted_set[i] != *it_got) {
        //                 std::cout << "\033[1;31m"
        //                           << "got ref " << *it_got << " but expected " << permuted_set[i]
        //                           << "(color_set: " << color_set_id << ")"
        //                           << "\033[0m" << std::endl;
        //                 return;
        //             }
        //         }
        //
        //         if (++num_checked_color_sets % 1000 == 0) {
        //             std::cout << "\rChecked " << num_checked_color_sets << "/" << num_color_sets
        //                       << " color sets" << std::flush;
        //         }
        //     }
        // };
        //
        // std::vector<std::thread> threads(m_build_config.num_threads);
        // uint64_t start = 0, curr_set_id = 0, curr_load = 0;
        // for (uint64_t thread_id = 0; thread_id != m_build_config.num_threads; ++thread_id) {
        //     while (curr_load < load_per_thread && curr_set_id < num_color_sets) {
        //         curr_load += idx.color_set(curr_set_id).size();
        //         ++curr_set_id;
        //     }
        //     threads[thread_id] = std::thread(exe, start, curr_set_id);
        //     start = curr_set_id;
        //     curr_load = 0;
        // }
        // for (auto& t : threads) {
        //     if (t.joinable()) t.join();
        // }
        //
        // std::cout << "\rChecked " << num_checked_color_sets << "/" << num_color_sets
        //           << " color sets" << std::endl;
        //
        // timer.stop();
        // std::cout << "** checking correctness took " << timer.elapsed() << " seconds / "
        //           << timer.elapsed() / 60 << " minutes" << std::endl;
        // essentials::logger("DONE!");
    }

private:
    build_configuration m_build_config;
    util::external_saver m_saver;

    std::string metacolor_set_file_name(const uint32_t id) const {
        return m_build_config.tmp_filename(std::format("metacolor_set_{}.bin", id));
    }
};

}  // namespace fulgor

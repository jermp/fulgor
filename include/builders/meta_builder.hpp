#pragma once

#include "stream.hpp"
#include "external/kmeans/include/kmeans.hpp"
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
    explicit permuter(build_configuration const& build_config)
        : m_build_config(build_config), m_num_partitions(0), m_max_partition_size(0) {}

    void compute_permutation(uint32_t num_colors) {
        m_num_colors = num_colors;
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            util::timed_phase timer_(" step 2.1. build sketches");
            constexpr uint64_t p = 10;  // use 2^p bytes per HLL sketch
            build_reference_sketches(
                num_colors, p, m_build_config.num_threads,
                m_build_config.tmp_filename(m_build_config.output_filename.filename()),
                m_build_config.tmp_filename("sketches.bin"));
        }

        {
            util::timed_phase timer_(" step 2.2. cluster sketches");
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
            m_permutation.resize(num_colors);
            assert(clustering_data.clusters.size() == num_colors);
            for (uint64_t i = 0; i != num_colors; ++i) {
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
        permuted_filenames.resize(m_num_colors);
        for (uint64_t i = 0; i != m_num_colors; ++i) {
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
    uint32_t m_num_colors;
    uint64_t m_num_partitions;
    uint64_t m_max_partition_size;
    std::vector<uint32_t> m_permutation;
    std::vector<uint32_t> m_partition_starts;
};

template <typename ColorSets>
struct meta_build_strategy {
    explicit meta_build_strategy(build_configuration const& build_config)
        : m_build_config(build_config)
        , m_saver(build_config.output_filename.string())
        , permuter_(build_config) {}

    void build_cdbg() {
        util::timed_phase timer("step 1. build colored compacted dBG.");
        cdbg_config.filenames_list = m_build_config.filenames_list;
        cdbg_config.out_basename =
            m_build_config.tmp_filename(m_build_config.output_filename.filename());
        cdbg_config.k = m_build_config.k;
        cdbg_config.m = m_build_config.m;
        cdbg_config.num_threads = m_build_config.num_threads;
        cdbg_config.max_ram_gb = m_build_config.ram_limit_in_GiB;
        uint64_t curr_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                               std::chrono::system_clock::now().time_since_epoch())
                               .count();
        cdbg_config.tmp_dir = m_build_config.tmp_dirname / std::format("cdbg_build_{}", curr_ms);
        cdbg::builder cdbg_builder(cdbg_config);

        cdbg_builder.build();
        m_build_config.num_colors = cdbg_builder.num_colors();
        m_num_color_sets = cdbg_builder.num_color_sets();
    }

    void build_color_sets() {
        util::timed_phase timer("step 2. compute partition and permutation");
        permuter_.compute_permutation(m_build_config.num_colors);

        const uint32_t num_colors = m_build_config.num_colors;
        const uint64_t num_partitions = permuter_.num_partitions();
        const uint64_t max_partition_size = permuter_.max_partition_size();
        if (m_build_config.verbose) {
            std::cout << "num_partitions = " << num_partitions << std::endl;
            std::cout << "max_partition_size = " << max_partition_size << std::endl;
        }

        {
            util::timed_phase timer(" step 2.3. build partial/meta color sets");

            std::atomic<uint64_t> num_integers_in_metacolor_sets = 0;
            typename ColorSets::builder color_sets_builder(
                num_colors, m_saver, permuter_.partition_starts(), m_build_config.tmp_dirname,
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
            std::atomic<uint64_t> num_processed_sets = 0;

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
            cdbg::unitigs_color_set_stream stream(cdbg_config.out_basename, max_queue_size);
            stream.start();

            auto process_color_set = [this, &color_sets_builder, &num_written_meta_sets,
                                      &num_integers_in_metacolor_sets, &write_mutex, &q,
                                      &flush_meta_queue, &stream, &num_processed_sets] {
                for (auto opt = stream.get(); opt != std::nullopt; opt = stream.get()) {
                    auto& [unitig_start, num_unitigs, cs_id, color_set] = opt.value();
                    permuter_.apply(color_set);

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

                    if (num_processed_sets.fetch_add(1) % 10000 == 0) {
                        auto progress = std::format("\r[encode-sets] {}/{} ({:.2f}%)",
                                                    num_processed_sets.load(), m_num_color_sets,
                                                    100. * num_processed_sets / m_num_color_sets);
                        std::cout << progress << std::flush;
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

            auto progress =
                std::format("\r[encode-sets] {}/{} ({:.2f}%)", num_processed_sets.load(),
                            m_num_color_sets, 100. * num_processed_sets / m_num_color_sets);
            std::cout << progress << std::endl;

            color_sets_builder.flush();
            flush_meta_queue();
            metacolor_sets_ofstream.close();

            color_sets_builder.init_meta_color_sets_builder(num_integers_in_metacolor_sets +
                                                            m_num_color_sets);

            std::vector<std::pair<uint32_t, uint32_t>> metacolor_set;
            metacolor_set.reserve(num_partitions);  // at most

            std::ifstream metacolor_set_in(metacolor_sets_filename, std::ios::binary);
            if (!metacolor_set_in.is_open()) throw std::runtime_error("error in opening file");

            for (uint64_t color_set_id = 0; color_set_id < m_num_color_sets; ++color_set_id) {
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

            try {
                std::remove(cdbg_config.cs_filename().c_str());
            } catch (std::exception const& e) {
                std::cerr << e.what() << std::endl;
            }
        }
    }

    void build_u2c() {
        util::timed_phase timer("step 3. copy u2c and build rank1_index");

        bits::bit_vector u2c;
        bits::rank9 u2c_rank1_index;
        essentials::load(u2c, cdbg_config.u2c_filename().c_str());
        u2c_rank1_index.build(u2c);
        m_saver.visit(u2c);
        m_saver.visit(u2c_rank1_index);

        std::cout << "m_u2c.num_bits() = " << u2c.num_bits() << std::endl;
        std::cout << "m_u2c_rank1_index.num_ones() = " << u2c_rank1_index.num_ones() << std::endl;

        try {  // remove u2c file
            std::remove(cdbg_config.u2c_filename().c_str());
        } catch (std::exception const& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    void build_kmer_dictionary() {
        util::timed_phase timer("step 4. build SSHash");

        sshash::build_configuration sshash_config;
        sshash_config.k = m_build_config.k;
        sshash_config.m = m_build_config.m;
        sshash_config.canonical = true;
        sshash_config.verbose = m_build_config.verbose;
        sshash_config.tmp_dirname = m_build_config.tmp_dirname;
        sshash_config.num_threads = m_build_config.num_threads;
        sshash_config.print();

        sshash::dictionary_type k2u;
        k2u.build(cdbg_config.fa_filename(), sshash_config);
        m_saver.visit(k2u);

        try {  // remove unitig file
            std::remove(cdbg_config.fa_filename().c_str());
        } catch (std::exception const& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    void build_filenames() {
        util::timed_phase timer("step 5. permute and write filenames");
        m_saver.visit(permuter_.permuted_filenames());
    }

    ~meta_build_strategy() {
        std::remove(cdbg_config.cs_filename().c_str());
        std::remove(cdbg_config.fa_filename().c_str());
        std::remove(cdbg_config.u2c_filename().c_str());
        std::remove(cdbg_config.metadata_filename().c_str());
    }

    util::external_saver& saver() { return m_saver; }

private:
    build_configuration m_build_config;
    cdbg::build_config cdbg_config;
    util::external_saver m_saver;
    permuter permuter_;

    uint64_t m_num_color_sets{};

    std::string metacolor_set_file_name(const uint32_t id) const {
        return m_build_config.tmp_filename(std::format("metacolor_set_{}.bin", id));
    }
};

}  // namespace fulgor

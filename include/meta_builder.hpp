#pragma once

#include "index.hpp"

namespace fulgor {

template <typename ColorClasses>
struct index<ColorClasses>::meta_builder {
    meta_builder() {}

    struct partition_endpoint {
        uint64_t begin, end;
    };

    meta_builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        uint64_t num_partitions = 0;
        uint64_t num_integers_in_metacolors = 0;
        uint64_t num_lists = 0;

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        essentials::logger("step 1. loading index to be partitioned...");
        essentials::load(m_index, m_build_config.index_filename_to_partition.c_str());
        essentials::logger("DONE");

        {
            essentials::logger("step 2.1. build colors");
            timer.start();

            std::vector<uint32_t> buffer;  // buffer the list of the current partition
            typename ColorClasses::builder colors_builder;

            {
                std::ifstream in(m_build_config.partitions_filename.c_str());
                uint64_t n = 0;
                while (in >> n) m_partition_size.push_back(n);
                assert(m_partition_size.size() > 1);
                assert(m_partition_size.front() == 0);
                num_partitions = m_partition_size.size() - 1;
                in.close();

                std::cout << "num_partitions = " << num_partitions << std::endl;

                uint64_t max_partition_size =
                    *std::max_element(m_partition_size.begin(), m_partition_size.end());
                buffer.reserve(max_partition_size);

                m_hashes.resize(num_partitions);
            }

            colors_builder.init_colors_builder(m_index.num_docs(), num_partitions);
            for (uint64_t partition_id = 0; partition_id != num_partitions; ++partition_id) {
                auto endpoints = partition_endpoints(partition_id);
                uint64_t num_docs_in_partition = endpoints.end - endpoints.begin;
                colors_builder.init_color_partition(partition_id, num_docs_in_partition);
            }

            uint64_t partition_id = 0;

            auto hash_and_compress = [&]() {
                auto hash = util::hash128(reinterpret_cast<char const*>(buffer.data()),
                                          buffer.size() * sizeof(uint32_t));
                if (auto it = m_hashes[partition_id].find(hash);
                    it == m_hashes[partition_id].cend()) {
                    uint32_t id = m_hashes[partition_id].size();
                    m_hashes[partition_id].insert({hash, id});
                    colors_builder.process_colors(partition_id, buffer.data(), buffer.size());
                }
                buffer.clear();
                num_integers_in_metacolors += 1;
            };

            uint64_t num_color_classes = m_index.num_color_classes();
            for (uint64_t color_class_id = 0; color_class_id != num_color_classes;
                 ++color_class_id) {
                auto it = m_index.colors(color_class_id);
                uint64_t list_size = it.size();
                partition_id = 0;
                partition_endpoint curr_partition = partition_endpoints(0);
                for (uint64_t i = 0; i != list_size; ++i, ++it) {
                    uint32_t ref_id = *it;
                    while (ref_id >= curr_partition.end) {
                        if (!buffer.empty()) hash_and_compress();
                        partition_id += 1;
                        curr_partition = partition_endpoints(partition_id);
                    }
                    assert(ref_id >= curr_partition.begin);
                    buffer.push_back(ref_id - curr_partition.begin);
                }
                if (!buffer.empty()) hash_and_compress();
            }

            uint64_t prev_id = 0;
            m_num_lists_in_partition.reserve(num_partitions);
            for (partition_id = 0; partition_id != num_partitions; ++partition_id) {
                uint64_t num_lists_in_partition = m_hashes[partition_id].size();
                num_lists += num_lists_in_partition;
                m_num_lists_in_partition.push_back(num_lists_in_partition);
                std::cout << "num_lists in partition-" << partition_id << ": "
                          << num_lists_in_partition << std::endl;
                for (auto& p : m_hashes[partition_id]) { p.second += prev_id; }
                prev_id += num_lists_in_partition;
            }

            std::cout << "total num_colors = " << num_lists << std::endl;

            essentials::logger("step 2.2. build metacolors");

            colors_builder.init_meta_colors_builder(
                num_integers_in_metacolors + m_index.num_color_classes(), num_lists,
                m_partition_size, m_num_lists_in_partition);

            std::vector<uint32_t> metacolors;
            metacolors.reserve(num_partitions);  // at most
            for (uint64_t color_class_id = 0; color_class_id != m_index.num_color_classes();
                 ++color_class_id) {
                get_metacolors(color_class_id, metacolors);
                colors_builder.process_metacolors(metacolors.data(), metacolors.size());
                metacolors.clear();
            }

            colors_builder.build(idx.m_ccs);

            timer.stop();
            std::cout << "** building and compressing colors/meta-colors took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 3. copy m_u2c, m_k2u, m_filenames");
            timer.start();
            idx.m_u2c = m_index.get_u2c();
            idx.m_k2u = m_index.get_dict();
            idx.m_filenames = m_index.get_filenames();
            timer.stop();
            std::cout << "** copying m_u2c, m_k2u, m_filenames took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        if (m_build_config.check) {
            essentials::logger("step 4. check correctness...");
            for (uint64_t color_class_id = 0; color_class_id != m_index.num_color_classes();
                 ++color_class_id) {
                auto fwd_it_exp = m_index.colors(color_class_id);
                auto fwd_it_got = idx.colors(color_class_id);
                uint64_t exp_size = fwd_it_exp.size();
                uint64_t got_size = fwd_it_exp.size();
                if (exp_size != got_size) {
                    std::cout << "got colors list of size " << got_size << " but expected "
                              << exp_size << std::endl;
                    return;
                }
                for (uint64_t i = 0; i != exp_size; ++i, ++fwd_it_exp, ++fwd_it_got) {
                    if (*fwd_it_exp != *fwd_it_got) {
                        std::cout << "got ref " << *fwd_it_got << " BUT expected " << *fwd_it_exp
                                  << std::endl;
                        return;
                    }
                }
            }
            essentials::logger("DONE!");
        }
    }

private:
    build_configuration m_build_config;
    index_type m_index;
    std::vector<uint32_t> m_partition_size;
    std::vector<uint32_t> m_num_lists_in_partition;
    std::vector<std::unordered_map<__uint128_t, uint32_t>> m_hashes;  // (hash,id)

    void get_metacolors(uint64_t color_class_id, std::vector<uint32_t>& out) const {
        assert(color_class_id < m_index.num_color_classes());
        assert(out.empty());

        static std::vector<uint32_t> buffer;
        assert(buffer.empty());

        uint64_t partition_id = 0;

        auto hash_and_write_metacolor = [&]() {
            auto hash = util::hash128(reinterpret_cast<char const*>(buffer.data()),
                                      buffer.size() * sizeof(uint32_t));
            auto it = m_hashes[partition_id].find(hash);
            assert(it != m_hashes[partition_id].cend());  // must be found
            uint32_t metacolor = (*it).second;
            out.push_back(metacolor);
            buffer.clear();
        };

        auto it = m_index.colors(color_class_id);
        uint64_t list_size = it.size();

        partition_endpoint curr_partition = partition_endpoints(0);
        for (uint64_t i = 0; i != list_size; ++i, ++it) {
            uint32_t ref_id = *it;
            while (ref_id >= curr_partition.end) {
                if (!buffer.empty()) hash_and_write_metacolor();
                partition_id += 1;
                curr_partition = partition_endpoints(partition_id);
            }
            assert(ref_id >= curr_partition.begin);
            buffer.push_back(ref_id - curr_partition.begin);
        }
        if (!buffer.empty()) hash_and_write_metacolor();
    }

    partition_endpoint partition_endpoints(uint64_t partition_id) const {
        assert(partition_id + 1 < m_partition_size.size());
        return {m_partition_size[partition_id], m_partition_size[partition_id + 1]};
    }
};
}  // namespace fulgor
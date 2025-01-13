#pragma once

#include "external/s_indexes/include/s_index.hpp"
#include "external/s_indexes/include/builder.hpp"
#include "external/s_indexes/include/util.hpp"

namespace fulgor{

template <typename ColorSets>
struct index<ColorSets>::sliced_builder{
    sliced_builder() {}

    sliced_builder(build_configuration const& build_config): m_build_config(build_config) {}

    void build(index &idx){
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index_already_built");

        index_type index;
        essentials::logger("step 1. loading index to be sliced...");
        essentials::load(index, m_build_config.index_filename_to_partition.c_str());
        essentials::logger("DONE");

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            essentials::logger("step 2. building sliced colors");
            std::vector<uint32_t> block;
            block.reserve(::sliced::constants::block_size);
            uint64_t universe_size = index.num_colors();
            typename ColorSets::builder colors_builder(universe_size);
            uint64_t num_threads = m_build_config.num_threads;

            m_offsets.reserve(index.num_color_sets()+ 2);
            m_offsets.push_back(universe_size);
            m_offsets.push_back(0);
            
            std::vector<std::vector<uint8_t>> thread_sequences(num_threads);
            std::vector<std::vector<uint64_t>> thread_offsets(num_threads);
            std::vector<uint64_t> thread_slices(num_threads, 0);

            for(uint64_t i = 0; i < num_threads; ++i) {
                thread_slices[i] = index.num_color_sets() / num_threads * i;
            }
            thread_slices.push_back(index.num_color_sets());

            auto exe = [&](uint64_t thread_id) {
                uint64_t start = thread_slices[thread_id], end = thread_slices[thread_id + 1];
                auto& offsets = thread_offsets[thread_id];
                auto& sequences = thread_sequences[thread_id];
                for(uint64_t color_set_id = start; color_set_id < end; ++color_set_id){
                    auto it = index.color_set(color_set_id);
                    uint64_t n = it.size();
                    vector<uint32_t> color_set(n);
                    for(uint64_t i = 0; i < n; ++i, ++it){
                        color_set[i] = *it;
                    }

                    ::sliced::encode_sequence(color_set.data(), n, sequences);
                    offsets.push_back(sequences.size());
                }
                cout << "~~ Thread " << thread_id << " finished!" << endl;

            };

            std::vector<std::thread> threads(num_threads);
            for(uint64_t thread_id = 0; thread_id < num_threads; thread_id++){
                threads[thread_id] = std::thread(exe, thread_id);
            }

            for(auto& t : threads){
                if (t.joinable()) t.join();
            }

            uint64_t prev_offset = 0;
            for(uint64_t thread_id = 0; thread_id < num_threads; thread_id++){
                m_sequences.reserve(m_sequences.size() + thread_sequences[thread_id].size());
                m_sequences.insert(m_sequences.end(), thread_sequences[thread_id].begin(), thread_sequences[thread_id].end());
                for(uint64_t& offset : thread_offsets[thread_id]){
                    m_offsets.push_back(prev_offset + offset);
                }
                prev_offset += thread_offsets[thread_id].back();
            }


            cout << m_offsets.size() << " " << m_sequences.size() << " | " << m_offsets.back() << endl;

            m_offsets.pop_back();
            colors_builder.push_offsets(m_offsets);
            colors_builder.push_sequences(m_sequences);
            colors_builder.build(idx.m_color_sets);
        }  

        {
            essentials::logger("step 3. copy u2c and k2u");
            timer.start();
            idx.m_u2c = index.get_u2c();
            idx.m_k2u = index.get_k2u();
            timer.stop();
            std::cout << "** copying u2c and k2u took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 4. building filenames");
            timer.start();
            idx.m_filenames = index.get_filenames();
            timer.stop();
            std::cout << "** building filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }
    }

private:
    build_configuration m_build_config;
    std::vector<uint64_t> m_offsets;
    std::vector<uint8_t> m_sequences;
};

}

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
            block.reserve(::sliced::constants::block_size); //TODO: constants::block_size
            uint64_t universe_size = index.num_colors();
            typename ColorSets::builder colors_builder(universe_size);

            m_offsets.push_back(universe_size);
            m_offsets.push_back(0);

            for(uint64_t color_set_id = 0; color_set_id < index.num_color_sets(); ++color_set_id){
                auto it = index.color_set(color_set_id);
                uint64_t n = it.size();
                vector<uint32_t> color_set(n);
                for(uint64_t i = 0; i < n; ++i, ++it){
                    color_set[i] = *it;
                }
                // TODO: see if pass is skippable

                ::sliced::encode_sequence(color_set.data(), n, m_sequences); //TODO: implment
                m_offsets.push_back(m_sequences.size());
            }

            m_offsets.pop_back();
            colors_builder.push_offsets(m_offsets);
            colors_builder.push_sequences(m_sequences);

            colors_builder.build(idx.m_color_sets);
        }  

    }

private:
    build_configuration m_build_config;
    std::vector<uint64_t> m_offsets;
    std::vector<uint8_t> m_sequences;
};

}

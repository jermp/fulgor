#pragma once

#include "external/kmeans/include/kmeans.hpp"
#include "include/index.hpp"
#include "external/cdbg-builder/include/util.hpp"
#include "external/cdbg-builder/include/builder.hpp"
#include <span>

namespace fulgor {

template <typename ColorSets>
struct hybrid_build_strategy {
    explicit hybrid_build_strategy(build_configuration const& build_config)
        : m_build_config(build_config) {}

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
    }

    void build_color_sets() {
        util::timed_phase timer("step 2. copy color sets");
        std::filesystem::rename(cdbg_config.cs_filename(), m_build_config.output_filename);
        m_saver.set_stream(
            std::ofstream(m_build_config.output_filename, std::ios::binary | std::ios::app));
    }

    void build_u2c() {
        util::timed_phase timer("step 3. copy u2c and build rank1_index");
        bits::bit_vector u2c;
        bits::rank9 u2c_rank1_index;
        essentials::load(u2c, cdbg_config.u2c_filename().c_str());
        u2c_rank1_index.build(u2c);
        m_saver.visit(u2c);
        m_saver.visit(u2c_rank1_index);

        std::cout << "m_u2c.num_bits() " << u2c.num_bits() << std::endl;
        std::cout << "m_u2c_rank1_index.num_ones() " << u2c_rank1_index.num_ones() << std::endl;

        try {  // remove unitig file
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

        filenames filenames;
        filenames.build_from_file(m_build_config.filenames_list);
        m_saver.visit(filenames);
    }

    ~hybrid_build_strategy() {
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

    uint64_t num_color_sets{};

    std::string metacolor_set_file_name(const uint32_t id) const {
        return m_build_config.tmp_filename(std::format("metacolor_set_{}.bin", id));
    }
};
}  // namespace fulgor

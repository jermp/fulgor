#pragma once

#include "external/kmeans/include/kmeans.hpp"
#include "include/index.hpp"
#include "external/cdbg-builder/include/util.hpp"
#include "external/cdbg-builder/include/builder.hpp"

#include <span>

namespace fulgor {

template <typename ColorSets>
struct index<ColorSets>::hybrid_builder {
    hybrid_builder(build_configuration const& build_config)
        : m_build_config(build_config), m_saver(build_config.output_filename) {}

    void build(index& idx) {
        if (idx.m_k2u.num_kmers() != 0) throw std::runtime_error("index already built");

        cdbg::build_config cdbg_build_config;
        cdbg_build_config.filenames_list = m_build_config.filenames_list;
        cdbg_build_config.out_basename =
            std::format("{}/{}", m_build_config.tmp_dirname.string(),
                        std::filesystem::path(m_build_config.output_filename).filename().string());
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
            util::timed_phase("step 1. build colored compacted dBG");

            cdbg_builder.build();
            m_build_config.num_colors = cdbg_builder.num_colors();
        }

        {
            util::timed_phase("step 2. copy color sets");

            std::ifstream cs_file(cdbg_build_config.cs_filename(), std::ios::binary);
            m_saver.append(cs_file, 12);      // write num_cols, sp_thresh, vd_thresh
            cs_file.seekg(8, std::ios::cur);  // skip num_color_sets
            m_saver.append(cs_file);
            cs_file.close();

            assert(cdbg_builder.num_unitigs() > 0);
            assert(cdbg_builder.num_unitigs() <= UINT32_MAX);

            std::cout << "num_unitigs " << cdbg_builder.num_unitigs() << std::endl;
            std::cout << "num_distinct_color_sets " << cdbg_builder.num_color_sets() << std::endl;
        }

        {
            util::timed_phase("step 3. build unitig-to-color map");

            bits::bit_vector u2c;
            bits::rank9 u2c_rank1_index;
            essentials::load(u2c, cdbg_build_config.u2c_filename().c_str());
            u2c_rank1_index.build(u2c);
            m_saver.visit(u2c);
            m_saver.visit(u2c_rank1_index);

            assert(u2c.num_bits() == cdbg_builder.num_unitigs());
            assert(u2c_rank1_index.num_ones() == cdbg_builder.num_color_sets());

            std::cout << "m_u2c.num_bits() " << u2c.num_bits() << std::endl;
            std::cout << "m_u2c_rank1_index.num_ones() " << u2c_rank1_index.num_ones() << std::endl;
        }

        {
            util::timed_phase("step 4. build SSHash");

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
        }

        {
            util::timed_phase("step 5. write filenames");

            filenames filenames;
            filenames.build_from_file(m_build_config.filenames_list);
            m_saver.visit(filenames);
        }

        {
            m_saver.write(constants::current_version_number::major);
            m_saver.write(constants::current_version_number::minor);
            m_saver.write(constants::current_version_number::patch);
        }
    }

    void check(index const& idx) {
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        essentials::logger("checking correctness...");
        timer.start();
        std::atomic<uint64_t> num_checked_unitigs(0);

        /*
        m_ccdbg.loop_through_unitigs(
            [&](ggcat::Slice<char> const unitig,         //
                ggcat::Slice<uint32_t> const color_set,  //
                bool same_color_set)               //
            {
                auto lookup_result = idx.m_k2u.lookup(unitig.data);
                const uint64_t unitig_id = lookup_result.string_id;
                const uint64_t color_set_id = idx.u2c(unitig_id);
                for (uint64_t i = 1; i != unitig.size - idx.m_k2u.k() + 1; ++i) {
                    const uint64_t got = idx.m_k2u.lookup(unitig.data + i).string_id;
                    if (got != unitig_id) {
                        std::cout << "\033[1;31m"
                                  << "got unitig_id " << got << " but expected " << unitig_id
                                  << "\033[0m" << std::endl;
                        return;
                    }
                }
                auto fwd_it = idx.m_color_sets.color_set(color_set_id);
                const uint64_t size = fwd_it.size();
                if (size != color_set.size) {
                    std::cout << "\033[1;31m [" << color_set_id << "] "
                              << "got color_set size " << size << " but expected " << color_set.size
                              << "\033[0m" << std::endl;
                    return;
                }
                for (uint64_t i = 0; i != size; ++i, ++fwd_it) {
                    const uint32_t ref = *fwd_it;
                    if (ref != color_set.data[i]) {
                        std::cout << "\033[1;31m"
                                  << "got ref " << ref << " but expected " << color_set.data[i]
                                  << "\033[0m" << std::endl;
                        return;
                    }
                }

                if (++num_checked_unitigs % 1000 == 0) {
                    std::cout << "\rChecked " << num_checked_unitigs << "/"
                              << idx.m_k2u.num_strings() << " unitigs" << std::flush;
                }
            },
            m_build_config.num_threads  //
        );
        */

        std::cout << "\rChecked " << num_checked_unitigs << "/" << idx.m_k2u.num_strings()
                  << " unitigs" << std::endl;

        timer.stop();
        std::cout << "** checking correctness took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        essentials::logger("CHECK DONE!");
    }

private:
    build_configuration m_build_config;
    util::external_saver m_saver;
};

}  // namespace fulgor

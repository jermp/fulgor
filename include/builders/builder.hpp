#pragma once

#include "external/kmeans/include/kmeans.hpp"
#include "include/index.hpp"
#include "include/GGCAT.hpp"

namespace fulgor {

template <typename ColorSets>
struct index<ColorSets>::builder {
    builder() {}

    builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.num_kmers() != 0) throw std::runtime_error("index already built");

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            essentials::logger("step 1. building colored compacted dBG from GGCAT...");
            timer.start();
            m_ccdbg.build(m_build_config);
            m_build_config.num_colors = m_ccdbg.num_colors();
            timer.stop();
            std::cout << "** building the ccdBG took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        std::string input_filename_for_sshash = m_build_config.tmp_dirname + "/" +
                                                util::filename(m_build_config.file_base_name) +
                                                ".sshash.fa";

        {
            essentials::logger("step 2. building unitig-to-color map and encoding color sets...");
            timer.start();

            uint64_t num_unitigs = 0;
            uint64_t num_distinct_color_sets = 0;

            typename ColorSets::builder builder(m_build_config.num_colors, m_build_config);

            const uint64_t num_threads = m_build_config.num_threads;
            kmeans::thread_pool threads(num_threads);
            std::atomic<uint64_t> total_builders_bytes = 0;

            bits::bit_vector::builder u2c_builder;

            /* write unitigs to fasta file for SSHash */
            std::ofstream out(input_filename_for_sshash.c_str());
            if (!out.is_open()) throw std::runtime_error("cannot open output file");

            m_ccdbg.loop_through_unitigs([&](ggcat::Slice<char> const unitig,
                                             ggcat::Slice<uint32_t> const color_set,
                                             bool same_color_set) {
                try {
                    if (!same_color_set) {
                        if (num_unitigs > 0) u2c_builder.set(num_unitigs - 1, 1);

                        std::vector<uint32_t> cs(color_set.data, color_set.data + color_set.size);
                        threads.enqueue([&builder, cs = std::move(cs), num_distinct_color_sets]() mutable {
                            // for (auto c : cs) {
                            //     std::cout << c << " ";
                            // }
                            // std::cout << std::endl;
                            builder.encode_color_set(std::move(cs), num_distinct_color_sets);
                        });
                        num_distinct_color_sets += 1;
                    }
                    u2c_builder.push_back(0);

                    /*
                        Rewrite unitigs in color-set order.
                        This is *not* the same order in which
                        unitigs are written in the ggcat.fa file.
                    */
                    out << ">\n";
                    out.write(unitig.data, unitig.size);
                    out << '\n';

                    num_unitigs += 1;

                } catch (std::exception const& e) {
                    std::cerr << e.what() << std::endl;
                    exit(1);
                }
            });
            threads.wait();

            out.close();

            assert(num_unitigs > 0);
            assert(num_unitigs < (uint64_t(1) << 32));

            std::cout << "num_unitigs " << num_unitigs << std::endl;
            std::cout << "num_distinct_color_sets " << num_distinct_color_sets << std::endl;

            builder.build(idx.m_color_sets);

            timer.stop();
            std::cout << "** encoding color sets took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();

            timer.start();
            u2c_builder.set(num_unitigs - 1, 1);
            u2c_builder.build(idx.m_u2c);
            idx.m_u2c_rank1_index.build(idx.m_u2c);
            assert(idx.m_u2c.num_bits() == num_unitigs);
            assert(idx.m_u2c_rank1_index.num_ones() == num_distinct_color_sets);

            std::cout << "m_u2c.num_bits() " << idx.m_u2c.num_bits() << std::endl;
            std::cout << "m_u2c_rank1_index.num_ones() " << idx.m_u2c_rank1_index.num_ones()
                      << std::endl;

            timer.stop();
            std::cout << "** building unitig-to-color map took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 3. building SSHash...");
            timer.start();

            sshash::build_configuration sshash_config;
            sshash_config.k = m_build_config.k;
            sshash_config.m = m_build_config.m;
            sshash_config.canonical = true;
            sshash_config.verbose = m_build_config.verbose;
            sshash_config.tmp_dirname = m_build_config.tmp_dirname;
            sshash_config.num_threads = m_build_config.num_threads;
            sshash_config.print();
            idx.m_k2u.build(input_filename_for_sshash, sshash_config);
            try {  // remove unitig file
                std::remove(input_filename_for_sshash.c_str());
            } catch (std::exception const& e) { std::cerr << e.what() << std::endl; }

            timer.stop();
            std::cout << "** building SSHash took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 4. writing filenames...");
            timer.start();
            idx.m_filenames.build(m_ccdbg.filenames());
            timer.stop();
            std::cout << "** writing filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }
    }

    void check(index const& idx) {
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        essentials::logger("checking correctness...");
        timer.start();
        std::atomic<uint64_t> num_checked_unitigs(0);

        m_ccdbg.loop_through_unitigs(
            [&](ggcat::Slice<char> const unitig,         //
                ggcat::Slice<uint32_t> const color_set,  //
                bool /* same_color_set */)               //
            {
                auto lookup_result = idx.m_k2u.lookup(unitig.data);
                const uint64_t unitig_id = lookup_result.string_id;
                const uint64_t color_set_id = idx.u2c(unitig_id);
                for (uint64_t i = 1; i != unitig.size - idx.m_k2u.k() + 1; ++i) {
                    uint64_t got = idx.m_k2u.lookup(unitig.data + i).string_id;
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
                    uint32_t ref = *fwd_it;
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

        std::cout << "\rChecked " << num_checked_unitigs << "/"
                  << idx.m_k2u.num_strings() << " unitigs" << std::endl;

        timer.stop();
        std::cout << "** checking correctness took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        essentials::logger("CHECK DONE!");
    }

private:
    build_configuration m_build_config;
    GGCAT m_ccdbg;
};

}  // namespace fulgor

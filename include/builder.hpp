#pragma once

#include "index.hpp"
#include "GGCAT.hpp"

namespace fulgor {

template <typename ColorSets>
struct index<ColorSets>::builder {
    builder() {}

    builder(build_configuration const& build_config) : m_build_config(build_config) {}

    void build(index& idx) {
        if (idx.m_k2u.size() != 0) throw std::runtime_error("index already built");

        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

        {
            essentials::logger("step 1. build colored compacted dBG");
            timer.start();
            m_ccdbg.build(m_build_config);
            m_build_config.num_colors = m_ccdbg.num_colors();
            timer.stop();
            std::cout << "** building the ccdBG took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 2. build m_u2c and m_color_sets");
            timer.start();

            uint64_t num_unitigs = 0;
            uint64_t num_distinct_colors = 0;

            pthash::bit_vector_builder u2c_builder;

            /* write unitigs to fasta file for SSHash */
            std::ofstream out((m_build_config.file_base_name + ".fa").c_str());
            if (!out.is_open()) throw std::runtime_error("cannot open output file");

            typename ColorSets::builder colors_builder(m_build_config.num_colors);

            m_ccdbg.loop_through_unitigs([&](ggcat::Slice<char> const unitig,
                                             ggcat::Slice<uint32_t> const colors, bool same_color) {
                try {
                    if (!same_color) {
                        num_distinct_colors += 1;
                        if (num_unitigs > 0) u2c_builder.set(num_unitigs - 1, 1);

                        /* compress colors */
                        colors_builder.process(colors.data, colors.size);
                    }
                    u2c_builder.push_back(0);

                    out << ">\n";
                    out.write(unitig.data, unitig.size);
                    out << '\n';

                    num_unitigs += 1;

                } catch (std::exception const& e) {
                    std::cerr << e.what() << std::endl;
                    exit(1);
                }
            });

            out.close();

            assert(num_unitigs > 0);
            assert(num_unitigs < (uint64_t(1) << 32));

            std::cout << "num_unitigs " << num_unitigs << std::endl;
            std::cout << "num_distinct_colors " << num_distinct_colors << std::endl;

            u2c_builder.set(num_unitigs - 1, 1);
            idx.m_u2c.build(&u2c_builder);
            assert(idx.m_u2c.size() == num_unitigs);
            assert(idx.m_u2c.num_ones() == num_distinct_colors);

            std::cout << "m_u2c.size() " << idx.m_u2c.size() << std::endl;
            std::cout << "m_u2c.num_ones() " << idx.m_u2c.num_ones() << std::endl;
            std::cout << "m_u2c.num_zeros() " << idx.m_u2c.num_zeros() << std::endl;

            colors_builder.build(idx.m_color_sets);

            timer.stop();
            std::cout << "** building m_u2c and m_color_sets took " << timer.elapsed()
                      << " seconds / " << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 3. build m_k2u");
            timer.start();

            sshash::build_configuration sshash_config;
            sshash_config.k = m_build_config.k;
            sshash_config.m = m_build_config.m;
            sshash_config.canonical_parsing = m_build_config.canonical_parsing;
            sshash_config.verbose = m_build_config.verbose;
            sshash_config.tmp_dirname = m_build_config.tmp_dirname;
            sshash_config.print();
            idx.m_k2u.build(m_build_config.file_base_name + ".fa", sshash_config);
            try {  // remove unitig file
                std::remove((m_build_config.file_base_name + ".fa").c_str());
            } catch (std::exception const& e) { std::cerr << e.what() << std::endl; }

            timer.stop();
            std::cout << "** building m_k2u took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        {
            essentials::logger("step 4. write filenames");
            timer.start();
            idx.m_filenames.build(m_ccdbg.filenames());
            timer.stop();
            std::cout << "** writing filenames took " << timer.elapsed() << " seconds / "
                      << timer.elapsed() / 60 << " minutes" << std::endl;
            timer.reset();
        }

        if (m_build_config.check)  //
        {
            essentials::logger("step 5. check correctness...");
            m_ccdbg.loop_through_unitigs(
                [&](ggcat::Slice<char> const unitig, ggcat::Slice<uint32_t> const colors,
                    bool /* same_color */)  //
                {
                    auto lookup_result = idx.m_k2u.lookup_advanced(unitig.data);
                    const uint64_t unitig_id = lookup_result.contig_id;
                    const uint64_t color_id = idx.u2c(unitig_id);
                    for (uint64_t i = 1; i != unitig.size - idx.m_k2u.k() + 1; ++i) {
                        uint64_t got = idx.m_k2u.lookup_advanced(unitig.data + i).contig_id;
                        if (got != unitig_id) {
                            std::cout << "got unitig_id " << got << " but expected " << unitig_id
                                      << std::endl;
                            return;
                        }
                    }
                    auto fwd_it = idx.m_color_sets.color_set(color_id);
                    const uint64_t size = fwd_it.size();
                    if (size != colors.size) {
                        std::cout << "got colors list of size " << size << " but expected "
                                  << colors.size << std::endl;
                        return;
                    }
                    for (uint64_t i = 0; i != size; ++i, ++fwd_it) {
                        uint32_t ref = *fwd_it;
                        if (ref != colors.data[i]) {
                            std::cout << "got ref " << ref << " but expected " << colors.data[i]
                                      << std::endl;
                            return;
                        }
                    }
                },
                m_build_config.num_threads  //
            );
            essentials::logger("DONE!");
        }
    }

private:
    build_configuration m_build_config;
    GGCAT m_ccdbg;
};

}  // namespace fulgor

#include "index.hpp"

namespace fulgor {

template <typename ColorClasses>
void index<ColorClasses>::build_from(ccdbg_builder const& builder) {
    if (m_k2u.size() != 0) throw std::runtime_error("index already built");

    {
        essentials::logger("step 1. build m_u2c, m_ccs, and m_k2u");

        uint64_t num_unitigs = 0;
        uint64_t num_distinct_colors = 0;

        pthash::bit_vector_builder bvb;  // for m_u2c

        sshash::build_configuration sshash_config;
        sshash_config.k = builder.config.k;
        sshash_config.m = builder.config.m;
        sshash_config.canonical_parsing = builder.config.canonical_parsing;
        sshash_config.verbose = builder.config.verbose;
        sshash_config.tmp_dirname = builder.config.tmp_dirname;
        sshash_config.print();
        sshash::builder k2u_builder(sshash_config);

        typename ColorClasses::builder colors_builder(builder.config);

        builder.ggcat->loop_through_unitigs([&](ggcat::Slice<char> const unitig,
                                                ggcat::Slice<uint32_t> const colors,
                                                bool same_color) {
            try {
                if (!same_color) {
                    num_distinct_colors += 1;
                    if (num_unitigs > 0) bvb.set(num_unitigs - 1, 1);

                    /* compress colors */
                    colors_builder.process(colors.data, colors.size);
                }
                bvb.push_back(0);

                k2u_builder.add_sequence(unitig.data, unitig.size);
                num_unitigs += 1;

            } catch (std::exception const& e) {
                std::cerr << e.what() << std::endl;
                exit(1);
            }
        });

        assert(num_unitigs < (uint64_t(1) << 32));
        std::cout << "num_unitigs " << num_unitigs << std::endl;
        std::cout << "num_distinct_colors " << num_distinct_colors << std::endl;

        m_u2c.build(&bvb);
        std::cout << "m_u2c.size() " << m_u2c.size() << std::endl;
        std::cout << "m_u2c.num_ones() " << m_u2c.num_ones() << std::endl;
        std::cout << "m_u2c.num_zeros() " << m_u2c.num_zeros() << std::endl;

        m_k2u.build_from(k2u_builder, sshash_config);

        colors_builder.build(m_ccs);
    }

    {
        essentials::logger("step 3. write filenames");
        m_filenames.build(builder.ggcat->filenames());
    }

    if (builder.config.check) {
        essentials::logger("step 4. check correctness...");
        builder.ggcat->loop_through_unitigs(
            [&](ggcat::Slice<char> const unitig, ggcat::Slice<uint32_t> const colors,
                bool /* same_color */) {
                auto lookup_result = m_k2u.lookup_advanced(unitig.data);
                uint32_t unitig_id = lookup_result.contig_id;
                uint32_t color_id = u2c(unitig_id);
                for (uint64_t i = 1; i != unitig.size - m_k2u.k() + 1; ++i) {
                    uint32_t got = m_k2u.lookup_advanced(unitig.data + i).contig_id;
                    if (got != unitig_id) {
                        std::cout << "got unitig_id " << got << " but expected " << unitig_id
                                  << std::endl;
                        return;
                    }
                }
                auto fwd_it = m_ccs.colors(color_id);
                uint64_t size = fwd_it.size();
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
            builder.config.num_threads);
        essentials::logger("DONE!");
    }
}

template <typename ColorClasses>
void index<ColorClasses>::print_stats() const {
    uint64_t total_bits = num_bits();
    std::cout << "total index size: " << essentials::convert(total_bits / 8, essentials::GB)
              << " GB" << '\n';
    std::cout << "SPACE BREAKDOWN:\n";
    std::cout << "  K2U: " << m_k2u.num_bits() / 8 << " bytes / "
              << essentials::convert(m_k2u.num_bits() / 8, essentials::GB) << " GB ("
              << (m_k2u.num_bits() * 100.0) / total_bits << "%)\n";
    std::cout << "  CCs: " << m_ccs.num_bits() / 8 << " bytes / "
              << essentials::convert(m_ccs.num_bits() / 8, essentials::GB) << " GB ("
              << (m_ccs.num_bits() * 100.0) / total_bits << "%)\n";
    uint64_t other_bits = m_u2c.bytes() * 8 + m_filenames.num_bits();
    std::cout << "  Other: " << other_bits / 8 << " bytes / "
              << essentials::convert(other_bits / 8, essentials::GB) << " GB ("
              << (other_bits * 100.0) / total_bits << "%)\n";
    std::cout << "    U2C: " << m_u2c.bytes() << " bytes / "
              << essentials::convert(m_u2c.bytes(), essentials::GB) << " GB ("
              << (m_u2c.bytes() * 8 * 100.0) / total_bits << "%)\n";
    std::cout << "    filenames: " << m_filenames.num_bits() / 8 << " bytes / "
              << essentials::convert(m_filenames.num_bits() / 8, essentials::GB) << " GB ("
              << (m_filenames.num_bits() * 100.0) / total_bits << "%)\n";

    uint64_t num_ints_in_ccs = 0;
    uint64_t num_ccs = m_ccs.num_color_classes();
    std::cout << "Color id range 0.." << num_docs() - 1 << '\n';
    std::cout << "Number of distinct color classes: " << num_ccs << '\n';
    for (uint64_t color_class_id = 0; color_class_id != num_ccs; ++color_class_id) {
        uint64_t list_size = m_ccs.colors(color_class_id).size();
        num_ints_in_ccs += list_size;
    }
    std::cout << "Number of ints in distinct color classes: " << num_ints_in_ccs << " ("
              << static_cast<double>(m_ccs.num_bits()) / num_ints_in_ccs << " bits/int)\n";
    std::cout << "k: " << m_k2u.k() << '\n';
    std::cout << "m: " << m_k2u.m() << " (minimizer length used in K2U)\n";
    std::cout << "Number of kmers in dBG: " << m_k2u.size() << " ("
              << static_cast<double>(m_k2u.num_bits()) / m_k2u.size() << " bits/kmer)\n";
    std::cout << "Number of unitigs in dBG: " << m_k2u.num_contigs() << std::endl;

    m_ccs.print_stats();
}

}  // namespace fulgor
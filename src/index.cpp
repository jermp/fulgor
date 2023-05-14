#include "index.hpp"

namespace fulgor {

template <typename ColorClasses>
void index<ColorClasses>::build(build_configuration const& build_config) {
    if (m_k2u.size() != 0) throw std::runtime_error("index already built");

    {
        essentials::logger("step 1. build m_u2c and m_ccs");

        uint64_t num_unitigs = 0;
        uint64_t num_distinct_colors = 0;

        pthash::bit_vector_builder bvb;  // for m_u2c

        /* write unitigs to fasta file for SSHash */
        std::ofstream out((build_config.file_base_name + ".fa").c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open output file");

        typename ColorClasses::builder colors_builder(build_config);

        build_config.ggcat->loop_through_unitigs([&](ggcat::Slice<char> const unitig,
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

        assert(num_unitigs < (uint64_t(1) << 32));
        std::cout << "num_unitigs " << num_unitigs << std::endl;
        std::cout << "num_distinct_colors " << num_distinct_colors << std::endl;

        m_u2c.build(&bvb);
        std::cout << "m_u2c.size() " << m_u2c.size() << std::endl;
        std::cout << "m_u2c.num_ones() " << m_u2c.num_ones() << std::endl;
        std::cout << "m_u2c.num_zeros() " << m_u2c.num_zeros() << std::endl;

        colors_builder.build(m_ccs);
    }

    {
        essentials::logger("step 2. build m_k2u");
        sshash::build_configuration sshash_build_config;
        sshash_build_config.k = build_config.k;
        sshash_build_config.m = build_config.m;
        sshash_build_config.canonical_parsing = build_config.canonical_parsing;
        sshash_build_config.verbose = build_config.verbose;
        sshash_build_config.tmp_dirname = build_config.tmp_dirname;
        sshash_build_config.print();
        m_k2u.build(build_config.file_base_name + ".fa", sshash_build_config);
    }

    {
        essentials::logger("step 3. write filenames");
        m_filenames.build(build_config.ggcat->filenames());
    }

    if (build_config.check) {
        essentials::logger("step 4. check correctness...");
        build_config.ggcat->loop_through_unitigs([&](ggcat::Slice<char> const unitig,
                                                     ggcat::Slice<uint32_t> const colors,
                                                     bool /* same_color */) {
            auto lookup_result = m_k2u.lookup_advanced(unitig.data);
            uint32_t unitig_id = lookup_result.contig_id;
            uint32_t color_id = u2c(unitig_id);
            for (uint64_t i = 1; i != unitig.size - m_k2u.k() + 1; ++i) {
                uint32_t got = m_k2u.lookup_advanced(unitig.data + i).contig_id;
                if (got != unitig_id) {
                    std::cout << "got unitig_id " << got << " but expected " << unitig_id
                              << std::endl;
                }
            }
            auto fwd_it = m_ccs.colors(color_id);
            uint64_t size = fwd_it.size();
            if (size != colors.size) {
                std::cout << "got colors list of size " << size << " but expected " << colors.size
                          << std::endl;
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
        });
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
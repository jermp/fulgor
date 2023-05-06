#include "index.hpp"

namespace fulgor {

template <typename ColorClasses>
void index<ColorClasses>::build(build_configuration const& build_config) {
    if (m_k2u.size() != 0) throw std::runtime_error("index already built");

    {
        essentials::logger("step 1. build m_k2u");
        sshash::build_configuration sshash_build_config;
        sshash_build_config.k = build_config.k;
        sshash_build_config.m = build_config.m;
        sshash_build_config.canonical_parsing = build_config.canonical_parsing;
        sshash_build_config.verbose = build_config.verbose;
        sshash_build_config.tmp_dirname = build_config.tmp_dirname;
        sshash_build_config.print();
        m_k2u.build(build_config.file_base_name + ".permuted.fa", sshash_build_config);
    }

    {
        essentials::logger("step 2. build m_u2c");
        mm::file_source<seg_id_t> in(build_config.file_base_name + ".map", mm::advice::sequential);
        uint64_t num_sequences = in.bytes() / (sizeof(seg_id_t) + sizeof(uint64_t));
        pthash::bit_vector_builder bvb(num_sequences);
        seg_id_t const* data = in.data();
        assert(num_sequences < (uint64_t(1) << 32));
        uint64_t prev_color_class_id = 0;
        for (uint64_t i = 0; i != num_sequences; ++i) {
            data += 1;  // skip seg_id
            uint64_t color_class_id = 0;
            if constexpr (sizeof(seg_id_t) == 4) {
                color_class_id = *reinterpret_cast<uint64_t const*>(data);
                data += 2;
            } else {
                assert(sizeof(seg_id_t) == 8);
                color_class_id = *data;
                data += 1;
            }
            if (color_class_id != prev_color_class_id) {
                assert(i > 0);
                bvb.set(i - 1, 1);
            }
            prev_color_class_id = color_class_id;
        }
        m_u2c.build(&bvb);
        std::cout << "m_u2c.size() " << m_u2c.size() << std::endl;
        std::cout << "m_u2c.num_ones() " << m_u2c.num_ones() << std::endl;
        std::cout << "m_u2c.num_zeros() " << m_u2c.num_zeros() << std::endl;
    }

    {
        essentials::logger("step 2.1. check correctness of m_u2c");
        mm::file_source<seg_id_t> in(build_config.file_base_name + ".map", mm::advice::sequential);
        uint64_t num_sequences = in.bytes() / (sizeof(seg_id_t) + sizeof(uint64_t));
        seg_id_t const* data = in.data();
        for (uint64_t i = 0; i != num_sequences; ++i) {
            data += 1;  // skip seg_id
            uint64_t color = 0;
            if constexpr (sizeof(seg_id_t) == 4) {
                color = *reinterpret_cast<uint64_t const*>(data);
                data += 2;
            } else {
                assert(sizeof(seg_id_t) == 8);
                color = *data;
                data += 1;
            }
            uint64_t got = m_u2c.rank(i);
            if (got != color) {
                std::cout << "Error: expected color " << color << " but got " << got << std::endl;
            }
        }
    }

    {
        essentials::logger("step 3. build colors");
        m_ccs.build(build_config);
    }

    {
        essentials::logger("step 4. write filenames");
        std::vector<std::string> filenames;
        uint64_t n = num_docs();
        filenames.reserve(n);
        std::ifstream in(build_config.file_base_name + ".json");
        if (!in.is_open()) throw std::runtime_error("error in opening file");
        std::string line;
        std::getline(in, line, ':');
        std::getline(in, line, ':');
        in.get();  // skip '"'
        for (uint64_t i = 0; i != n; ++i) {
            in.get();  // skip ' '
            char delim = ',';
            if (i == n - 1) {
                delim = '"';  // last
            }
            std::getline(in, line, delim);
            filenames.push_back(line);
        }
        m_filenames.build(filenames);
        // m_filenames.print();
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
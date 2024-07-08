#include "include/index.hpp"

namespace fulgor {

template <typename ColorClasses>
void index<ColorClasses>::print_stats() const {
    const uint64_t total_bits = num_bits();
    auto const& k2u = get_k2u();
    auto const& u2c = get_u2c();
    auto const& ccs = get_color_sets();
    auto const& filenames = get_filenames();

    std::cout << "total index size: " << total_bits / 8 << " [B] -- "
              << essentials::convert(total_bits / 8, essentials::GB) << " [GB]" << '\n';
    std::cout << "SPACE BREAKDOWN:\n";
    std::cout << "  K2U: " << k2u.num_bits() / 8 << " bytes / "
              << essentials::convert(k2u.num_bits() / 8, essentials::GB) << " GB ("
              << (k2u.num_bits() * 100.0) / total_bits << "%)\n";
    std::cout << "  CCs: " << ccs.num_bits() / 8 << " bytes / "
              << essentials::convert(ccs.num_bits() / 8, essentials::GB) << " GB ("
              << (ccs.num_bits() * 100.0) / total_bits << "%)\n";
    uint64_t other_bits = u2c.bytes() * 8 + filenames.num_bits();
    std::cout << "  Other: " << other_bits / 8 << " bytes / "
              << essentials::convert(other_bits / 8, essentials::GB) << " GB ("
              << (other_bits * 100.0) / total_bits << "%)\n";
    std::cout << "    U2C: " << u2c.bytes() << " bytes / "
              << essentials::convert(u2c.bytes(), essentials::GB) << " GB ("
              << (u2c.bytes() * 8 * 100.0) / total_bits << "%)\n";
    std::cout << "    filenames: " << filenames.num_bits() / 8 << " bytes / "
              << essentials::convert(filenames.num_bits() / 8, essentials::GB) << " GB ("
              << (filenames.num_bits() * 100.0) / total_bits << "%)\n";

    uint64_t num_ints_in_ccs = 0;
    uint64_t num_ccs = ccs.num_color_sets();
    std::cout << "Color id range 0.." << num_docs() - 1 << '\n';
    std::cout << "Number of distinct color classes: " << num_ccs << '\n';
    for (uint64_t color_set_id = 0; color_set_id != num_ccs; ++color_set_id) {
        uint64_t list_size = ccs.color_set(color_set_id).size();
        num_ints_in_ccs += list_size;
    }
    std::cout << "Number of ints in distinct color classes: " << num_ints_in_ccs << " ("
              << static_cast<double>(ccs.num_bits()) / num_ints_in_ccs << " bits/int)\n";
    std::cout << "k: " << k2u.k() << '\n';
    std::cout << "m: " << k2u.m() << " (minimizer length used in K2U)\n";
    std::cout << "Number of kmers in dBG: " << k2u.size() << " ("
              << static_cast<double>(k2u.num_bits()) / k2u.size() << " bits/kmer)\n";
    std::cout << "Number of unitigs in dBG: " << k2u.num_contigs() << std::endl;

    ccs.print_stats();
}

// template <typename ColorClasses>
// void index<ColorClasses>::dump_colors(std::ofstream& os) const {
//     os << "num_references " << num_docs() << '\n';
//     auto const& ccs = get_color_classes();
//     ccs.dump(os);
// }

template <typename ColorClasses>
void index<ColorClasses>::dump(std::string const& basename) const {
    /* metadata file */
    std::ofstream metadata(basename + ".metadata.txt");
    if (!metadata.is_open()) throw std::runtime_error("cannot open output file");
    metadata << "num_references=" << num_docs() << '\n';
    metadata << "num_unitigs=" << num_unitigs() << '\n';
    metadata << "num_color_classes=" << num_color_sets() << '\n';
    metadata.close();

    /* unitigs file */
    std::ofstream unitigs(basename + ".unitigs.fa");
    if (!unitigs.is_open()) throw std::runtime_error("cannot open output file");
    const uint64_t u = num_unitigs();
    const uint64_t kmer_length = k();
    for (uint64_t unitig_id = 0; unitig_id != u; ++unitig_id) {
        auto it = m_k2u.at_contig_id(unitig_id);
        const uint64_t color_id = u2c(unitig_id);
        unitigs << "> unitig_id=" << unitig_id << " color_id=" << color_id << '\n';
        auto [_, kmer] = it.next();
        unitigs << kmer;
        while (it.has_next()) {
            auto [_, kmer] = it.next();
            unitigs << kmer[kmer_length - 1];  // overlaps!
        }
        unitigs << '\n';
    }
    unitigs.close();

    /* colors file */
    std::ofstream colors(basename + ".colors.txt");
    if (!colors.is_open()) throw std::runtime_error("cannot open output file");
    auto const& ccs = get_color_sets();
    const uint64_t n = num_color_sets();
    for (uint64_t color_id = 0; color_id != n; ++color_id) {
        auto it = ccs.color_set(color_id);
        const uint32_t size = it.size();
        colors << "color_id=" << color_id << " size=" << size << ' ';
        for (uint32_t j = 0; j != size; ++j) {
            colors << it.value();
            it.next();
            if (j != size - 1) colors << ' ';
        }
        colors << '\n';
    }
    colors.close();
}

}  // namespace fulgor
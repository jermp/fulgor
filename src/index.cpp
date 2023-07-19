#include "index.hpp"

namespace fulgor {

template <typename ColorClasses>
void index<ColorClasses>::print_stats() const {
    uint64_t total_bits = num_bits();
    auto const& k2u = get_dict();
    auto const& ccs = color_classes();
    auto const& u2c = get_u2c();
    auto const& filenames = get_filenames();

    std::cout << "total index size: " << essentials::convert(total_bits / 8, essentials::GB)
              << " GB" << '\n';
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
    uint64_t num_ccs = ccs.num_color_classes();
    std::cout << "Color id range 0.." << num_docs() - 1 << '\n';
    std::cout << "Number of distinct color classes: " << num_ccs << '\n';
    for (uint64_t color_class_id = 0; color_class_id != num_ccs; ++color_class_id) {
        uint64_t list_size = ccs.colors(color_class_id).size();
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

}  // namespace fulgor
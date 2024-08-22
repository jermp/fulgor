#include "include/index.hpp"

namespace fulgor {

template <typename ColorSets>
void index<ColorSets>::print_stats() const {
    const uint64_t total_bits = num_bits();
    assert(total_bits > 0);
    auto const& k2u = get_k2u();
    auto const& u2c = get_u2c();
    auto const& color_sets = get_color_sets();
    auto const& filenames = get_filenames();

    std::cout << "total index size: " << total_bits / 8 << " [B] -- "
              << essentials::convert(total_bits / 8, essentials::GB) << " [GB]" << '\n';
    std::cout << "SPACE BREAKDOWN:\n";
    std::cout << "  dBG (SSHash): " << k2u.num_bits() / 8 << " bytes / "
              << essentials::convert(k2u.num_bits() / 8, essentials::GB) << " GB ("
              << (k2u.num_bits() * 100.0) / total_bits << "%)\n";
    std::cout << "  Color sets: " << color_sets.num_bits() / 8 << " bytes / "
              << essentials::convert(color_sets.num_bits() / 8, essentials::GB) << " GB ("
              << (color_sets.num_bits() * 100.0) / total_bits << "%)\n";
    uint64_t other_bits = u2c.bytes() * 8 + filenames.num_bits();
    std::cout << "  Other: " << other_bits / 8 << " bytes / "
              << essentials::convert(other_bits / 8, essentials::GB) << " GB ("
              << (other_bits * 100.0) / total_bits << "%)\n";
    std::cout << "    Map from unitig_id to color_set_id: " << u2c.bytes() << " bytes / "
              << essentials::convert(u2c.bytes(), essentials::GB) << " GB ("
              << (u2c.bytes() * 8 * 100.0) / total_bits << "%)\n";
    std::cout << "    filenames: " << filenames.num_bits() / 8 << " bytes / "
              << essentials::convert(filenames.num_bits() / 8, essentials::GB) << " GB ("
              << (filenames.num_bits() * 100.0) / total_bits << "%)\n";

    uint64_t num_ints_in_color_sets = 0;
    uint64_t num_color_sets = color_sets.num_color_sets();
    std::cout << "Color id range 0.." << num_colors() - 1 << '\n';
    std::cout << "Number of distinct color sets: " << num_color_sets << '\n';
    for (uint64_t color_set_id = 0; color_set_id != num_color_sets; ++color_set_id) {
        uint64_t list_size = color_sets.color_set(color_set_id).size();
        num_ints_in_color_sets += list_size;
    }
    std::cout << "Number of ints in distinct color sets: " << num_ints_in_color_sets << " ("
              << static_cast<double>(color_sets.num_bits()) / num_ints_in_color_sets
              << " bits/int)\n";
    std::cout << "k: " << k2u.k() << '\n';
    std::cout << "m: " << k2u.m() << " (minimizer length used in K2U)\n";
    std::cout << "Number of kmers in dBG: " << k2u.size() << " ("
              << static_cast<double>(k2u.num_bits()) / k2u.size() << " bits/kmer)\n";
    std::cout << "Number of unitigs in dBG: " << k2u.num_contigs() << std::endl;

    color_sets.print_stats();
}

template <typename ColorSets>
void index<ColorSets>::dump(std::string const& basename) const {
    /* metadata file */
    std::ofstream metadata_file(basename + ".metadata.txt");
    if (!metadata_file.is_open()) throw std::runtime_error("cannot open output file");
    metadata_file << "k=" << k() << '\n';
    metadata_file << "num_colors=" << num_colors() << '\n';
    metadata_file << "num_unitigs=" << num_unitigs() << '\n';
    metadata_file << "num_color_sets=" << num_color_sets() << '\n';
    metadata_file.close();

    /* unitigs file */
    std::ofstream unitigs_file(basename + ".unitigs.fa");
    if (!unitigs_file.is_open()) throw std::runtime_error("cannot open output file");
    const uint64_t u = num_unitigs();
    const uint64_t kmer_length = k();
    for (uint64_t unitig_id = 0; unitig_id != u; ++unitig_id) {
        auto it = m_k2u.at_contig_id(unitig_id);
        const uint64_t color_set_id = u2c(unitig_id);
        unitigs_file << "> unitig_id=" << unitig_id << " color_set_id=" << color_set_id << '\n';
        auto [_, kmer] = it.next();
        unitigs_file << kmer;
        while (it.has_next()) {
            auto [_, kmer] = it.next();
            unitigs_file << kmer[kmer_length - 1];  // overlaps!
        }
        unitigs_file << '\n';
    }
    unitigs_file.close();

    /* color_sets file */
    std::ofstream color_sets_file(basename + ".color_sets.txt");
    if (!color_sets_file.is_open()) throw std::runtime_error("cannot open output file");
    auto const& color_sets = get_color_sets();
    const uint64_t n = num_color_sets();
    for (uint64_t color_set_id = 0; color_set_id != n; ++color_set_id) {
        auto it = color_sets.color_set(color_set_id);
        const uint32_t size = it.size();
        color_sets_file << "color_set_id=" << color_set_id << " size=" << size << ' ';
        for (uint32_t j = 0; j != size; ++j) {
            color_sets_file << it.value();
            it.next();
            if (j != size - 1) color_sets_file << ' ';
        }
        color_sets_file << '\n';
    }
    color_sets_file.close();
}

}  // namespace fulgor
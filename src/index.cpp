#include "include/index.hpp"

namespace fulgor {

template <typename ColorSets>
void index<ColorSets>::print_stats() const {
    const uint64_t total_bits = num_bits();
    assert(total_bits > 0);
    auto const& k2u = get_k2u();
    auto const& u2c = get_u2c();
    auto const& u2c_rank1_index = get_u2c_rank1_index();
    auto const& color_sets = get_color_sets();
    auto const& filenames = get_filenames();

    uint64_t k2u_bytes = k2u.num_bits() / 8;
    uint64_t color_bytes = color_sets.num_bits() / 8;
    uint64_t u2c_bytes = u2c.num_bytes() + u2c_rank1_index.num_bytes();
    uint64_t filenames_bytes = filenames.num_bits() / 8;
    uint64_t other_bytes = u2c_bytes + filenames_bytes;

    uint64_t max_bytes = total_bits / 8;
    double max_gb = essentials::convert(max_bytes, essentials::GB);
    int byte_width = std::to_string(max_bytes).length();
    int gb_width = std::to_string(static_cast<uint64_t>(max_gb)).length() + 4;
    std::cout << std::format("{:<42} {:>{}} B / {:>{}.3f} GB\n",
                             "Total index size:", total_bits / 8, byte_width,
                             essentials::convert(total_bits / 8, essentials::GB), gb_width);
    std::cout << "SPACE BREAKDOWN:\n";
    std::cout << std::format("{:<42} {:>{}} B / {:>{}.3f} GB ({:>6.2f}%)\n",
                             "  dBG (SSHash):", k2u_bytes, byte_width,
                             essentials::convert(k2u_bytes, essentials::GB), gb_width,
                             k2u.num_bits() * 100.0 / total_bits);
    std::cout << std::format("{:<42} {:>{}} B / {:>{}.3f} GB ({:>6.2f}%)\n",
                             "  Color sets:", color_bytes, byte_width,
                             essentials::convert(color_bytes, essentials::GB), gb_width,
                             color_sets.num_bits() * 100.0 / total_bits);
    std::cout << std::format("{:<42} {:>{}} B / {:>{}.3f} GB ({:>6.2f}%)\n",
                             "  Other:", other_bytes, byte_width,
                             essentials::convert(other_bytes, essentials::GB), gb_width,
                             other_bytes * 8.0 * 100.0 / total_bits);
    std::cout << std::format("{:<42} {:>{}} B / {:>{}.3f} GB ({:>6.2f}%)\n",
                             "    Map from unitig_id to color_set_id:", u2c_bytes, byte_width,
                             essentials::convert(u2c_bytes, essentials::GB), gb_width,
                             u2c_bytes * 8.0 * 100.0 / total_bits);
    std::cout << std::format("{:<42} {:>{}} B / {:>{}.3f} GB ({:>6.2f}%)\n",
                             "    filenames:", filenames_bytes, byte_width,
                             essentials::convert(filenames_bytes, essentials::GB), gb_width,
                             filenames.num_bits() * 100.0 / total_bits);

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
    std::cout << "m: " << k2u.m() << " (minimizer length used in SSHash)\n";
    std::cout << "Number of kmers in dBG: " << k2u.num_kmers() << " ("
              << static_cast<double>(k2u.num_bits()) / k2u.num_kmers() << " bits/kmer)\n";
    std::cout << "Number of unitigs in dBG: " << k2u.num_strings() << std::endl;

    color_sets.print_stats();
}

template <typename ColorSets>
void index<ColorSets>::dump(build_configuration const& build_config) const  //
{
    /* metadata file */
    std::filesystem::path metadata_fn = build_config.output_filename;
    metadata_fn.replace_extension("metadata.txt");
    std::filesystem::path unitigs_fn = build_config.output_filename;
    unitigs_fn.replace_extension("unitigs.fa");
    std::filesystem::path color_sets_fn = build_config.output_filename;
    color_sets_fn.replace_extension("color_sets.fa");
    std::filesystem::path filenames_fn = build_config.output_filename;
    filenames_fn.replace_extension("filenames.txt");

    essentials::logger(std::format("writing '{}'...", metadata_fn.string()));
    std::ofstream metadata_file(metadata_fn);
    if (!metadata_file.is_open()) throw std::runtime_error("cannot open output file");
    metadata_file << "k=" << k() << '\n';
    metadata_file << "num_kmers=" << num_kmers() << '\n';
    metadata_file << "num_colors=" << num_colors() << '\n';
    metadata_file << "num_unitigs=" << num_unitigs() << '\n';
    metadata_file << "num_color_sets=" << num_color_sets() << '\n';
    metadata_file.close();

    /* filenames file */
    essentials::logger(std::format("writing '{}'...", filenames_fn.string()));
    std::ofstream filenames_file(filenames_fn);
    if (!filenames_file.is_open()) throw std::runtime_error("cannot open output file");
    for (uint64_t i = 0; i != num_colors(); ++i) filenames_file << filename(i) << '\n';
    filenames_file.close();

    /* unitigs file */
    essentials::logger(std::format("writing '{}'...", unitigs_fn.string()));
    std::ofstream unitigs_file(unitigs_fn);
    if (!unitigs_file.is_open()) throw std::runtime_error("cannot open output file");
    const uint64_t u = num_unitigs();
    const uint64_t kmer_length = k();
    std::string kmer(kmer_length, 0);
    for (uint64_t unitig_id = 0; unitig_id != u; ++unitig_id) {
        auto it = m_k2u.at_string_id(unitig_id);
        const uint64_t color_set_id = u2c(unitig_id);
        unitigs_file << "> color_set_id=" << color_set_id << '\n';
        auto [_, uint_kmer] = it.next();
        sshash::util::uint_kmer_to_string<kmer_type>(uint_kmer, kmer.data(), kmer_length);
        unitigs_file << kmer;
        while (it.has_next()) {
            auto [_, uint_kmer] = it.next();
            unitigs_file << kmer_type::uint64_to_char(uint_kmer.at(kmer_length - 1));  // overlaps!
        }
        unitigs_file << '\n';
    }
    unitigs_file.close();

    /* color_sets file */
    essentials::logger(std::format("writing '{}'...", color_sets_fn.string()));
    std::ofstream color_sets_file(color_sets_fn);
    if (!color_sets_file.is_open()) throw std::runtime_error("cannot open output file");
    auto const& color_sets = get_color_sets();
    const uint64_t n = num_color_sets();
    for (uint64_t color_set_id = 0; color_set_id != n; ++color_set_id) {
        auto it = color_sets.color_set(color_set_id);
        const uint32_t size = it.size();
        color_sets_file << "size=" << size << ' ';
        for (uint32_t j = 0; j != size; ++j) {
            color_sets_file << it.value();
            it.next();
            if (j != size - 1) color_sets_file << ' ';
        }
        color_sets_file << '\n';
    }
    color_sets_file.close();

    essentials::logger("DONE");
}
template <typename ColorSets>
void index<ColorSets>::loader::load_metadata() {
    util::timed_phase timer("step 1. read metadata");

    std::filesystem::path metadata_fn = m_config.base_filename;
    metadata_fn.replace_extension("metadata.txt");
    std::ifstream in(metadata_fn.c_str());
    if (!in.is_open()) throw std::runtime_error("cannot open metadata file");

    std::string line;
    while (std::getline(in, line)) {
        const size_t delimiter_pos = line.find('=');
        assert(delimiter_pos != std::string::npos);
        const std::string_view key(line.data(), delimiter_pos);
        char const* value_ptr = line.c_str() + delimiter_pos + 1;
        if (key == "k") {
            m_k = static_cast<uint32_t>(std::strtoul(value_ptr, nullptr, 10));
        } else if (key == "num_kmers") {
            m_num_kmers = std::strtoull(value_ptr, nullptr, 10);
        } else if (key == "num_colors") {
            m_num_colors = static_cast<uint32_t>(std::strtoul(value_ptr, nullptr, 10));
        } else if (key == "num_unitigs") {
            m_num_unitigs = std::strtoull(value_ptr, nullptr, 10);
        } else if (key == "num_color_sets") {
            m_num_color_sets = std::strtoull(value_ptr, nullptr, 10);
        }
    }
    in.close();

    if (m_config.verbose) {
        std::cout << "k=" << m_k << ", num_kmers=" << m_num_kmers << ", num_colors=" << m_num_colors
                  << ", num_unitigs=" << m_num_unitigs << ", num_color_sets=" << m_num_color_sets
                  << std::endl;
    }

    assert(m_num_unitigs > 0);
    assert(m_num_unitigs < (uint64_t(1) << 32));
}

template <typename ColorSets>
void index<ColorSets>::loader::load_u2c() {
    util::timed_phase timer("step 3. build unitig-to-color map");

    bits::bit_vector::builder u2c_builder;
    u2c_builder.resize(m_num_unitigs, 0);

    std::filesystem::path unitigs_fn = m_config.base_filename;
    unitigs_fn.replace_extension("unitigs.fa");
    std::ifstream in(unitigs_fn.c_str());
    if (!in.is_open()) throw std::runtime_error("cannot open unitigs file");

    uint64_t prev = static_cast<uint64_t>(-1);
    uint64_t count = 0;
    const std::string target = "color_set_id=";
    const uint64_t target_length = target.length();
    std::string line;
    for (uint64_t i = 0; i != m_num_unitigs; ++i) {
        std::getline(in, line);  // read header
        size_t pos = line.find(target);
        assert(pos != std::string::npos);
        char const* p = line.c_str() + pos + target_length;
        uint64_t color_set_id = std::strtoull(p, nullptr, 10);
        if (color_set_id != prev) {
            count += 1;
            if (i > 0) u2c_builder.set(i - 1, 1);
        }
        prev = color_set_id;
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // skip unitig sequence
    }
    assert(count == m_num_color_sets);
    (void)count;
    in.close();

    bits::bit_vector u2c;
    bits::rank9 u2c_rank1_index;
    u2c_builder.set(m_num_unitigs - 1, 1);
    u2c_builder.build(u2c);
    u2c_rank1_index.build(u2c);
    assert(u2c.num_bits() == m_num_unitigs);
    assert(u2c_rank1_index.num_ones() == m_num_color_sets);

    m_saver.visit(u2c);
    m_saver.visit(u2c_rank1_index);

    std::cout << "m_u2c.num_bits() " << u2c.num_bits() << std::endl;
    std::cout << "m_u2c_rank1_index.num_ones() " << u2c_rank1_index.num_ones() << std::endl;
}

template <typename ColorSets>
void index<ColorSets>::loader::load_kmer_dictionary() {
    util::timed_phase timer("step 4. build SSHash");

    std::filesystem::path unitigs_fn = m_config.base_filename;
    unitigs_fn.replace_extension("unitigs.fa");

    sshash::build_configuration sshash_config;
    sshash_config.k = m_k;
    sshash_config.m = m_config.m;
    sshash_config.verbose = m_config.verbose;
    sshash_config.tmp_dirname = m_config.tmp_dirname;
    sshash_config.num_threads = m_config.num_threads;
    sshash_config.print();

    sshash::dictionary_type k2u;
    k2u.build(unitigs_fn, sshash_config);
    m_saver.visit(k2u);
}

template <typename ColorSets>
void index<ColorSets>::loader::load_filenames() {
    util::timed_phase timer("step 5. permute and write filenames");

    std::filesystem::path filenames_fn = m_config.base_filename;
    filenames_fn.replace_extension("filenames.txt");

    filenames filenames;
    filenames.build_from_file(filenames_fn);
    m_saver.visit(filenames);
}

template <>
inline void index<hybrid>::loader::load_color_sets() {
    util::timed_phase timer("step 2. encode color sets");

    std::filesystem::path color_sets_fn = m_config.base_filename;
    color_sets_fn.replace_extension("color_sets.fa");

    hybrid::builder color_sets_builder(m_num_colors, m_saver, m_config);

    std::ifstream in(color_sets_fn);
    if (!in.is_open()) throw std::runtime_error("cannot open color sets file");

    std::string line;
    std::vector<uint32_t> v;
    const std::string target = "size=";
    const uint64_t target_length = target.length();
    for (uint64_t i = 0; i != m_num_color_sets; ++i) {
        std::getline(in, line);
        size_t size_pos = line.find(target);
        assert(size_pos != std::string::npos);
        char const* p = line.c_str() + size_pos + target_length;
        char* endptr = nullptr;
        uint64_t color_set_size = std::strtoul(p, &endptr, 10);
        assert(color_set_size > 0);
        p = endptr;
        v.clear();
        v.reserve(color_set_size);
        for (uint64_t i = 0; i != color_set_size; ++i) {
            v.push_back(static_cast<uint32_t>(std::strtoul(p, &endptr, 10)));
            p = endptr;
        }
        color_sets_builder.encode(v);
    }

    in.close();
    color_sets_builder.build();
}

}  // namespace fulgor
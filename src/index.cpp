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

    std::cout << "total index size: " << total_bits / 8 << " [B] -- "
              << essentials::convert(total_bits / 8, essentials::GB) << " [GB]" << '\n';
    std::cout << "SPACE BREAKDOWN:\n";
    std::cout << "  dBG (SSHash): " << k2u.num_bits() / 8 << " bytes / "
              << essentials::convert(k2u.num_bits() / 8, essentials::GB) << " GB ("
              << (k2u.num_bits() * 100.0) / total_bits << "%)\n";
    std::cout << "  Color sets: " << color_sets.num_bits() / 8 << " bytes / "
              << essentials::convert(color_sets.num_bits() / 8, essentials::GB) << " GB ("
              << (color_sets.num_bits() * 100.0) / total_bits << "%)\n";
    uint64_t other_bits =
        (u2c.num_bytes() + u2c_rank1_index.num_bytes()) * 8 + filenames.num_bits();
    std::cout << "  Other: " << other_bits / 8 << " bytes / "
              << essentials::convert(other_bits / 8, essentials::GB) << " GB ("
              << (other_bits * 100.0) / total_bits << "%)\n";
    std::cout << "    Map from unitig_id to color_set_id: "
              << u2c.num_bytes() + u2c_rank1_index.num_bytes() << " bytes / "
              << essentials::convert(u2c.num_bytes() + u2c_rank1_index.num_bytes(), essentials::GB)
              << " GB ("
              << ((u2c.num_bytes() + u2c_rank1_index.num_bytes()) * 8 * 100.0) / total_bits
              << "%)\n";
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
    essentials::logger("writing '" + build_config.file_base_name + ".metadata.txt'...");
    std::ofstream metadata_file(build_config.file_base_name + ".metadata.txt");
    if (!metadata_file.is_open()) throw std::runtime_error("cannot open output file");
    metadata_file << "k=" << k() << '\n';
    metadata_file << "num_kmers=" << num_kmers() << '\n';
    metadata_file << "num_colors=" << num_colors() << '\n';
    metadata_file << "num_unitigs=" << num_unitigs() << '\n';
    metadata_file << "num_color_sets=" << num_color_sets() << '\n';
    metadata_file.close();

    /* filenames file */
    essentials::logger("writing '" + build_config.file_base_name + ".filenames.txt'...");
    std::ofstream filenames_file(build_config.file_base_name + ".filenames.txt");
    if (!filenames_file.is_open()) throw std::runtime_error("cannot open output file");
    for (uint64_t i = 0; i != num_colors(); ++i) filenames_file << filename(i) << '\n';
    filenames_file.close();

    /* unitigs file */
    essentials::logger("writing '" + build_config.file_base_name + ".unitigs.fa'...");
    std::ofstream unitigs_file(build_config.file_base_name + ".unitigs.fa");
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
    essentials::logger("writing '" + build_config.file_base_name + ".color_sets.txt'...");
    std::ofstream color_sets_file(build_config.file_base_name + ".color_sets.txt");
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
void index<ColorSets>::load(build_configuration const& build_config)  //
{
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;

    uint64_t k = 0;
    uint64_t num_kmers = 0;
    uint64_t num_colors = 0;
    uint64_t num_unitigs = 0;
    uint64_t num_color_sets = 0;

    {
        essentials::logger("step 1. reading metadata...");

        std::string metadata_fn = build_config.file_base_name + ".metadata.txt";
        std::ifstream in(metadata_fn.c_str());
        if (!in.is_open()) throw std::runtime_error("cannot open metadata file");

        std::string line;
        while (std::getline(in, line)) {
            size_t delimiter_pos = line.find('=');
            assert(delimiter_pos != std::string::npos);
            std::string_view key(line.data(), delimiter_pos);
            char const* value_ptr = line.c_str() + delimiter_pos + 1;
            if (key == "k") {
                k = static_cast<uint32_t>(std::strtoul(value_ptr, nullptr, 10));
            } else if (key == "num_kmers") {
                num_kmers = std::strtoull(value_ptr, nullptr, 10);
            } else if (key == "num_colors") {
                num_colors = static_cast<uint32_t>(std::strtoul(value_ptr, nullptr, 10));
            } else if (key == "num_unitigs") {
                num_unitigs = std::strtoull(value_ptr, nullptr, 10);
            } else if (key == "num_color_sets") {
                num_color_sets = std::strtoull(value_ptr, nullptr, 10);
            }
        }

        in.close();

        if (build_config.verbose) {
            std::cout << "k=" << k << ", num_kmers=" << num_kmers << ", num_colors=" << num_colors
                      << ", num_unitigs=" << num_unitigs << ", num_color_sets=" << num_color_sets
                      << std::endl;
        }

        assert(num_unitigs > 0);
        assert(num_unitigs < (uint64_t(1) << 32));

        essentials::logger("DONE");
    }

    {
        essentials::logger("step 2. building unitig-to-color map...");
        timer.start();

        bits::bit_vector::builder u2c_builder;
        u2c_builder.resize(num_unitigs, 0);

        std::string unitigs_fn = build_config.file_base_name + ".unitigs.fa";
        std::ifstream in(unitigs_fn.c_str());
        if (!in.is_open()) throw std::runtime_error("cannot open unitigs file");

        uint64_t prev = uint64_t(-1);
        uint64_t count = 0;
        const std::string target = "color_set_id=";
        const uint64_t target_length = target.length();
        std::string line;
        for (uint64_t i = 0; i != num_unitigs; ++i) {
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
        assert(count == num_color_sets);
        (void)count;
        in.close();

        u2c_builder.set(num_unitigs - 1, 1);
        u2c_builder.build(m_u2c);
        m_u2c_rank1_index.build(m_u2c);
        assert(m_u2c.num_bits() == num_unitigs);
        assert(m_u2c_rank1_index.num_ones() == num_color_sets);

        std::cout << "m_u2c.num_bits() " << m_u2c.num_bits() << std::endl;
        std::cout << "m_u2c_rank1_index.num_ones() " << m_u2c_rank1_index.num_ones() << std::endl;

        timer.stop();
        std::cout << "** building unitig-to-color map took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        timer.reset();
    }

    {
        essentials::logger("step 3. encoding color sets...");
        timer.start();

        typename ColorSets::builder color_sets_builder(num_colors);
        const uint64_t num_bits = essentials::GiB * 8 * 8;
        color_sets_builder.reserve_num_bits(num_bits);

        std::string color_sets_fn = build_config.file_base_name + ".color_sets.txt";
        std::ifstream in(color_sets_fn);
        if (!in.is_open()) throw std::runtime_error("cannot open color sets file");

        std::string line;
        std::vector<uint32_t> v;
        const std::string target = "size=";
        const uint64_t target_length = target.length();
        for (uint64_t i = 0; i != num_color_sets; ++i) {
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
            color_sets_builder.encode_color_set(v.data(), v.size());
        }

        in.close();
        color_sets_builder.build(m_color_sets);

        timer.stop();
        std::cout << "** encoding color sets took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        timer.reset();
    }

    {
        essentials::logger("step 3. building SSHash...");
        timer.start();
        std::string unitigs_fn = build_config.file_base_name + ".unitigs.fa";
        sshash::build_configuration sshash_config;
        sshash_config.k = k;
        sshash_config.m = build_config.m;
        sshash_config.canonical = true;
        sshash_config.verbose = build_config.verbose;
        sshash_config.tmp_dirname = build_config.tmp_dirname;
        sshash_config.num_threads = build_config.num_threads;
        sshash_config.print();
        m_k2u.build(unitigs_fn, sshash_config);
        timer.stop();
        std::cout << "** building SSHash took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        timer.reset();
    }

    {
        essentials::logger("step 4. reading filenames...");
        timer.start();
        std::string filenames_fn = build_config.file_base_name + ".filenames.txt";
        std::ifstream in(filenames_fn.c_str());
        if (!in.is_open()) throw std::runtime_error("cannot open filenames file");
        std::vector<std::string> filenames;
        filenames.reserve(num_colors);
        std::string filename;
        for (uint64_t i = 0; i != num_colors; ++i) {
            in >> filename;
            filenames.push_back(filename);
        }
        in.close();
        m_filenames.build(filenames);
        timer.stop();
        std::cout << "** building filenames took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        timer.reset();
    }

    if (build_config.verbose) print_stats();
}

}  // namespace fulgor
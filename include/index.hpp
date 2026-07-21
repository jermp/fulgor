#pragma once

#include "external/sshash/include/dictionary_types.hpp"
#include "external/sshash/external/pthash/external/bits/include/integer_codes.hpp"
#include "external/sshash/external/pthash/external/bits/include/bit_vector.hpp"
#include "external/sshash/external/pthash/external/bits/include/rank9.hpp"

#include "filenames.hpp"
#include "util.hpp"

namespace fulgor {

using kmer_type = sshash::default_kmer_t;
using sshash_type = sshash::dictionary_type;

template <typename ColorSets>
struct index {
    typedef ColorSets color_sets_type;

    // struct builder;
    struct hybrid_builder;
    struct meta_builder;
    struct differential_builder;
    struct meta_differential_builder;

    index()
        : m_vnum(constants::current_version_number::major,  //
                 constants::current_version_number::minor,  //
                 constants::current_version_number::patch)  //
    {}

    typename color_sets_type::iterator_type color_set(uint64_t color_set_id) const {
        assert(color_set_id < num_color_sets());
        return m_color_sets.color_set(color_set_id);
    }

    /* from unitig_id to color_set_id */
    uint64_t u2c(uint64_t unitig_id) const { return m_u2c_rank1_index.rank1(m_u2c, unitig_id); }

    void fetch_color_set_ids(std::string const& sequence,
                             std::vector<uint32_t>& color_set_ids) const;
    void pseudoalign_full_intersection(std::vector<uint32_t>& color_set_ids,  //
                                       std::vector<uint32_t>& results,
                                       std::vector<uint32_t>& tmp) const;  //
    void pseudoalign_threshold_union(std::string const& sequence,          //
                                     std::vector<uint32_t>& results,       //
                                     const double threshold) const;        //

    void kmer_conservation(std::string const& sequence,                                           //
                           std::vector<kmer_conservation_triple>& kmer_conservation_info) const;  //

    void kmer_matches(std::string const& sequence,                            //
                      bits::bit_vector::builder& positive_kmers_in_sequence,  //
                      std::vector<count_type>& counts) const;                 //

    std::string_view filename(uint64_t color) const {
        assert(color < num_colors());
        return m_filenames[color];
    }

    void print_stats() const;
    void dump(build_configuration const& build_config) const;
    void load(build_configuration const& build_config);

    uint64_t k() const { return m_k2u.k(); }
    uint64_t num_kmers() const { return m_k2u.num_kmers(); }
    uint64_t num_colors() const { return m_color_sets.num_colors(); }
    uint64_t num_unitigs() const { return m_k2u.num_strings(); }
    uint64_t num_color_sets() const { return m_color_sets.num_color_sets(); }

    sshash_type const& get_k2u() const { return m_k2u; }
    bits::bit_vector const& get_u2c() const { return m_u2c; }
    bits::rank9 const& get_u2c_rank1_index() const { return m_u2c_rank1_index; }
    ColorSets const& get_color_sets() const { return m_color_sets; }
    filenames const& get_filenames() const { return m_filenames; }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    uint64_t num_bits() const {
        return m_k2u.num_bits() +
               (sizeof(m_vnum) + m_u2c.num_bytes() + m_u2c_rank1_index.num_bytes()) * 8 +
               m_color_sets.num_bits() + m_filenames.num_bits();
    }

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_color_sets);
        visitor.visit(t.m_u2c);
        visitor.visit(t.m_u2c_rank1_index);
        visitor.visit(t.m_k2u);
        visitor.visit(t.m_filenames);
        visitor.visit(t.m_vnum);
        util::check_version_number(t.m_vnum);
    }

    ColorSets m_color_sets;
    bits::bit_vector m_u2c;
    bits::rank9 m_u2c_rank1_index;
    sshash_type m_k2u;
    filenames m_filenames;
    essentials::version_number m_vnum;
};

/*
template <typename ColorSets>
struct index<ColorSets>::builder {
    explicit builder(const build_configuration& build_config)
        : m_build_config(build_config), m_saver(build_config.output_filename) {}

    void build_cdbg();
    void build_color_sets();

    void build_u2c() {
        essentials::logger("step 2. building unitig-to-color map and encoding color sets...");
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        timer.start();

        std::ifstream cs_file(cdbg_build_config.cs_filename(), std::ios::binary);
        m_saver.append(cs_file, 12);      // write num_cols, sp_thresh, vd_thresh
        cs_file.seekg(8, std::ios::cur);  // skip num_color_sets
        m_saver.append(cs_file);
        cs_file.close();

        assert(cdbg_builder.num_unitigs() > 0);
        assert(cdbg_builder.num_unitigs() <= UINT32_MAX);

        std::cout << "num_unitigs " << cdbg_builder.num_unitigs() << std::endl;
        std::cout << "num_distinct_color_sets " << cdbg_builder.num_color_sets() << std::endl;

        timer.stop();
        std::cout << "** encoding color sets took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        timer.reset();

        timer.start();

        bits::bit_vector u2c;
        bits::rank9 u2c_rank1_index;
        essentials::load(u2c, cdbg_build_config.u2c_filename().c_str());
        u2c_rank1_index.build(u2c);
        m_saver.visit(u2c);
        m_saver.visit(u2c_rank1_index);

        assert(u2c.num_bits() == cdbg_builder.num_unitigs());
        assert(u2c_rank1_index.num_ones() == cdbg_builder.num_color_sets());

        std::cout << "m_u2c.num_bits() " << u2c.num_bits() << std::endl;
        std::cout << "m_u2c_rank1_index.num_ones() " << u2c_rank1_index.num_ones() << std::endl;

        timer.stop();
        std::cout << "** building unitig-to-color map took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        timer.reset();
    }

    void build_kmer_index(std::filesystem::path fasta_filename) {
        essentials::logger("step 3. building SSHash...");
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        timer.start();

        sshash::build_configuration sshash_config;
        sshash_config.k = m_build_config.k;
        sshash_config.m = m_build_config.m;
        sshash_config.canonical = true;
        sshash_config.verbose = m_build_config.verbose;
        sshash_config.tmp_dirname = m_build_config.tmp_dirname;
        sshash_config.num_threads = m_build_config.num_threads;
        sshash_config.print();

        sshash::dictionary_type k2u;
        k2u.build(fasta_filename, sshash_config);
        m_saver.visit(k2u);

        try {  // remove unitig file
            std::remove(fasta_filename.c_str());
        } catch (std::exception const& e) {
            std::cerr << e.what() << std::endl;
        }

        timer.stop();
        std::cout << "** building SSHash took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
        timer.reset();
    }

    void build_filenames() {
        essentials::logger("step 4. writing filenames...");
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
        timer.start();

        filenames filenames;
        filenames.build_from_file(m_build_config.filenames_list);
        m_saver.visit(filenames);

        timer.stop();
        std::cout << "** writing filenames took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
    }

    void build() {
        build_color_sets();
        build_u2c();
        build_kmer_index();
        build_filenames();

        m_saver.write(constants::current_version_number::major);
        m_saver.write(constants::current_version_number::minor);
        m_saver.write(constants::current_version_number::patch);
    }

private:
    build_configuration m_build_config;
    util::external_saver m_saver;
};
*/

}  // namespace fulgor

/*
#include <concepts>
#include <filesystem>

template <typename T>
concept IndexBuildingStrategy = requires(T strategy, std::filesystem::path path) {
    { strategy.build_cdbg() }        -> std::same_as<void>;
    { strategy.build_color_sets() }  -> std::same_as<void>;
    { strategy.build_u2c() }         -> std::same_as<void>;
    { strategy.build_filenames() }   -> std::same_as<void>;
    { strategy.build_kmer_index(path) } -> std::same_as<void>;

    // The strategy provides access to its internal saver
    { strategy.saver() } -> std::derived_from<util::external_saver>;
};

 *

template <IndexBuildingStrategy Strategy>
class index_builder {
public:
    // Accept the specific building strategy via aggregate construction or forward values
    explicit index_builder(build_configuration config, Strategy strategy)
        : m_build_config(std::move(config)), m_strategy(std::move(strategy)) {}

    void build(std::filesystem::path fasta_filename) {
        // Safe, rigid execution order
        m_strategy.build_cdbg();
        m_strategy.build_color_sets();
        m_strategy.build_u2c();
        m_strategy.build_kmer_index(fasta_filename);
        m_strategy.build_filenames();

        // Metadata footer handled consistently
        auto& saver = m_strategy.saver();
        saver.write(constants::current_version_number::major);
        saver.write(constants::current_version_number::minor);
        saver.write(constants::current_version_number::patch);
    }

private:
    build_configuration m_build_config;
    Strategy m_strategy;
};

*

struct modern_color_sets_strategy {
    explicit modern_color_sets_strategy(const build_configuration& config)
        : m_saver(config.output_filename) {}

    util::external_saver& saver() { return m_saver; }

    void build_cdbg() { }
    void build_color_sets() {}
    void build_u2c() { }
    void build_filenames() {}
    void build_kmer_index(std::filesystem::path fasta_filename) {}

private:
util::external_saver m_saver;
};

 *

build_configuration config = load_config();

// The compiler uses CTAD (Class Template Argument Deduction) to automatically
// deduce the strategy type: index_builder<modern_color_sets_strategy>
index_builder builder(config, modern_color_sets_strategy{config});

// Executes safely, efficiently, and in order
builder.build("data.fasta");
 */
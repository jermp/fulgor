#pragma once

#include <concepts>
#include <utility>

#include "external/sshash/include/dictionary_types.hpp"
#include "external/sshash/external/pthash/external/bits/include/integer_codes.hpp"
#include "external/sshash/external/pthash/external/bits/include/bit_vector.hpp"
#include "external/sshash/external/pthash/external/bits/include/rank9.hpp"

#include "filenames.hpp"
#include "util.hpp"

namespace fulgor {

using kmer_type = sshash::default_kmer_t;
using sshash_type = sshash::dictionary_type;

template <typename T>
concept IndexBuildingStrategy =
    std::constructible_from<T, const build_configuration&> && requires(T strategy) {
        { strategy.build_cdbg() } -> std::same_as<void>;
        { strategy.build_color_sets() } -> std::same_as<void>;
        { strategy.build_u2c() } -> std::same_as<void>;
        { strategy.build_filenames() } -> std::same_as<void>;
        { strategy.build_kmer_dictionary() } -> std::same_as<void>;

        { strategy.saver() } -> std::convertible_to<util::external_saver&>;
    };

struct unsupported_builder {};
template <typename ColorSets>
struct builder_strategy {
    using type = unsupported_builder;
};

template <typename ColorSets>
struct index {
    typedef ColorSets color_sets_type;

    using builder_strategy_t = builder_strategy<ColorSets>::type;
    struct builder {
        explicit builder(const build_configuration& config)
            requires IndexBuildingStrategy<builder_strategy_t>
            : m_strategy(config) {}

        void build() {
            m_strategy.build_cdbg();
            m_strategy.build_color_sets();
            m_strategy.build_u2c();
            m_strategy.build_kmer_dictionary();
            m_strategy.build_filenames();

            auto&& saver = m_strategy.saver();
            saver.write(constants::current_version_number::major);
            saver.write(constants::current_version_number::minor);
            saver.write(constants::current_version_number::patch);
        }

    private:
        builder_strategy_t m_strategy;
    };

    struct loader {
        explicit loader(build_configuration config) : m_config(std::move(config)) {
            m_saver.set_stream(
                std::ofstream(m_config.output_filename, std::ios::binary | std::ios::trunc));
        }

        void load() {
            load_metadata();
            load_color_sets();
            load_u2c();
            load_kmer_dictionary();
            load_filenames();

            m_saver.write(constants::current_version_number::major);
            m_saver.write(constants::current_version_number::minor);
            m_saver.write(constants::current_version_number::patch);
        }

    private:
        build_configuration m_config;
        util::external_saver m_saver;

        uint64_t m_k;
        uint64_t m_num_kmers;
        uint64_t m_num_colors;
        uint64_t m_num_unitigs;
        uint64_t m_num_color_sets;

        void load_metadata();
        void load_color_sets();
        void load_u2c();
        void load_kmer_dictionary();
        void load_filenames();
    };

    // struct hybrid_builder;
    // struct meta_builder;
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

    std::string_view filename(const uint64_t color) const {
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

}  // namespace fulgor
#pragma once

#include "external/sshash/include/dictionary.hpp"
#include "external/sshash/external/pthash/external/bits/include/integer_codes.hpp"
#include "external/sshash/external/pthash/external/bits/include/bit_vector.hpp"
#include "external/sshash/external/pthash/external/bits/include/rank9.hpp"

#include "filenames.hpp"
#include "util.hpp"

namespace fulgor {

using kmer_type = sshash::default_kmer_t;
using sshash_type = sshash::dictionary<kmer_type>;

template <typename ColorSets>
struct index {
    typedef ColorSets color_sets_type;

    struct builder;
    struct meta_builder;
    struct differential_builder;
    struct meta_differential_builder;

    index()
        : m_vnum(constants::current_version_number::x,  //
                 constants::current_version_number::y,  //
                 constants::current_version_number::z)  //
    {}

    typename color_sets_type::iterator_type color_set(uint64_t color_set_id) const {
        assert(color_set_id < num_color_sets());
        return m_color_sets.color_set(color_set_id);
    }

    /* from unitig_id to color_set_id */
    uint64_t u2c(uint64_t unitig_id) const { return m_u2c_rank1_index.rank1(m_u2c, unitig_id); }

    void pseudoalign_full_intersection(std::string const& sequence,            //
                                       std::vector<uint32_t>& results) const;  //

    void pseudoalign_threshold_union(std::string const& sequence,     //
                                     std::vector<uint32_t>& results,  //
                                     const double threshold) const;   //

    void kmer_conservation(std::string const& sequence,                                           //
                           std::vector<kmer_conservation_triple>& kmer_conservation_info) const;  //

    std::string_view filename(uint64_t color) const {
        assert(color < num_colors());
        return m_filenames[color];
    }

    void print_stats() const;
    void dump(std::string const& basename) const;

    uint64_t k() const { return m_k2u.k(); }
    uint64_t num_colors() const { return m_color_sets.num_colors(); }
    uint64_t num_unitigs() const { return m_k2u.num_contigs(); }
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
        visitor.visit(t.m_vnum);
        util::check_version_number(t.m_vnum);
        visitor.visit(t.m_k2u);
        visitor.visit(t.m_u2c);
        visitor.visit(t.m_u2c_rank1_index);
        visitor.visit(t.m_color_sets);
        visitor.visit(t.m_filenames);
    }

    essentials::version_number m_vnum;
    sshash_type m_k2u;
    bits::bit_vector m_u2c;
    bits::rank9 m_u2c_rank1_index;
    ColorSets m_color_sets;
    filenames m_filenames;
};

}  // namespace fulgor

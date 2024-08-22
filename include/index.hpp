#pragma once

#include "external/sshash/include/dictionary.hpp"
#include "ranked_bit_vector.hpp"
#include "filenames.hpp"
#include "util.hpp"

namespace fulgor {

template <typename ColorSets>
struct index {
    typedef ColorSets color_sets_type;

    struct builder;
    struct meta_builder;
    struct differential_builder;
    struct meta_differential_builder;

    typename color_sets_type::iterator_type color_set(uint64_t color_set_id) const {
        assert(color_set_id < num_color_sets());
        return m_color_sets.color_set(color_set_id);
    }

    /* from unitig_id to color_set_id */
    uint64_t u2c(uint64_t unitig_id) const { return m_u2c.rank(unitig_id); }

    void pseudoalign_full_intersection(std::string const& sequence,
                                       std::vector<uint32_t>& results) const;
    void pseudoalign_threshold_union(std::string const& sequence, std::vector<uint32_t>& results,
                                     const double threshold) const;

    void intersect_unitigs(std::vector<uint64_t>& unitig_ids,
                           std::vector<uint32_t>& color_set) const;

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

    sshash::dictionary const& get_k2u() const { return m_k2u; }
    ranked_bit_vector const& get_u2c() const { return m_u2c; }
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
        return m_k2u.num_bits() + m_u2c.bytes() * 8 + m_color_sets.num_bits() +
               m_filenames.num_bits();
    }

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_k2u);
        visitor.visit(t.m_u2c);
        visitor.visit(t.m_color_sets);
        visitor.visit(t.m_filenames);
    }

    sshash::dictionary m_k2u;  // map: kmer to unitig-id
    ranked_bit_vector m_u2c;   // map: unitig-id to color-class-id
    ColorSets m_color_sets;
    filenames m_filenames;
};

}  // namespace fulgor

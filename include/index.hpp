#pragma once

#include "external/sshash/include/dictionary.hpp"
#include "ranked_bit_vector.hpp"
#include "filenames.hpp"
#include "util.hpp"

namespace fulgor {

template <typename ColorClasses>
struct index {
    typedef ColorClasses color_classes_type;

    struct builder;
    struct meta_builder;
    struct differential_builder;
    struct meta_differential_builder;

    typename color_classes_type::iterator_type color_set(uint64_t color_set_id) const {
        assert(color_set_id < num_color_sets());
        return m_ccs.color_set(color_set_id);
    }

    /* from unitig_id to color_set_id */
    uint64_t u2c(uint64_t unitig_id) const { return m_u2c.rank(unitig_id); }

    void pseudoalign_full_intersection(std::string const& sequence,
                                       std::vector<uint32_t>& results) const;
    void pseudoalign_threshold_union(std::string const& sequence, std::vector<uint32_t>& results,
                                     const double threshold) const;

    void intersect_unitigs(std::vector<uint64_t>& unitig_ids, std::vector<uint32_t>& color_set) const;

    std::string_view filename(uint64_t doc_id) const {
        assert(doc_id < num_docs());
        return m_filenames.filename(doc_id);
    }

    void print_stats() const;
    void dump_colors(std::ofstream& os) const;

    uint64_t k() const { return m_k2u.k(); }
    uint64_t num_docs() const { return m_ccs.num_docs(); }
    uint64_t num_color_sets() const { return m_ccs.num_color_sets(); }

    sshash::dictionary const& get_k2u() const { return m_k2u; }
    ranked_bit_vector const& get_u2c() const { return m_u2c; }
    ColorClasses const& get_color_sets() const { return m_ccs; }
    filenames const& get_filenames() const { return m_filenames; }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_k2u);
        visitor.visit(m_u2c);
        visitor.visit(m_ccs);
        visitor.visit(m_filenames);
    }

    uint64_t num_bits() const {
        return m_k2u.num_bits() + m_u2c.bytes() * 8 + m_ccs.num_bits() + m_filenames.num_bits();
    }

private:
    sshash::dictionary m_k2u;  // map: kmer to unitig-id
    ranked_bit_vector m_u2c;   // map: unitig-id to color-class-id
    ColorClasses m_ccs;
    filenames m_filenames;
};

}  // namespace fulgor

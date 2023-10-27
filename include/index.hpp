#pragma once

#include "../external/sshash/include/dictionary.hpp"
#include "ranked_bit_vector.hpp"
#include "filenames.hpp"
#include "util.hpp"

namespace fulgor {

template <typename ColorClasses>
struct index {
    typedef ColorClasses color_classes_type;

    void build_from(ccdbg_builder const& cb);

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_k2u);
        visitor.visit(m_u2c);
        visitor.visit(m_ccs);
        visitor.visit(m_filenames);
    }

    uint64_t num_docs() const { return m_ccs.num_docs(); }

    ColorClasses const& color_classes() const { return m_ccs; }

    /* from unitig_id to color_class_id */
    uint64_t u2c(uint64_t unitig_id) const { return m_u2c.rank(unitig_id); }

    void pseudoalign_full_intersection(std::string const& sequence,
                                       std::vector<uint32_t>& results) const;
    void pseudoalign_threshold_union(std::string const& sequence, std::vector<uint32_t>& results,
                                     const double threshold) const;

    // instead of mapping a read to the actual references, just return 
    // the list of color ids that must be intersected to obtain that
    // list of references.
    void pseudoalign_full_intersection_color_ids(std::string const& sequence,
                                                 std::vector<uint32_t>& results) const;

    void intersect_unitigs_to_ids(std::vector<uint32_t>& unitig_ids, std::vector<uint32_t>& colors) const;

    void intersect_unitigs(std::vector<uint32_t>& unitig_ids, std::vector<uint32_t>& colors) const;

    void intersect_color_ids(const std::vector<uint32_t>& color_ids, const size_t skip, std::vector<uint32_t>& colors) const;

    uint64_t num_bits() const {
        return m_k2u.num_bits() + m_u2c.bytes() * 8 + m_ccs.num_bits() + m_filenames.num_bits();
    }

    std::string_view filename(uint64_t doc_id) const {
        assert(doc_id < num_docs());
        return m_filenames.filename(doc_id);
    }

    void print_stats() const;

    size_t k() const { return m_k2u.k(); }
    sshash::dictionary const* get_dict() const { return &m_k2u; }
    pthash::bit_vector const& contigs() const { return m_k2u.strings(); }

private:
    sshash::dictionary m_k2u;
    ranked_bit_vector m_u2c;
    ColorClasses m_ccs;
    filenames m_filenames;
};

}  // namespace fulgor

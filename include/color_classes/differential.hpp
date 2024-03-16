#pragma once

namespace fulgor {

template <typename ColorClasses>
struct differential {
    static const bool meta_colored = false;
    static const bool differential_colored = true;

    struct builder {
        builder() {
            // Giulio: perhaps here some
            // m_bvb.reserve() ?
        }

        void encode_reference(std::vector<uint32_t> const& reference) {
            util::write_delta(m_bvb, reference.size());

            // Giulio: take the gaps and then encode with delta, like:
            //   gap = color_id - (prev + 1)
            //   write_delta(gap)
            // and the first integer written with no gap, of course.
            for (auto color_id : reference) { util::write_delta(m_bvb, color_id); }
        }

        void encode_list(std::vector<uint32_t> const& reference,
                         typename ColorClasses::forward_iterator it) {
            // Giulio:
            // same comment as above: take gaps + delta
            // you need to keep track of the last written integer

            uint64_t i = 0, j = 0;
            while (i < it.size() && j < reference.size()) {
                if (*it == reference[j]) {
                    i += 1;
                    ++it;
                    j += 1;
                } else if (*it < reference[j]) {
                    util::write_delta(m_bvb, *it);
                    i += 1;
                    ++it;
                } else {
                    util::write_delta(m_bvb, reference[j]);
                    j += 1;
                }
            }
            while (i < it.size()) {
                util::write_delta(m_bvb, *it);
                i += 1;
                ++it;
            }

            while (j < reference.size()) {
                util::write_delta(m_bvb, reference[j]);
                j += 1;
            }
        }

    private:
        //
        // Giulio: these two would be needed to delimit
        // where each list begins in the bits of m_bvb.
        //
        // Additionally, we would need a bitvector + rank
        // to compute for each differential list its corresponding
        // reference. (See the notes I sent you.)
        std::vector<uint64_t> m_reference_offsets;
        std::vector<uint64_t> m_list_offsets;

        bit_vector_builder m_bvb;
    };

    struct forward_iterator {};

    typedef forward_iterator iterator_type;
};
}  // namespace fulgor
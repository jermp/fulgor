#pragma once

namespace fulgor {

template <typename ColorClasses>
struct differential {
    static const bool meta_colored = false;
    static const bool differential_colored = true;

    struct builder {
        builder() {
            m_list_offsets.push_back(0);
            m_reference_offsets.push_back(0);
        }

        void encode_reference(std::vector<uint32_t> const& reference) {
            uint64_t size = reference.size();
            util::write_delta(m_bvb, size);

            if (size == 0) return;

            uint32_t prev_val = reference[0];
            util::write_delta(m_bvb, prev_val);

            for (uint64_t i = 1; i < size; ++i) {
                uint32_t val = reference[i];
                assert(val >= prev_val + 1);
                util::write_delta(m_bvb, val - (prev_val + 1));
                prev_val = val;
            }
            m_reference_offsets.push_back(m_bvb.num_bits());
        }

        void encode_list(std::vector<uint32_t> const& reference,
                         typename ColorClasses::forward_iterator it) {
            std::vector<uint32_t> edit_list;
            uint64_t ref_size = reference.size();
            uint64_t it_size = it.size();
            edit_list.reserve(ref_size + it_size);

            uint64_t i = 0, j = 0;
            while (i < it_size && j < ref_size) {
                if (*it == reference[j]) {
                    i += 1;
                    ++it;
                    j += 1;
                } else if (*it < reference[j]) {
                    edit_list.push_back(*it);
                    i += 1;
                    ++it;
                } else {
                    edit_list.push_back(reference[j]);
                    j += 1;
                }
            }
            while (i < it_size) {
                edit_list.push_back(*it);
                i += 1;
                ++it;
            }
            while (j < ref_size) {
                edit_list.push_back(reference[j]);
                j += 1;
            }

            uint64_t size = edit_list.size();
            util::write_delta(m_bvb, size);

            if (size == 0) return;

            uint32_t prev_val = edit_list[0];
            util::write_delta(m_bvb, prev_val);

            for (uint64_t i = 1; i < size; ++i) {
                uint32_t val = edit_list[i];
                assert(val >= prev_val + 1);
                util::write_delta(m_bvb, val - (prev_val + 1));
                prev_val = val;
            }

            m_list_offsets.push_back(m_bvb.num_bits());
        }

    private:
        std::vector<uint64_t> m_reference_offsets;
        std::vector<uint64_t> m_list_offsets;

        bit_vector_builder m_bvb;
    };

    struct forward_iterator {};

    typedef forward_iterator iterator_type;
};
}  // namespace fulgor
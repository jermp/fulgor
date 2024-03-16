#pragma once

namespace fulgor {

template <typename ColorClasses>
struct differential {
    static const bool meta_colored = false;
    static const bool differential_colored = true;

    struct builder {

        void encode_reference(std::vector<uint32_t> const& reference){
            util::write_delta(m_bvb, reference.size());
            for(auto color_id : reference){
                util::write_delta(m_bvb, color_id);
            }
        }

        void encode_list(std::vector<uint32_t> const& reference, typename ColorClasses::forward_iterator it){
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
                    util::write_delta(m_bvb, reference[j]);  // note: should be -R[j]
                    j += 1;
                }
            }
            while (i < it.size()) {
                util::write_delta(m_bvb, *it);
                i += 1;
                ++it;
            }

            while (j < reference.size()) {
                util::write_delta(m_bvb, reference[j]);  // note: should be -R[j]
                j += 1;
            }
        }


    private:

        bit_vector_builder m_bvb;
    };

    struct forward_iterator {

    };

    typedef forward_iterator iterator_type;

};
}
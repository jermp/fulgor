#pragma once

namespace fulgor {

struct hybrid {
    static const index_t type = index_t::HYBRID;

    struct builder {
        builder() : m_num_lists(0) {}
        builder(uint64_t num_colors) { init(num_colors); }

        void init(uint64_t num_colors) {
            m_num_colors = num_colors;

            /* if list contains < sparse_set_threshold_size ints, code it with gaps+delta */
            m_sparse_set_threshold_size = 0.25 * m_num_colors;

            /* if list contains > very_dense_set_threshold_size ints, code it as a complementary set
               with gaps+delta */
            m_very_dense_set_threshold_size = 0.75 * m_num_colors;
            /* otherwise: code it as a bitmap of m_num_colors bits */

            std::cout << "m_num_colors: " << m_num_colors << std::endl;
            std::cout << "m_sparse_set_threshold_size " << m_sparse_set_threshold_size << std::endl;
            std::cout << "m_very_dense_set_threshold_size " << m_very_dense_set_threshold_size
                      << std::endl;

            m_bvb.reserve(320 * essentials::GB);
            m_offsets.push_back(0);

            m_num_lists = 0;
            m_num_total_integers = 0;
        }

        void process(uint32_t const* colors, uint64_t list_size) {
            /* encode list_size */
            bits::util::write_delta(m_bvb, list_size);
            if (list_size < m_sparse_set_threshold_size) {
                uint32_t prev_val = colors[0];
                bits::util::write_delta(m_bvb, prev_val);
                for (uint64_t i = 1; i != list_size; ++i) {
                    uint32_t val = colors[i];
                    assert(val >= prev_val + 1);
                    bits::util::write_delta(m_bvb, val - (prev_val + 1));
                    prev_val = val;
                }
            } else if (list_size < m_very_dense_set_threshold_size) {
                bits::bit_vector::builder bvb_ints;
                bvb_ints.resize(m_num_colors);
                for (uint64_t i = 0; i != list_size; ++i) bvb_ints.set(colors[i]);
                m_bvb.append(bvb_ints);
            } else {
                bool first = true;
                uint32_t val = 0;
                uint32_t prev_val = -1;
                uint32_t written = 0;
                for (uint64_t i = 0; i != list_size; ++i) {
                    uint32_t x = colors[i];
                    while (val < x) {
                        if (first) {
                            bits::util::write_delta(m_bvb, val);
                            first = false;
                            ++written;
                        } else {
                            assert(val >= prev_val + 1);
                            bits::util::write_delta(m_bvb, val - (prev_val + 1));
                            ++written;
                        }
                        prev_val = val;
                        ++val;
                    }
                    assert(val == x);
                    val++;  // skip x
                }
                while (val < m_num_colors) {
                    assert(val >= prev_val + 1);
                    bits::util::write_delta(m_bvb, val - (prev_val + 1));
                    prev_val = val;
                    ++val;
                    ++written;
                }
                assert(val == m_num_colors);
                /* complementary_list_size = m_num_colors - list_size */
                assert(m_num_colors - list_size <= m_num_colors);
                assert(written == m_num_colors - list_size);
            }
            m_offsets.push_back(m_bvb.num_bits());
            m_num_total_integers += list_size;
            m_num_lists += 1;
            if (m_num_lists % 500000 == 0) {
                std::cout << "  processed " << m_num_lists << " lists" << std::endl;
            }
        }

        void append(hybrid::builder& hb) {
            if (hb.m_num_lists == 0) return;

            m_bvb.append(hb.m_bvb);
            assert(m_offsets.size() > 0);
            uint64_t delta = m_offsets.back();

            m_offsets.reserve(m_offsets.size() + hb.m_offsets.size());
            for (uint64_t offset_id = 1; offset_id < hb.m_offsets.size(); offset_id++) {
                m_offsets.push_back(hb.m_offsets[offset_id] + delta);
            }
            m_num_lists += hb.m_num_lists;
            m_num_total_integers += hb.m_num_total_integers;
            assert(m_num_lists == m_offsets.size() - 1);
        }

        void build(hybrid& h) {
            h.m_num_colors = m_num_colors;
            h.m_sparse_set_threshold_size = m_sparse_set_threshold_size;
            h.m_very_dense_set_threshold_size = m_very_dense_set_threshold_size;

            std::cout << "processed " << m_num_lists << " color sets" << std::endl;
            std::cout << "m_num_total_integers " << m_num_total_integers << std::endl;
            assert(m_num_lists == m_offsets.size() - 1);

            h.m_offsets.encode(m_offsets.begin(), m_offsets.size(), m_offsets.back());
            m_bvb.build(h.m_color_sets);

            std::cout << "  total bits for ints = " << 8 * h.m_color_sets.num_bytes() << std::endl;
            std::cout << "  total bits per offsets = " << 8 * h.m_offsets.num_bytes() << std::endl;
            std::cout << "  total bits = "
                      << 8 * (h.m_color_sets.num_bytes() + h.m_offsets.num_bytes()) << std::endl;
            std::cout << "  offsets: " << (8.0 * h.m_offsets.num_bytes()) / m_num_total_integers
                      << " bits/int" << std::endl;
            std::cout << "  color sets: "
                      << (8.0 * h.m_color_sets.num_bytes()) / m_num_total_integers << " bits/int"
                      << std::endl;
        }

        void clear() {
            m_offsets.clear();
            m_bvb.clear();
            init(m_num_colors);
        }

    private:
        uint32_t m_num_colors;
        uint32_t m_sparse_set_threshold_size;
        uint32_t m_very_dense_set_threshold_size;
        uint64_t m_num_lists;
        uint64_t m_num_total_integers;

        bits::bit_vector::builder m_bvb;
        std::vector<uint64_t> m_offsets;
    };

    struct forward_iterator {
        forward_iterator() {}

        forward_iterator(hybrid const* ptr, uint64_t begin)
            : m_ptr(ptr)
            , m_bitmap_begin(begin)
            , m_color_sets_begin(begin)
            , m_num_colors(ptr->m_num_colors) {
            rewind();
        }

        void rewind() {
            m_pos_in_list = 0;
            m_pos_in_comp_list = 0;
            m_comp_list_size = 0;
            m_comp_val = -1;
            m_prev_val = -1;
            m_curr_val = 0;
            m_it = (m_ptr->m_color_sets).get_iterator_at(m_color_sets_begin);
            m_size = bits::util::read_delta(m_it);
            /* set m_type and read the first value */
            if (m_size < m_ptr->m_sparse_set_threshold_size) {
                m_type = list_type::delta_gaps;
                m_curr_val = bits::util::read_delta(m_it);
            } else if (m_size < m_ptr->m_very_dense_set_threshold_size) {
                m_type = list_type::bitmap;
                m_bitmap_begin = m_it.position();  // after m_size
                m_it.skip_to(m_bitmap_begin);
                uint64_t pos = m_it.next();
                assert(pos >= m_bitmap_begin);
                m_curr_val = pos - m_bitmap_begin;
            } else {
                m_type = list_type::complement_delta_gaps;
                m_comp_list_size = m_num_colors - m_size;
                if (m_comp_list_size > 0) m_comp_val = bits::util::read_delta(m_it);
                next_comp_val();
            }
        }

        /* this is needed to annul the next_comp_val() done in the constructor
           if we want to iterate through the complemented set */
        void reinit_for_complemented_set_iteration() {
            assert(m_type == list_type::complement_delta_gaps);
            m_pos_in_comp_list = 0;
            m_prev_val = -1;
            m_curr_val = 0;
            m_it = (m_ptr->m_color_sets).get_iterator_at(m_color_sets_begin);
            bits::util::read_delta(m_it); /* skip m_size */
            if (m_comp_list_size > 0) {
                m_comp_val = bits::util::read_delta(m_it);
            } else {
                m_comp_val = m_num_colors;
            }
        }

        uint64_t value() const { return m_curr_val; }
        uint64_t comp_value() const { return m_comp_val; }
        uint64_t operator*() const { return value(); }

        void next() {
            if (m_type == list_type::complement_delta_gaps) {
                ++m_curr_val;
                if (m_curr_val >= m_num_colors) {  // saturate
                    m_curr_val = m_num_colors;
                    return;
                }
                next_comp_val();
            } else if (m_type == list_type::delta_gaps) {
                m_pos_in_list += 1;
                if (m_pos_in_list >= m_size) {  // saturate
                    m_curr_val = m_num_colors;
                    return;
                }
                m_prev_val = m_curr_val;
                m_curr_val = bits::util::read_delta(m_it) + (m_prev_val + 1);
            } else {
                assert(m_type == list_type::bitmap);
                m_pos_in_list += 1;
                if (m_pos_in_list >= m_size) {  // saturate
                    m_curr_val = m_num_colors;
                    return;
                }
                uint64_t pos = m_it.next();
                assert(pos >= m_bitmap_begin);
                m_curr_val = pos - m_bitmap_begin;
            }
        }

        void next_comp() {
            ++m_pos_in_comp_list;
            if (m_pos_in_comp_list >= m_comp_list_size) {  // saturate
                m_comp_val = m_num_colors;
                return;
            }
            m_prev_val = m_comp_val;
            m_comp_val = bits::util::read_delta(m_it) + (m_prev_val + 1);
        }

        void operator++() { next(); }

        /* update the state of the iterator to the element
           which is greater-than or equal-to lower_bound */
        void next_geq(const uint64_t lower_bound) {
            assert(lower_bound <= num_colors());
            if (m_type == list_type::complement_delta_gaps) {
                if (value() > lower_bound) return;
                next_geq_comp_val(lower_bound);
                m_curr_val = lower_bound + (m_comp_val == lower_bound);
            } else {
                while (value() < lower_bound) next();
            }
            assert(value() >= lower_bound);
        }

        uint32_t size() const { return m_size; }
        uint32_t num_colors() const { return m_num_colors; }
        int type() const { return m_type; }

    private:
        hybrid const* m_ptr;
        uint64_t m_bitmap_begin;
        uint64_t m_color_sets_begin;
        uint32_t m_num_colors;
        int m_type;

        bits::bit_vector::iterator m_it;
        uint32_t m_pos_in_list;
        uint32_t m_size;

        uint32_t m_pos_in_comp_list;
        uint32_t m_comp_list_size;

        uint32_t m_comp_val;
        uint32_t m_prev_val;
        uint32_t m_curr_val;

        void next_comp_val() {
            while (m_curr_val == m_comp_val) {
                ++m_curr_val;
                ++m_pos_in_comp_list;
                if (m_pos_in_comp_list >= m_comp_list_size) break;
                m_prev_val = m_comp_val;
                m_comp_val = bits::util::read_delta(m_it) + (m_prev_val + 1);
            }
        }

        void next_geq_comp_val(const uint64_t lower_bound) {
            while (m_comp_val < lower_bound) {
                ++m_pos_in_comp_list;
                if (m_pos_in_comp_list >= m_comp_list_size) break;
                m_prev_val = m_comp_val;
                m_comp_val = bits::util::read_delta(m_it) + (m_prev_val + 1);
            }
        }
    };

    typedef forward_iterator iterator_type;

    forward_iterator color_set(uint64_t color_set_id) const {
        assert(color_set_id < num_color_sets());
        uint64_t begin = m_offsets.access(color_set_id);
        return forward_iterator(this, begin);
    }

    uint32_t num_colors() const { return m_num_colors; }
    uint64_t num_color_sets() const { return m_offsets.size() - 1; }

    uint64_t num_bits() const {
        return (sizeof(m_num_colors) + sizeof(m_sparse_set_threshold_size) +
                sizeof(m_very_dense_set_threshold_size) + m_offsets.num_bytes() +
                m_color_sets.num_bytes()) *
               8;
    }

    void print_stats() const;

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_num_colors);
        visitor.visit(t.m_sparse_set_threshold_size);
        visitor.visit(t.m_very_dense_set_threshold_size);
        visitor.visit(t.m_offsets);
        visitor.visit(t.m_color_sets);
    }

    uint32_t m_num_colors;
    uint32_t m_sparse_set_threshold_size;
    uint32_t m_very_dense_set_threshold_size;

    bits::elias_fano<false, false> m_offsets;
    bits::bit_vector m_color_sets;
};

}  // namespace fulgor

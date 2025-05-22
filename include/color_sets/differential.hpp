#pragma once

namespace fulgor {

struct differential {
    static const index_t type = index_t::DIFF;

    struct builder {
        builder() : m_prev_cluster_id(0) {
            m_list_offsets.push_back(0);
            m_representative_offsets.push_back(0);
        }

        void init_colors_builder(uint64_t num_colors) {
            m_num_colors = num_colors;
            m_num_total_integers = 0;
            m_num_lists = 0;
        }

        void encode_representative(std::vector<uint32_t> const& representative) {
            uint64_t size = representative.size();
            bits::util::write_delta(m_bvb, size);
            m_num_total_integers += size + 1;  // size plus size number
            m_num_lists += 1;
            if (size > 0) {
                uint32_t prev_val = representative[0];
                bits::util::write_delta(m_bvb, prev_val);
                for (uint64_t i = 1; i < size; ++i) {
                    uint32_t val = representative[i];
                    assert(val >= prev_val + 1);
                    bits::util::write_delta(m_bvb, val - (prev_val + 1));
                    prev_val = val;
                }
            }
            m_representative_offsets.push_back(m_bvb.num_bits());
        }

        void encode_list(uint64_t cluster_id, std::vector<uint32_t> const& representative,
                         uint64_t it_size, function<void()> next, function<uint64_t()> get)  //
        {
            std::vector<uint32_t> differential_list;
            uint64_t ref_size = representative.size();
            differential_list.reserve(ref_size + it_size);

            if (cluster_id != m_prev_cluster_id) {
                m_prev_cluster_id = cluster_id;
                m_clusters.set(m_clusters.num_bits() - 1);
            }
            m_clusters.push_back(false);

            uint64_t i = 0, j = 0;
            while (i < it_size && j < ref_size) {
                if (get() == representative[j]) {
                    i += 1;
                    j += 1;
                    next();
                } else if (get() < representative[j]) {
                    differential_list.push_back(get());
                    i += 1;
                    next();
                } else {
                    differential_list.push_back(representative[j]);
                    j += 1;
                }
            }
            while (i < it_size) {
                differential_list.push_back(get());
                next();
                i += 1;
            }
            while (j < ref_size) {
                differential_list.push_back(representative[j]);
                j += 1;
            }

            uint64_t size = differential_list.size();
            bits::util::write_delta(m_bvb, size);
            bits::util::write_delta(m_bvb, it_size);

            // size plus differential_list size plus original list size
            m_num_total_integers += size + 2;

            m_num_lists += 1;

            if (size > 0) {
                uint32_t prev_val = differential_list[0];
                bits::util::write_delta(m_bvb, prev_val);
                for (uint64_t pos = 1; pos < size; ++pos) {
                    uint32_t val = differential_list[pos];
                    assert(val >= prev_val + 1);
                    bits::util::write_delta(m_bvb, val - (prev_val + 1));
                    prev_val = val;
                }
            }

            uint64_t last_offset = m_representative_offsets[m_representative_offsets.size() - 1];
            m_list_offsets.push_back(m_bvb.num_bits() - last_offset);
        }

        void build(differential& d) {
            d.m_num_colors = m_num_colors;
            m_bvb.build(d.m_color_sets);
            m_clusters.build(d.m_clusters);
            d.m_clusters_rank1_index.build(d.m_clusters);
            d.m_representative_offsets.encode(m_representative_offsets.begin(),  //
                                              m_representative_offsets.size(),   //
                                              m_representative_offsets.back());
            d.m_list_offsets.encode(m_list_offsets.begin(),  //
                                    m_list_offsets.size(),   //
                                    m_list_offsets.back());

            std::cout << "processed " << m_num_lists << " color sets\n";
            std::cout << "m_num_total_integers " << m_num_total_integers << '\n';

            std::cout << "  total bits for ints = " << 8 * d.m_color_sets.num_bytes() << '\n';
            std::cout << "  total bits per offset = "
                      << 8 * (d.m_list_offsets.num_bytes() + d.m_representative_offsets.num_bytes())
                      << " (differences: " << 8 * d.m_list_offsets.num_bytes()
                      << ", representatives: " << 8 * d.m_representative_offsets.num_bytes()
                      << ")\n";
            std::cout << "  offsets: "
                      << 8.0 *
                             (d.m_list_offsets.num_bytes() +
                              d.m_representative_offsets.num_bytes()) /
                             m_num_total_integers
                      << " bits/int\n";
            std::cout << "  color sets: "
                      << (8.0 * d.m_color_sets.num_bytes()) / m_num_total_integers << " bits/int\n";
        }

    private:
        bits::bit_vector::builder m_bvb, m_clusters;
        uint64_t m_num_total_integers, m_num_lists;

        uint64_t m_num_colors;
        uint64_t m_prev_cluster_id;
        std::vector<uint64_t> m_representative_offsets, m_list_offsets;
    };

    struct forward_iterator {
        forward_iterator() {}

        forward_iterator(differential const* ptr, uint64_t list_begin,
                         uint64_t representative_begin)
            : m_ptr(ptr)
            , m_differential_list_begin(list_begin)
            , m_representative_begin(representative_begin) {
            rewind();
        }

        void rewind() {
            init();
            update_curr_val();
        }

        void full_rewind() { init(); }

        uint32_t size() const { return m_size; }

        uint64_t value() const { return m_curr_val; }
        uint64_t operator*() const { return value(); }

        void next() {
            if (m_pos_in_representative >= m_representative_size &&
                m_pos_in_differential_list >= m_differential_list_size) {
                m_curr_val = num_colors();
                return;
            }
            if (m_pos_in_representative >= m_representative_size ||
                m_curr_differential_val < m_curr_representative_val) {
                next_differential_val();
            } else if (m_pos_in_differential_list >= m_differential_list_size ||
                       m_curr_representative_val < m_curr_differential_val) {
                next_representative_val();
            }
            update_curr_val();
        }
        void operator++() { next(); }

        void next_geq(const uint64_t lower_bound) {
            assert(lower_bound <= num_colors());
            while (value() < lower_bound) next();
            assert(value() >= lower_bound);
        }

        uint32_t num_colors() const { return m_ptr->m_num_colors; }
        uint64_t differential_list_size() const { return m_differential_list_size; }
        int encoding_type() const { return encoding_t::symmetric_difference; }

        uint64_t representative_begin() const { return m_representative_begin; }

        void next_representative_val() {
            m_pos_in_representative += 1;
            m_prev_representative_val = m_curr_representative_val;
            if (m_pos_in_representative < m_representative_size) {
                m_curr_representative_val =
                    m_prev_representative_val + bits::util::read_delta(m_representative_it) + 1;
            } else {
                m_curr_representative_val = num_colors();
            }
        }

        uint32_t representative_val() const { return m_curr_representative_val; }

        void next_differential_val() {
            m_pos_in_differential_list += 1;
            m_prev_differential_val = m_curr_differential_val;
            if (m_pos_in_differential_list < m_differential_list_size) {
                m_curr_differential_val =
                    m_prev_differential_val + bits::util::read_delta(m_differential_list_it) + 1;
            } else {
                m_curr_differential_val = num_colors();
            }
        }

        uint32_t differential_val() const { return m_curr_differential_val; }

    private:
        differential const* m_ptr;
        uint64_t m_differential_list_begin, m_representative_begin;
        uint64_t m_representative_size, m_differential_list_size;
        uint64_t m_pos_in_differential_list, m_pos_in_representative;
        uint32_t m_curr_representative_val, m_curr_differential_val;
        uint32_t m_prev_representative_val, m_prev_differential_val;
        uint32_t m_curr_val;
        uint32_t m_size;
        bits::bit_vector::iterator m_representative_it, m_differential_list_it;

        void init() {
            m_differential_list_it =  //
                (m_ptr->m_color_sets).get_iterator_at(m_differential_list_begin);
            m_representative_it =  //
                (m_ptr->m_color_sets).get_iterator_at(m_representative_begin);

            m_differential_list_size = bits::util::read_delta(m_differential_list_it);
            m_representative_size = bits::util::read_delta(m_representative_it);
            m_size = bits::util::read_delta(m_differential_list_it);

            m_curr_differential_val = m_differential_list_size == 0
                                          ? num_colors()
                                          : bits::util::read_delta(m_differential_list_it);
            m_prev_differential_val = 0;
            m_curr_representative_val = m_representative_size == 0
                                            ? num_colors()
                                            : bits::util::read_delta(m_representative_it);
            m_prev_representative_val = 0;

            m_pos_in_differential_list = 0;
            m_pos_in_representative = 0;
        }

        void update_curr_val() {
            while (m_curr_representative_val == m_curr_differential_val &&
                   m_pos_in_representative <= m_representative_size &&
                   m_pos_in_differential_list <= m_differential_list_size) {
                next_differential_val();
                next_representative_val();
            }
            m_curr_val = min(m_curr_differential_val, m_curr_representative_val);
        }
    };

    typedef forward_iterator iterator_type;

    forward_iterator color_set(uint64_t color_id) const {
        assert(color_id < num_color_sets());
        uint64_t last_representative = m_representative_offsets.access(num_partitions());
        uint64_t list_begin = m_list_offsets.access(color_id) + last_representative;
        uint64_t representative_begin =
            m_representative_offsets.access(m_clusters_rank1_index.rank1(m_clusters, color_id));
        return forward_iterator(this, list_begin, representative_begin);
    }

    uint64_t num_color_sets() const { return m_list_offsets.size() - 1; }
    uint64_t num_partitions() const { return m_clusters_rank1_index.num_ones() + 1; }
    uint64_t num_colors() const { return m_num_colors; }

    uint64_t num_bits() const {
        return (sizeof(m_num_colors) + m_representative_offsets.num_bytes() +
                m_list_offsets.num_bytes() + m_color_sets.num_bytes() + m_clusters.num_bytes() +
                m_clusters_rank1_index.num_bytes()) *
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
        visitor.visit(t.m_representative_offsets);
        visitor.visit(t.m_list_offsets);
        visitor.visit(t.m_color_sets);
        visitor.visit(t.m_clusters);
        visitor.visit(t.m_clusters_rank1_index);
    }

    uint32_t m_num_colors;
    bits::elias_fano<false, false> m_representative_offsets, m_list_offsets;
    bits::bit_vector m_color_sets;
    bits::bit_vector m_clusters;
    bits::rank9 m_clusters_rank1_index;
};

}  // namespace fulgor

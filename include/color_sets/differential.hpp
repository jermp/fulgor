#pragma once

namespace fulgor {

struct differential {
    static const index_t type = index_t::DIFF;

    struct builder {
        builder()
            : m_num_total_integers(0)
            , m_num_sets(0) { }
        builder(uint32_t num_colors)
            : m_num_total_integers(0)
            , m_num_sets(0)
            , m_num_colors(num_colors) { }

        void init_color_sets_builder(uint64_t num_colors) {
            m_num_colors = num_colors;
            m_num_total_integers = 0;
            m_num_sets = 0;
        }

        void reserve_num_bits(uint64_t num_bits) { m_bvb.reserve(num_bits); }

        void process_partition(std::vector<uint32_t> const& representative){
            m_representative_offsets.push_back(m_bvb.num_bits());
            m_curr_representative = representative;
            if (m_clusters.num_bits() > 0) m_clusters.set(m_clusters.num_bits() - 1);

            uint64_t size = representative.size();

            bits::util::write_delta(m_bvb, size);
            m_num_total_integers += size + 1;  // size plus size number
            m_num_sets += 1;

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
        }

        template<typename Iterator>
        void process_color_set(Iterator& it){
            m_color_set_offsets.push_back(m_bvb.num_bits());
            uint64_t it_size = it.size();
            uint64_t rp_size = m_curr_representative.size();
            std::vector<uint32_t> diff_set;
            diff_set.reserve(m_num_colors);

            m_clusters.push_back(false);

            uint64_t i = 0, j = 0;
            while (i < it_size && j < rp_size) {
                if (*it == m_curr_representative[j]) {
                    i += 1;
                    j += 1;
                    ++it;
                } else if (*it < m_curr_representative[j]) {
                    diff_set.push_back(*it);
                    i += 1;
                    ++it;
                } else {
                    diff_set.push_back(m_curr_representative[j]);
                    j += 1;
                }
            }
            while (i < it_size) {
                diff_set.push_back(*it);
                ++it;
                i += 1;
            }
            while (j < rp_size) {
                diff_set.push_back(m_curr_representative[j]);
                j += 1;
            }

            uint64_t size = diff_set.size();
            bits::util::write_delta(m_bvb, size);
            bits::util::write_delta(m_bvb, it_size);

            // size plus diff_set size plus original set size
            m_num_total_integers += size + 2;

            m_num_sets += 1;

            if (size > 0) {
                uint32_t prev_val = diff_set[0];
                bits::util::write_delta(m_bvb, prev_val);
                for (uint64_t pos = 1; pos < size; ++pos) {
                    uint32_t val = diff_set[pos];
                    assert(val >= prev_val + 1);
                    bits::util::write_delta(m_bvb, val - (prev_val + 1));
                    prev_val = val;
                }
            }

        }

        void append(differential::builder const& db) {
            if (db.m_color_set_offsets.size() - 1 == 0) return;
            uint64_t delta = m_bvb.num_bits();
            m_bvb.append(db.m_bvb);
            m_num_total_integers += db.m_num_total_integers;
            m_num_sets += db.m_num_sets;
            m_clusters.set(m_clusters.num_bits() - 1);
            m_clusters.append(db.m_clusters);

            assert(m_representative_offsets.size() > 0);
            m_representative_offsets.reserve(m_representative_offsets.size() + db.m_representative_offsets.size());
            for (uint64_t i = 0; i != db.m_representative_offsets.size(); ++i){
                m_representative_offsets.push_back(db.m_representative_offsets[i] + delta);
            }

            assert(m_color_set_offsets.size() > 0);
            m_color_set_offsets.reserve(m_color_set_offsets.size() + db.m_color_set_offsets.size());
            for (uint64_t i = 0; i != db.m_color_set_offsets.size(); ++i){
                m_color_set_offsets.push_back(db.m_color_set_offsets[i] + delta);
            }

        }

        void build(differential& d) {
            d.m_num_colors = m_num_colors;
            m_bvb.build(d.m_color_sets);

            m_clusters.set(m_clusters.num_bits() - 1);
            m_clusters.build(d.m_clusters);
            d.m_clusters_rank1_index.build(d.m_clusters);
            std::cout << "Processing representatives" << std::endl;
            d.m_representative_offsets.encode(m_representative_offsets.begin(),  //
                                              m_representative_offsets.size(),   //
                                              m_representative_offsets.back());
            std::cout << "Processing differential sets" << std::endl;
            d.m_color_set_offsets.encode(m_color_set_offsets.begin(),  //
                                         m_color_set_offsets.size(),   //
                                         m_color_set_offsets.back());

            std::cout << "processed " << m_num_sets << " color sets\n";
            std::cout << "m_num_total_integers " << m_num_total_integers << '\n';

            std::cout << "  total bits for ints = " << 8 * d.m_color_sets.num_bytes() << '\n';
            std::cout << "  total bits per offset = "
                      << 8 * (d.m_color_set_offsets.num_bytes() +
                              d.m_representative_offsets.num_bytes())
                      << " (differences: " << 8 * d.m_color_set_offsets.num_bytes()
                      << ", representatives: " << 8 * d.m_representative_offsets.num_bytes()
                      << ")\n";
            std::cout << "  offsets: "
                      << 8.0 *
                             (d.m_color_set_offsets.num_bytes() +
                              d.m_representative_offsets.num_bytes()) /
                             m_num_total_integers
                      << " bits/int\n";
            std::cout << "  color sets: "
                      << (8.0 * d.m_color_sets.num_bytes()) / m_num_total_integers << " bits/int\n";
        }

    private:
        bits::bit_vector::builder m_bvb, m_clusters;
        uint64_t m_num_total_integers, m_num_sets;

        uint64_t m_num_colors;
        std::vector<uint64_t> m_representative_offsets, m_color_set_offsets;

        std::vector<uint32_t> m_curr_representative;
    };

    struct forward_iterator {
        forward_iterator() {}

        forward_iterator(differential const* ptr, uint64_t set_begin, uint64_t representative_begin)
            : m_ptr(ptr)
            , m_differential_set_begin(set_begin)
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
                m_pos_in_differential_set >= m_differential_set_size) {
                m_curr_val = num_colors();
                return;
            }
            if (m_pos_in_representative >= m_representative_size ||
                m_curr_differential_val < m_curr_representative_val) {
                next_differential_val();
            } else if (m_pos_in_differential_set >= m_differential_set_size ||
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
        uint64_t differential_set_size() const { return m_differential_set_size; }
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
            m_pos_in_differential_set += 1;
            m_prev_differential_val = m_curr_differential_val;
            if (m_pos_in_differential_set < m_differential_set_size) {
                m_curr_differential_val =
                    m_prev_differential_val + bits::util::read_delta(m_differential_set_it) + 1;
            } else {
                m_curr_differential_val = num_colors();
            }
        }

        uint32_t differential_val() const { return m_curr_differential_val; }

    private:
        differential const* m_ptr;
        uint64_t m_differential_set_begin, m_representative_begin;
        uint64_t m_representative_size, m_differential_set_size;
        uint64_t m_pos_in_differential_set, m_pos_in_representative;
        uint32_t m_curr_representative_val, m_curr_differential_val;
        uint32_t m_prev_representative_val, m_prev_differential_val;
        uint32_t m_curr_val;
        uint32_t m_size;
        bits::bit_vector::iterator m_representative_it, m_differential_set_it;

        void init() {
            m_differential_set_it =  //
                (m_ptr->m_color_sets).get_iterator_at(m_differential_set_begin);
            m_representative_it =  //
                (m_ptr->m_color_sets).get_iterator_at(m_representative_begin);

            m_differential_set_size = bits::util::read_delta(m_differential_set_it);
            m_representative_size = bits::util::read_delta(m_representative_it);
            m_size = bits::util::read_delta(m_differential_set_it);

            m_curr_differential_val = m_differential_set_size == 0
                                          ? num_colors()
                                          : bits::util::read_delta(m_differential_set_it);
            m_prev_differential_val = 0;
            m_curr_representative_val = m_representative_size == 0
                                            ? num_colors()
                                            : bits::util::read_delta(m_representative_it);
            m_prev_representative_val = 0;

            m_pos_in_differential_set = 0;
            m_pos_in_representative = 0;
        }

        void update_curr_val() {
            while (m_curr_representative_val == m_curr_differential_val &&
                   m_pos_in_representative <= m_representative_size &&
                   m_pos_in_differential_set <= m_differential_set_size) {
                next_differential_val();
                next_representative_val();
            }
            m_curr_val = min(m_curr_differential_val, m_curr_representative_val);
        }
    };

    typedef forward_iterator iterator_type;

    forward_iterator color_set(uint64_t color_id) const {
        assert(color_id < num_color_sets());
        uint64_t set_begin = m_color_set_offsets.access(color_id);
        uint64_t representative_begin =
            m_representative_offsets.access(m_clusters_rank1_index.rank1(m_clusters, color_id));
        return forward_iterator(this, set_begin, representative_begin);
    }

    uint64_t num_color_sets() const { return m_color_set_offsets.size(); }
    uint64_t num_partitions() const { return m_representative_offsets.size(); }
    uint64_t num_colors() const { return m_num_colors; }

    uint64_t num_bits() const {
        return (sizeof(m_num_colors) + m_representative_offsets.num_bytes() +
                m_color_set_offsets.num_bytes() + m_color_sets.num_bytes() +
                m_clusters.num_bytes() + m_clusters_rank1_index.num_bytes()) *
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
        visitor.visit(t.m_color_set_offsets);
        visitor.visit(t.m_color_sets);
        visitor.visit(t.m_clusters);
        visitor.visit(t.m_clusters_rank1_index);
    }

    uint32_t m_num_colors;
    bits::elias_fano<false, false> m_representative_offsets, m_color_set_offsets;
    bits::bit_vector m_color_sets;
    bits::bit_vector m_clusters;
    bits::rank9 m_clusters_rank1_index;
};

}  // namespace fulgor

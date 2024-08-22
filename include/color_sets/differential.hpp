#pragma once

namespace fulgor {

struct differential {
    static const bool meta_colored = false;
    static const bool differential_colored = true;

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
            util::write_delta(m_bvb, size);
            m_num_total_integers += size + 1;  // size plus size number
            m_num_lists += 1;

            if (size > 0) {
                uint32_t prev_val = representative[0];
                util::write_delta(m_bvb, prev_val);

                for (uint64_t i = 1; i < size; ++i) {
                    uint32_t val = representative[i];
                    assert(val >= prev_val + 1);
                    util::write_delta(m_bvb, val - (prev_val + 1));
                    prev_val = val;
                }
            }
            m_representative_offsets.push_back(m_bvb.num_bits());
        }

        void encode_list(uint64_t cluster_id, std::vector<uint32_t> const& representative,
                         uint64_t it_size, function<void()> next, function<uint64_t()> get) {
            std::vector<uint32_t> differential_list;
            uint64_t ref_size = representative.size();
            differential_list.reserve(ref_size + it_size);

            if (cluster_id != m_prev_cluster_id) {
                m_prev_cluster_id = cluster_id;
                m_clusters.set(m_clusters.size() - 1);
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
            util::write_delta(m_bvb, size);
            util::write_delta(m_bvb, it_size);
            m_num_total_integers +=
                size + 2;  // size plus differential_list size plus original list size
            m_num_lists += 1;

            if (size > 0) {
                uint32_t prev_val = differential_list[0];
                util::write_delta(m_bvb, prev_val);

                for (uint64_t pos = 1; pos < size; ++pos) {
                    uint32_t val = differential_list[pos];
                    assert(val >= prev_val + 1);
                    util::write_delta(m_bvb, val - (prev_val + 1));
                    prev_val = val;
                }
            }

            uint64_t last_offset = m_representative_offsets[m_representative_offsets.size() - 1];
            m_list_offsets.push_back(m_bvb.num_bits() - last_offset);
        }

        void build(differential& d) {
            d.m_num_colors = m_num_colors;
            d.m_colors.swap(m_bvb.bits());
            d.m_clusters.build(&m_clusters);

            d.m_representative_offsets.encode(m_representative_offsets.begin(),
                                              m_representative_offsets.size(),
                                              m_representative_offsets.back());
            d.m_list_offsets.encode(m_list_offsets.begin(), m_list_offsets.size(),
                                    m_list_offsets.back());

            std::cout << "processed " << m_num_lists << " lists\n";
            std::cout << "m_num_total_integers " << m_num_total_integers << '\n';

            std::cout << "  total bits for ints = " << d.m_colors.size() * 64 << '\n';
            std::cout << "  total bits per offset = "
                      << d.m_list_offsets.num_bits() + d.m_representative_offsets.num_bits()
                      << " (lists: " << d.m_list_offsets.num_bits()
                      << ", representatives: " << d.m_representative_offsets.num_bits() << ")\n";
            std::cout << "  offsets: "
                      << static_cast<double>(d.m_list_offsets.num_bits() +
                                             d.m_representative_offsets.num_bits()) /
                             m_num_total_integers
                      << " bits/int\n";
            std::cout << "  lists: "
                      << static_cast<double>(d.m_colors.size() * 64) / m_num_total_integers
                      << " bits/int\n";
        }

    private:
        bit_vector_builder m_bvb;
        pthash::bit_vector_builder m_clusters;
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
            m_differential_list_it = bit_vector_iterator(
                (m_ptr->m_colors).data(), (m_ptr->m_colors).size(), m_differential_list_begin);
            m_representative_it = bit_vector_iterator(
                (m_ptr->m_colors).data(), (m_ptr->m_colors).size(), m_representative_begin);
            m_differential_list_size = util::read_delta(m_differential_list_it);
            m_representative_size = util::read_delta(m_representative_it);
            m_size = util::read_delta(m_differential_list_it);

            m_curr_differential_val = m_differential_list_size == 0
                                          ? num_colors()
                                          : util::read_delta(m_differential_list_it);
            m_prev_differential_val = 0;
            m_curr_representative_val =
                m_representative_size == 0 ? num_colors() : util::read_delta(m_representative_it);
            m_prev_representative_val = 0;

            m_pos_in_differential_list = 0;
            m_pos_in_representative = 0;
            update_curr_val();
        }

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

        int type() const { return list_type::differential_list; }

    private:
        differential const* m_ptr;
        uint64_t m_differential_list_begin, m_representative_begin;
        uint64_t m_representative_size, m_differential_list_size;
        uint64_t m_pos_in_differential_list, m_pos_in_representative;
        uint32_t m_curr_representative_val, m_curr_differential_val;
        uint32_t m_prev_representative_val, m_prev_differential_val;
        uint32_t m_curr_val;
        uint32_t m_size;
        bit_vector_iterator m_representative_it, m_differential_list_it;

        void next_representative_val() {
            m_pos_in_representative += 1;
            m_prev_representative_val = m_curr_representative_val;
            if (m_pos_in_representative < m_representative_size) {
                m_curr_representative_val =
                    m_prev_representative_val + util::read_delta(m_representative_it) + 1;
            } else {
                m_curr_representative_val = num_colors();
            }
        }

        void next_differential_val() {
            m_pos_in_differential_list += 1;
            m_prev_differential_val = m_curr_differential_val;
            if (m_pos_in_differential_list < m_differential_list_size) {
                m_curr_differential_val =
                    m_prev_differential_val + util::read_delta(m_differential_list_it) + 1;
            } else {
                m_curr_differential_val = num_colors();
            }
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
        uint64_t representative_begin = m_representative_offsets.access(m_clusters.rank(color_id));
        return forward_iterator(this, list_begin, representative_begin);
    }

    uint64_t num_color_sets() const { return m_list_offsets.size() - 1; }
    uint64_t num_partitions() const { return m_clusters.num_ones() + 1; }
    uint64_t num_colors() const { return m_num_colors; }

    uint64_t num_bits() const {
        return sizeof(m_num_colors) * 8 + m_representative_offsets.num_bits() +
               m_list_offsets.num_bits() + essentials::vec_bytes(m_colors) * 8 +
               m_clusters.bytes() * 8;
    }

    void print_stats() const {
        std::cout << "Color statistics:\n";
        std::cout << "  Number of partitions: " << num_partitions() << std::endl;

        uint64_t num_bits_representative_offsets = m_representative_offsets.num_bits();
        uint64_t num_bits_list_offsets = m_list_offsets.num_bits();
        uint64_t num_bits_colors = essentials::vec_bytes(m_colors) * 8;

        uint64_t num_clusters = m_clusters.size();
        uint64_t num_representatives = 0;
        uint64_t num_differential_lists = 0;
        uint64_t num_metadata = 0;

        uint64_t num_colors_tenth = num_colors() / 10;

        std::vector<uint64_t> distribution(11, 0);

        for (uint64_t representative_id = 0; representative_id < num_partitions();
             representative_id++) {
            uint64_t representative_begin = m_representative_offsets.access(representative_id);
            auto it = bit_vector_iterator(m_colors.data(), m_colors.size(), representative_begin);
            uint64_t prev_position = it.position();

            uint64_t size = util::read_delta(it);
            num_metadata += it.position() - prev_position;
            prev_position = it.position();

            for (uint64_t i = 0; i < size; i++) {
                util::read_delta(it);
                num_representatives += it.position() - prev_position;
                prev_position = it.position();
            }
        }
        uint64_t last_representative = m_representative_offsets.access(num_partitions());
        for (uint64_t color_id = 0; color_id < num_color_sets(); color_id++) {
            uint64_t list_begin = m_list_offsets.access(color_id) + last_representative;
            auto it = bit_vector_iterator(m_colors.data(), m_colors.size(), list_begin);
            uint64_t prev_position = it.position();

            uint64_t size = util::read_delta(it);
            num_metadata += it.position() - prev_position;
            prev_position = it.position();

            util::read_delta(it);  // original list size
            num_metadata += it.position() - prev_position;
            prev_position = it.position();

            for (uint64_t i = 0; i < size; i++) {
                util::read_delta(it);
                uint64_t delta_size = it.position() - prev_position;
                num_differential_lists += delta_size;

                prev_position = it.position();
            }
            uint64_t q = 0;
            if (num_colors_tenth != 0) {
                q = size / (num_colors_tenth) > 10 ? 10 : size / (num_colors_tenth);
            }

            distribution[q]++;
        }

        assert(num_bits() > 0);
        assert(num_bits_colors > 0);

        std::cout << "  representative offsets: " << num_bits_representative_offsets / 8
                  << " bytes (" << (num_bits_representative_offsets * 100.0) / num_bits() << "%)"
                  << std::endl;
        std::cout << "  differential list offsets: " << num_bits_list_offsets / 8 << " bytes ("
                  << (num_bits_list_offsets * 100.0) / num_bits() << "%)" << std::endl;
        std::cout << "  clusters: " << num_clusters / 8 << " bytes ("
                  << (num_clusters * 100.0) / num_bits() << "%)" << std::endl;
        std::cout << "  differential colors: " << num_bits_colors / 8 << " bytes ("
                  << (num_bits_colors * 100.0) / num_bits() << "%)" << std::endl;
        std::cout << "    representatives: " << num_representatives / 8 << " bytes ("
                  << (num_representatives * 100.0) / num_bits_colors << "%)" << std::endl;
        std::cout << "    differential lists: " << num_differential_lists / 8 << " bytes ("
                  << (num_differential_lists * 100.0) / num_bits_colors << "%)" << std::endl;
        std::cout << "    metadata: " << num_metadata / 8 << " bytes ("
                  << (num_metadata * 100.0) / num_bits_colors << "%)" << std::endl;
        std::cout << "  differential lists size distribution:" << std::endl;
        for (uint64_t partition = 0; partition < 11; partition++) {
            std::cout << distribution[partition] << " ";
        }
        std::cout << std::endl;
    }

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
        visitor.visit(t.m_colors);
        visitor.visit(t.m_clusters);
    }

    uint32_t m_num_colors;

    sshash::ef_sequence<false> m_representative_offsets, m_list_offsets;

    std::vector<uint64_t> m_colors;
    ranked_bit_vector m_clusters;

    std::vector<uint64_t> read_representative_set(uint64_t begin) const {
        auto it = bit_vector_iterator(m_colors.data(), m_colors.size(), begin);
        uint64_t size = util::read_delta(it);

        if (size == 0) return {};

        std::vector<uint64_t> set(size);
        set[0] = util::read_delta(it);
        for (uint64_t i = 1; i != size; ++i) { set[i] = set[i - 1] + util::read_delta(it) + 1; }
        return set;
    }

    std::vector<uint64_t> read_differential_set(uint64_t begin) const {
        auto it = bit_vector_iterator(m_colors.data(), m_colors.size(), begin);
        uint64_t size = util::read_delta(it);
        util::read_delta(it);  // skip color_set size

        if (size == 0) return {};

        std::vector<uint64_t> set(size);
        set[0] = util::read_delta(it);
        for (uint64_t i = 1; i != size; ++i) { set[i] = set[i - 1] + util::read_delta(it) + 1; }
        return set;
    }
};

}  // namespace fulgor

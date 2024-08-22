#pragma once

namespace fulgor {

struct meta_differential {
    static const bool meta_colored = true;
    static const bool differential_colored = true;

    struct partition_endpoint {
        template <typename Visitor>
        void visit(Visitor& visitor) {
            visitor.visit(docid_lower_bound);
            visitor.visit(num_lists);
        }
        uint64_t docid_lower_bound;
        uint64_t num_lists;
    };

    struct builder {
        builder() : m_prev_docs(0), m_prev_partition_set_id(0) {
            m_partition_sets_offsets.push_back(0);
            m_relative_colors_offsets.push_back(0);
        }

        void init(uint64_t num_colors, uint64_t num_partitions) {
            m_num_colors = num_colors;
            m_partition_endpoints.reserve(num_partitions);
            m_partial_colors.reserve(num_partitions);
        }

        void init_meta_color_partition_sets(uint64_t num_sets) {
            m_num_partition_sets = num_sets;
            m_partition_sets_offsets.reserve(num_sets);
        }

        void process_meta_color_partition_set(vector<uint64_t>& partition_set) {
            uint64_t size = partition_set.size();
            uint64_t prev_val = partition_set[0];

            util::write_delta(m_partition_sets, size);
            util::write_delta(m_partition_sets, prev_val);
            for (uint64_t i = 1; i < size; i++) {
                assert(prev_val < partition_set[i]);
                util::write_delta(m_partition_sets, partition_set[i] - prev_val);
                prev_val = partition_set[i];
            }

            m_partition_sets_offsets.push_back(m_partition_sets.num_bits());
        }

        void process_partition(differential& d) {
            m_partial_colors.push_back(d);
            m_partition_endpoints.push_back({m_prev_docs, d.num_color_sets()});
            m_prev_docs += d.num_colors();
        }

        void process_metacolors(uint64_t partition_set_id, vector<uint64_t>& partition_set,
                                vector<uint64_t>& relative_colors) {
            assert(partition_set.size() == relative_colors.size());
            if (partition_set_id != m_prev_partition_set_id) {
                m_prev_partition_set_id = partition_set_id;
                m_partition_sets_partitions.set(m_partition_sets_partitions.size() - 1);
            }
            m_partition_sets_partitions.push_back(false);

            uint64_t size = partition_set.size();
            for (uint64_t i = 0; i < size; i++) {
                uint64_t partition_id = partition_set[i];
                uint64_t relative_id = relative_colors[i];

                uint64_t partition_size = m_partition_endpoints[partition_id].num_lists;
                m_relative_colors.append_bits(relative_id, msb(partition_size));
            }
            m_relative_colors_offsets.push_back(m_relative_colors.num_bits());
        }

        void build(meta_differential& m) {
            m.m_num_colors = m_num_colors;
            m.m_num_partition_sets = m_num_partition_sets;
            m.m_partial_colors.swap(m_partial_colors);
            m.m_relative_colors.swap(m_relative_colors.bits());
            m.m_partition_sets.swap(m_partition_sets.bits());

            m.m_partition_sets_partitions.build(&m_partition_sets_partitions);
            m.m_partition_sets_offsets.encode(m_partition_sets_offsets.begin(),
                                              m_partition_sets_offsets.size(),
                                              m_partition_sets_offsets.back());
            m.m_relative_colors_offsets.encode(m_relative_colors_offsets.begin(),
                                               m_relative_colors_offsets.size(),
                                               m_relative_colors_offsets.back());
            m.m_partition_endpoints.swap(m_partition_endpoints);
        }

    private:
        uint8_t msb(uint64_t n) { return 64 - __builtin_clzll(n); }

        std::vector<differential> m_partial_colors;
        bit_vector_builder m_relative_colors;
        bit_vector_builder m_partition_sets;

        pthash::bit_vector_builder m_partition_sets_partitions;

        std::vector<uint64_t> m_partition_sets_offsets;
        std::vector<uint64_t> m_relative_colors_offsets;

        uint64_t m_num_colors, m_prev_docs;
        uint64_t m_num_partition_sets;

        uint64_t m_prev_partition_set_id;

        std::vector<partition_endpoint> m_partition_endpoints;
    };

    struct forward_iterator {
        forward_iterator(meta_differential const* ptr, uint64_t begin_partition_set,
                         uint64_t begin_rel)
            : m_ptr(ptr), m_begin_partition_set(begin_partition_set), m_begin_rel(begin_rel) {
            rewind();
            assert(m_meta_color_list_size > 0);
        }

        void rewind() {
            init();
            change_partition();
        }

        void init() {
            m_num_lists_before = 0;
            m_pos_in_meta_color = 0;
            m_pos_in_partial_color = 0;
            m_curr_partition_id = 0;
            m_partition_set_id =
                bit_vector_iterator((m_ptr->m_partition_sets).data(),
                                    (m_ptr->m_partition_sets).size(), m_begin_partition_set);
            m_relative_colors_it = bit_vector_iterator(
                (m_ptr->m_relative_colors).data(), (m_ptr->m_relative_colors).size(), m_begin_rel);
            m_meta_color_list_size = util::read_delta(m_partition_set_id);
        }

        uint64_t value() const { return m_curr_val; }
        uint64_t operator*() const { return value(); }

        bool has_next() const { return m_pos_in_partial_color != m_curr_partition_size; }

        void next() {
            if (m_pos_in_partial_color == m_curr_partition_size - 1) {
                if (m_pos_in_meta_color == meta_color_list_size() - 1) {  // saturate
                    m_curr_val = num_colors();
                    return;
                }
                m_pos_in_meta_color += 1;
                change_partition();
            } else {
                next_in_partition();
            }
        }
        void operator++() { next(); }

        /* update the state of the iterator to the element
           which is greater-than or equal-to lower_bound */
        void next_geq(const uint64_t lower_bound) {
            assert(lower_bound <= num_colors());
            while (value() < lower_bound) next();
            assert(value() >= lower_bound);
        }

        void next_in_partition() {
            m_pos_in_partial_color += 1;
            ++m_curr_partition_it;
            update_curr_val();
        }

        void change_partition() {
            read_partition_id();
            update_partition();
        }

        void next_partition_id() {
            m_pos_in_meta_color += 1;
            if (m_pos_in_meta_color == meta_color_list_size()) {
                m_curr_partition_id = num_partitions();
                return;
            }
            read_partition_id();
        }

        void read_partition_id() {
            uint64_t delta = util::read_delta(m_partition_set_id);
            for (uint64_t i = 0; i < delta; i++) {
                m_num_lists_before +=
                    m_ptr->m_partition_endpoints[m_curr_partition_id + i].num_lists;
            }
            m_curr_partition_id += delta;
            uint8_t relative_color_size =
                msb(m_ptr->m_partition_endpoints[m_curr_partition_id].num_lists);
            m_curr_relative_color = m_relative_colors_it.take(relative_color_size);
        }

        void next_geq_partition_id(const uint32_t lower_bound) {
            assert(lower_bound <= num_partitions());
            while (partition_id() < lower_bound) next_partition_id();
            assert(partition_id() >= lower_bound);
        }

        void update_partition() {
            m_docid_lower_bound =
                m_ptr->m_partition_endpoints[m_curr_partition_id].docid_lower_bound;

            m_pos_in_partial_color = 0;
            m_curr_partition_it =
                m_ptr->m_partial_colors[m_curr_partition_id].color_set(m_curr_relative_color);
            m_curr_partition_size = m_curr_partition_it.size();

            update_curr_val();
        }

        uint64_t size() const {
            uint64_t size = 0;
            auto partition_set_it =
                bit_vector_iterator((m_ptr->m_partition_sets).data(),
                                    (m_ptr->m_partition_sets).size(), m_begin_partition_set);
            auto rel_it = bit_vector_iterator((m_ptr->m_relative_colors).data(),
                                              (m_ptr->m_relative_colors).size(), m_begin_rel);
            uint64_t partition_id = 0;
            util::read_delta(partition_set_it);  // remove size
            for (uint64_t partial_color_id = 0; partial_color_id < m_meta_color_list_size;
                 partial_color_id++) {
                partition_id += util::read_delta(partition_set_it);
                uint8_t relative_color_size =
                    msb(m_ptr->m_partition_endpoints[partition_id].num_lists);
                uint64_t relative_color = rel_it.take(relative_color_size);
                size += m_ptr->m_partial_colors[partition_id].color_set(relative_color).size();
            }
            return size;
        }

        uint32_t partition_id() const { return m_curr_partition_id; }
        uint32_t partition_upper_bound() const {
            return m_docid_lower_bound + m_curr_partition_it.num_colors();
        }
        uint32_t meta_color() const { return m_num_lists_before + m_curr_relative_color; }
        uint32_t num_colors() const { return m_ptr->num_colors(); }
        uint32_t num_partitions() const { return m_ptr->num_partitions(); }
        uint64_t meta_color_list_size() const { return m_meta_color_list_size; }

    private:
        meta_differential const* m_ptr;
        differential::iterator_type m_curr_partition_it;
        bit_vector_iterator m_partition_set_id, m_relative_colors_it;
        uint64_t m_meta_color_list_size;
        uint64_t m_begin_partition_set, m_begin_rel;
        uint64_t m_pos_in_meta_color, m_pos_in_partial_color;
        uint64_t m_curr_relative_color;
        uint64_t m_curr_partition_id, m_curr_partition_size;
        uint64_t m_curr_val;
        uint64_t m_docid_lower_bound;
        uint64_t m_num_lists_before;

        void update_curr_val() { m_curr_val = m_docid_lower_bound + *m_curr_partition_it; }

        uint8_t msb(uint64_t n) const { return 64 - __builtin_clzll(n); }
    };

    typedef forward_iterator iterator_type;

    forward_iterator color_set(uint64_t color_set_id) const {
        assert(color_set_id < num_color_sets());
        uint64_t begin_partition_set =
            m_partition_sets_offsets.access(m_partition_sets_partitions.rank(color_set_id));
        uint64_t begin_rel = m_relative_colors_offsets.access(color_set_id);
        return forward_iterator(this, begin_partition_set, begin_rel);
    }

    uint32_t num_colors() const { return m_num_colors; }

    /* num. meta color lists */
    uint64_t num_color_sets() const { return m_relative_colors_offsets.size() - 1; }

    /* num. partial color sets */
    uint64_t num_partitions() const { return m_partition_endpoints.size(); }

    uint64_t num_partition_sets() const { return m_num_partition_sets; }

    uint64_t num_bits() const {
        uint64_t partial_colors_size = 0;
        for (auto d : m_partial_colors) { partial_colors_size += d.num_bits(); }
        return sizeof(size_t) * 8 + sizeof(m_num_colors) * 8 + sizeof(m_num_partition_sets) * 8 +
               m_partition_sets_offsets.num_bits() + m_relative_colors_offsets.num_bits() +
               essentials::vec_bytes(m_partition_endpoints) * 8 + partial_colors_size +
               essentials::vec_bytes(m_relative_colors) * 8 +
               essentials::vec_bytes(m_partition_sets) * 8 +
               m_partition_sets_partitions.bytes() * 8;
    }

    void print_stats() const {
        std::cout << "Color statistics:\n";
        std::cout << "  Number of partitions: " << num_partitions() << '\n';
        std::cout << "  Number of partition sets: " << num_partition_sets() << '\n';
        uint64_t num_bits_colors = 0;

        uint64_t meta_colors_size =
            (m_relative_colors_offsets.num_bits() + m_partition_sets_offsets.num_bits()) / 8 +
            essentials::vec_bytes(m_relative_colors) + essentials::vec_bytes(m_partition_sets) +
            m_partition_sets_partitions.size();

        for (auto const& c : m_partial_colors) { num_bits_colors += c.num_bits(); }

        assert(num_bits() > 0);
        std::cout << "  partial colors: " << num_bits_colors / 8 << " bytes ("
                  << (num_bits_colors * 100.0) / num_bits() << "%)\n";
        std::cout << "  meta colors: " << meta_colors_size << " bytes ("
                  << (meta_colors_size * 8 * 100.0) / num_bits() << "%)\n";
        std::cout << "  other: " << essentials::vec_bytes(m_partition_endpoints) << " bytes ("
                  << ((essentials::vec_bytes(m_partition_endpoints) * 8) * 100.0) / num_bits()
                  << "%)\n";
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
        visitor.visit(t.m_num_partition_sets);
        visitor.visit(t.m_partition_sets_offsets);
        visitor.visit(t.m_relative_colors_offsets);
        visitor.visit(t.m_partition_endpoints);
        visitor.visit(t.m_partial_colors);
        visitor.visit(t.m_relative_colors);
        visitor.visit(t.m_partition_sets);
        visitor.visit(t.m_partition_sets_partitions);
    }

    uint32_t m_num_colors;
    uint32_t m_num_partition_sets;
    sshash::ef_sequence<false> m_partition_sets_offsets, m_relative_colors_offsets;
    std::vector<partition_endpoint> m_partition_endpoints;
    std::vector<differential> m_partial_colors;
    std::vector<uint64_t> m_relative_colors;
    std::vector<uint64_t> m_partition_sets;
    ranked_bit_vector m_partition_sets_partitions;
};

}  // namespace fulgor

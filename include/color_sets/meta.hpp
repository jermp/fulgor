#pragma once

namespace fulgor {

template <typename ColorSets>
struct meta {
    static const bool meta_colored = true;
    static const bool differential_colored = false;

    struct partition_endpoint {
        template <typename Visitor>
        void visit(Visitor& visitor) {
            visitor.visit(docid_lower_bound);
            visitor.visit(num_lists_before);
        }
        uint32_t docid_lower_bound;
        uint32_t num_lists_before;
    };

    struct builder {
        builder() : m_offset(0) { m_meta_colors_offsets.push_back(0); }

        void init_meta_colors_builder(uint64_t num_integers_in_metacolors, uint64_t num_color_lists,
                                      std::vector<uint32_t> const& partition_sizes,
                                      std::vector<uint32_t> const& num_lists_in_partitions) {
            m_meta_colors_builder.resize(num_integers_in_metacolors,
                                         std::ceil(std::log2(num_color_lists)));
            m_partition_endpoints.reserve(num_lists_in_partitions.size());
            assert(partition_sizes.front() == 0);
            m_partition_endpoints.push_back({partition_sizes[0], 0});
            for (uint32_t i = 0, val = 0; i != num_lists_in_partitions.size(); ++i) {
                val += num_lists_in_partitions[i];
                m_partition_endpoints.push_back({partition_sizes[i + 1], val});
            }
        }

        void init_colors_builder(uint64_t num_colors, uint64_t num_partitions) {
            m_num_colors = num_colors;
            m_colors_builders.resize(num_partitions);
        }

        void init_color_partition(uint64_t partition_id, uint64_t num_colors_in_partition) {
            assert(partition_id < m_colors_builders.size());
            m_colors_builders[partition_id].init(num_colors_in_partition);
        }

        void process_colors(uint64_t partition_id, uint32_t const* colors, uint64_t list_size) {
            assert(partition_id < m_colors_builders.size());
            m_colors_builders[partition_id].process(colors, list_size);
        }

        void process_metacolors(uint32_t const* metacolors, uint64_t list_size) {
            assert(list_size < (1ULL << m_meta_colors_builder.width()));
            m_meta_colors_builder.push_back(list_size);
            for (uint64_t i = 0; i != list_size; ++i) {
                m_meta_colors_builder.push_back(metacolors[i]);
            }
            m_offset += list_size + 1;
            m_meta_colors_offsets.push_back(m_offset);
        }

        void build(meta& m) {
            m.m_num_colors = m_num_colors;
            m_meta_colors_builder.build(m.m_meta_colors);
            m.m_colors.resize(m_colors_builders.size());
            for (uint64_t i = 0; i != m_colors_builders.size(); ++i) {
                m_colors_builders[i].build(m.m_colors[i]);
            }
            m.m_meta_colors_offsets.encode(m_meta_colors_offsets.begin(),
                                           m_meta_colors_offsets.size(),
                                           m_meta_colors_offsets.back());
            m.m_partition_endpoints.swap(m_partition_endpoints);
        }

    private:
        pthash::compact_vector::builder m_meta_colors_builder;
        std::vector<typename ColorSets::builder> m_colors_builders;

        uint64_t m_num_colors;
        uint64_t m_offset;
        std::vector<uint64_t> m_meta_colors_offsets;

        std::vector<partition_endpoint> m_partition_endpoints;
    };

    struct forward_iterator {
        forward_iterator(meta<ColorSets> const* ptr, uint64_t begin)
            : m_ptr(ptr), m_begin(begin), m_meta_color_list_size((m_ptr->m_meta_colors)[m_begin]) {
            rewind();
        }

        void rewind() {
            init();
            assert(m_meta_color_list_size > 0);
            change_partition();
        }

        void init() {
            m_pos_in_meta_color_list = 0;
            m_partition_id = 0;
            m_partition_lower_bound = 0;
        }

        uint64_t value() const { return m_curr_val; }
        uint64_t operator*() const { return value(); }

        bool has_next() const { return m_pos_in_curr_partition != m_curr_partition_size; }
        void next_in_partition() {
            m_pos_in_curr_partition += 1;
            m_curr_partition_it.next();
            update_curr_val();
        }

        void next() {
            if (m_pos_in_curr_partition == m_curr_partition_size - 1) {
                if (m_pos_in_meta_color_list == meta_color_list_size() - 1) {  // saturate
                    m_curr_val = num_colors();
                    return;
                }
                m_pos_in_meta_color_list += 1;
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

        /* Warning: this might be slow. */
        uint32_t size() const {
            uint64_t n = 0;
            for (uint32_t i = 0, partition_id = 0; i != meta_color_list_size(); ++i) {
                uint32_t meta_color = (m_ptr->m_meta_colors)[m_begin + 1 + i];
                partition_id = update_partition_id(meta_color, partition_id);
                uint32_t num_lists_before =
                    (m_ptr->m_partition_endpoints)[partition_id].num_lists_before;
                n +=
                    (m_ptr->m_colors)[partition_id].color_set(meta_color - num_lists_before).size();
            }
            return n;
        }

        uint32_t meta_color() const { return m_curr_meta_color; }

        void read_partition_id() {
            m_curr_meta_color = (m_ptr->m_meta_colors)[m_begin + 1 + m_pos_in_meta_color_list];
            m_partition_id = update_partition_id(m_curr_meta_color, m_partition_id);
        }

        void next_partition_id() {
            m_pos_in_meta_color_list += 1;
            if (m_pos_in_meta_color_list == meta_color_list_size()) {  // saturate
                m_partition_id = num_partitions();
                return;
            }
            read_partition_id();
        }

        void next_geq_partition_id(const uint32_t lower_bound) {
            assert(lower_bound <= num_partitions());
            while (partition_id() < lower_bound) next_partition_id();
            assert(partition_id() >= lower_bound);
        }

        void update_partition() {
            /* update partition lower/upper bound */
            auto const& endpoints = m_ptr->m_partition_endpoints;
            m_partition_lower_bound = endpoints[m_partition_id].docid_lower_bound;
            m_partition_upper_bound = endpoints[m_partition_id + 1].docid_lower_bound;

            uint32_t num_lists_before = endpoints[m_partition_id].num_lists_before;
            m_curr_partition_it =
                (m_ptr->m_colors)[m_partition_id].color_set(m_curr_meta_color - num_lists_before);
            m_curr_partition_size = m_curr_partition_it.size();
            assert(m_curr_partition_size > 0);
            m_pos_in_curr_partition = 0;

            update_curr_val();
        }

        void change_partition() {
            read_partition_id();
            update_partition();
        }

        uint32_t partition_id() const { return m_partition_id; }
        uint32_t meta_color_list_size() const { return m_meta_color_list_size; }
        uint32_t num_colors() const { return m_ptr->num_colors(); }
        uint32_t num_partitions() const { return m_ptr->num_partitions(); }
        uint32_t partition_lower_bound() const { return m_partition_lower_bound; }
        uint32_t partition_upper_bound() const { return m_partition_upper_bound; }
        uint32_t num_lists_before() const {
            return m_ptr->m_partition_endpoints[m_partition_id].num_lists_before;
        }

    private:
        meta<ColorSets> const* m_ptr;
        typename ColorSets::iterator_type m_curr_partition_it;
        uint64_t m_begin;
        uint32_t m_curr_meta_color, m_curr_val;
        uint32_t m_meta_color_list_size, m_pos_in_meta_color_list;
        uint32_t m_curr_partition_size, m_pos_in_curr_partition;
        uint32_t m_partition_id;
        uint32_t m_partition_lower_bound, m_partition_upper_bound;

        void update_curr_val() {
            m_curr_val = m_curr_partition_it.value() + m_partition_lower_bound;
        }

        uint32_t update_partition_id(const uint32_t meta_color, uint32_t partition_id) const {
            auto const& endpoints = m_ptr->m_partition_endpoints;
            while (partition_id + 1 < endpoints.size() and
                   meta_color >= endpoints[partition_id + 1].num_lists_before) {
                partition_id += 1;
            }
            assert(partition_id < m_ptr->num_partitions());
            return partition_id;
        }
    };

    typedef forward_iterator iterator_type;

    forward_iterator color_set(uint64_t color_set_id) const {
        assert(color_set_id < num_color_sets());
        uint64_t begin = m_meta_colors_offsets.access(color_set_id);
        return forward_iterator(this, begin);
    }

    std::vector<ColorSets> partial_colors() const { return m_colors; }

    uint32_t num_colors() const { return m_num_colors; }

    /* num. meta color lists */
    uint64_t num_color_sets() const { return m_meta_colors_offsets.size() - 1; }

    /* num. partial color sets */
    uint64_t num_partitions() const { return m_partition_endpoints.size() - 1; }

    uint64_t num_bits() const {
        uint64_t num_bits_colors = sizeof(size_t) * 8;  // for std::vector::size
        for (auto const& c : m_colors) num_bits_colors += c.num_bits();
        return m_meta_colors_offsets.num_bits() + num_bits_colors +
               (m_meta_colors.bytes() + essentials::vec_bytes(m_partition_endpoints) +
                sizeof(m_num_colors)) *
                   8;
    }

    void print_stats() const {
        std::cout << "Color statistics:\n";
        std::cout << "  Number of partitions: " << num_partitions() << '\n';
        uint64_t num_bits_colors = 0;

        uint64_t num_partial_colors_very_dense = 0;
        uint64_t num_partial_colors_dense = 0;
        uint64_t num_partial_colors_sparse = 0;
        uint64_t num_total_partial_colors = 0;

        for (auto const& c : m_colors) {
            // c.print_stats();

            uint64_t n = c.num_color_sets();
            num_total_partial_colors += n;
            for (uint64_t i = 0; i != n; ++i) {
                auto it = c.color_set(i);
                if (it.type() == list_type::complement_delta_gaps) {
                    ++num_partial_colors_very_dense;
                } else if (it.type() == list_type::bitmap) {
                    ++num_partial_colors_dense;
                } else {
                    assert(it.type() == list_type::delta_gaps);
                    ++num_partial_colors_sparse;
                }
            }

            num_bits_colors += c.num_bits();
        }

        assert(num_total_partial_colors > 0);
        assert(num_bits() > 0);

        std::cout << "  num_partial_colors_very_dense = " << num_partial_colors_very_dense << " / "
                  << num_total_partial_colors << " ("
                  << (num_partial_colors_very_dense * 100.0) / num_total_partial_colors << "%)"
                  << std::endl;
        std::cout << "  num_partial_colors_dense = " << num_partial_colors_dense << " / "
                  << num_total_partial_colors << " ("
                  << (num_partial_colors_dense * 100.0) / num_total_partial_colors << "%)"
                  << std::endl;
        std::cout << "  num_partial_colors_sparse = " << num_partial_colors_sparse << " / "
                  << num_total_partial_colors << " ("
                  << (num_partial_colors_sparse * 100.0) / num_total_partial_colors << "%)"
                  << std::endl;

        std::cout << "  partial colors: " << num_bits_colors / 8 << " bytes ("
                  << (num_bits_colors * 100.0) / num_bits() << "%)\n";
        std::cout << "  meta colors: "
                  << m_meta_colors.bytes() + m_meta_colors_offsets.num_bits() / 8 << " bytes ("
                  << ((m_meta_colors.bytes() * 8 + m_meta_colors_offsets.num_bits()) * 100.0) /
                         num_bits()
                  << "%)\n";
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
        visitor.visit(t.m_meta_colors);
        visitor.visit(t.m_meta_colors_offsets);
        visitor.visit(t.m_colors);
        visitor.visit(t.m_partition_endpoints);
    }

    uint32_t m_num_colors;
    pthash::compact_vector m_meta_colors;
    sshash::ef_sequence<false> m_meta_colors_offsets;
    std::vector<ColorSets> m_colors;
    std::vector<partition_endpoint> m_partition_endpoints;
};

}  // namespace fulgor

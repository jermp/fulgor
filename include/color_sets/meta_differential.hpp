#pragma once

namespace fulgor {

struct meta_differential {
    static const index_t type = index_t::META_DIFF;

    struct partition_endpoint {
        template <typename Visitor>
        void visit(Visitor& visitor) {
            visitor.visit(min_color);
            visitor.visit(num_color_sets);
        }
        uint64_t min_color;
        uint64_t num_color_sets;
    };

    struct builder {
        builder() : m_prev_docs(0), m_prev_partition_set_id(0) {
            m_partition_sets_offsets.push_back(0);
            m_relative_colors_offsets.push_back(0);
        }

        void init(uint64_t num_colors, uint64_t num_partitions) {
            m_num_colors = num_colors;
            m_partition_endpoints.reserve(num_partitions);
            m_partial_color_sets.reserve(num_partitions);
        }

        void init_meta_color_partition_sets(uint64_t num_sets) {
            m_num_partition_sets = num_sets;
            m_partition_sets_offsets.reserve(num_sets);
        }

        void process_meta_color_partition_set(vector<uint64_t>& partition_set) {
            uint64_t size = partition_set.size();
            uint64_t prev_val = partition_set[0];

            bits::util::write_delta(m_partition_sets, size);
            bits::util::write_delta(m_partition_sets, prev_val);
            for (uint64_t i = 1; i < size; i++) {
                assert(prev_val < partition_set[i]);
                bits::util::write_delta(m_partition_sets, partition_set[i] - prev_val);
                prev_val = partition_set[i];
            }

            m_partition_sets_offsets.push_back(m_partition_sets.num_bits());
        }

        void process_partition(differential& d) {
            m_partial_color_sets.push_back(d);
            m_partition_endpoints.push_back({m_prev_docs, d.num_color_sets()});
            m_prev_docs += d.num_colors();
        }

        void process_metacolor_set(uint64_t partition_set_id,          //
                                   vector<uint64_t>& partition_set,    //
                                   vector<uint64_t>& relative_colors)  //
        {
            assert(partition_set.size() == relative_colors.size());
            if (partition_set_id != m_prev_partition_set_id) {
                m_prev_partition_set_id = partition_set_id;
                m_partition_sets_partitions.set(m_partition_sets_partitions.num_bits() - 1);
            }
            m_partition_sets_partitions.push_back(false);

            uint64_t size = partition_set.size();
            for (uint64_t i = 0; i < size; i++) {
                uint64_t partition_id = partition_set[i];
                uint64_t relative_id = relative_colors[i];
                uint64_t partition_size = m_partition_endpoints[partition_id].num_color_sets;
                m_relative_colors.append_bits(relative_id, bits::util::msbll(partition_size) + 1);
            }
            m_relative_colors_offsets.push_back(m_relative_colors.num_bits());
        }

        void build(meta_differential& m) {
            m.m_num_colors = m_num_colors;
            m.m_num_partition_sets = m_num_partition_sets;
            m.m_partial_color_sets.swap(m_partial_color_sets);

            m_relative_colors.build(m.m_relative_colors);
            m_partition_sets.build(m.m_partition_sets);
            m_partition_sets_partitions.build(m.m_partition_sets_partitions);
            m.m_partition_sets_partitions_rank1_index.build(m.m_partition_sets_partitions);

            m.m_partition_sets_offsets.encode(m_partition_sets_offsets.begin(),
                                              m_partition_sets_offsets.size(),
                                              m_partition_sets_offsets.back());
            m.m_relative_colors_offsets.encode(m_relative_colors_offsets.begin(),
                                               m_relative_colors_offsets.size(),
                                               m_relative_colors_offsets.back());

            m.m_partition_endpoints.swap(m_partition_endpoints);
        }

    private:
        std::vector<differential> m_partial_color_sets;

        bits::bit_vector::builder m_relative_colors;
        bits::bit_vector::builder m_partition_sets;
        bits::bit_vector::builder m_partition_sets_partitions;

        std::vector<uint64_t> m_partition_sets_offsets;
        std::vector<uint64_t> m_relative_colors_offsets;

        uint64_t m_num_colors, m_prev_docs;
        uint64_t m_num_partition_sets;

        uint64_t m_prev_partition_set_id;

        std::vector<partition_endpoint> m_partition_endpoints;
    };

    struct forward_iterator {
        forward_iterator(meta_differential const* ptr,  //
                         uint64_t begin_partition_set, uint64_t begin_rel)
            : m_ptr(ptr)
            , m_begin_partition_set(begin_partition_set)
            , m_begin_rel(begin_rel)  //
        {
            rewind();
            assert(m_meta_color_set_size > 0);
        }

        void rewind() {
            init();
            change_partition();
        }

        void init() {
            m_num_color_sets_before = 0;
            m_pos_in_meta_color = 0;
            m_pos_in_partial_color = 0;
            m_curr_partition_id = 0;
            m_partition_set_id = (m_ptr->m_partition_sets).get_iterator_at(m_begin_partition_set);
            m_relative_colors_it = (m_ptr->m_relative_colors).get_iterator_at(m_begin_rel);
            m_meta_color_set_size = bits::util::read_delta(m_partition_set_id);
        }

        uint64_t value() const { return m_curr_val; }
        uint64_t operator*() const { return value(); }

        bool has_next() const { return m_pos_in_partial_color != m_curr_partition_size; }

        void next() {
            if (m_pos_in_partial_color == m_curr_partition_size - 1) {
                if (m_pos_in_meta_color == meta_color_set_size() - 1) {  // saturate
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
            if (m_pos_in_meta_color == meta_color_set_size()) {
                m_curr_partition_id = num_partitions();
                return;
            }
            read_partition_id();
        }

        void read_partition_id() {
            uint64_t delta = bits::util::read_delta(m_partition_set_id);
            for (uint64_t i = 0; i < delta; i++) {
                m_num_color_sets_before +=
                    m_ptr->m_partition_endpoints[m_curr_partition_id + i].num_color_sets;
            }
            m_curr_partition_id += delta;
            uint8_t relative_color_size =
                bits::util::msbll(
                    m_ptr->m_partition_endpoints[m_curr_partition_id].num_color_sets) +
                1;
            m_curr_relative_color = m_relative_colors_it.take(relative_color_size);
        }

        void next_geq_partition_id(const uint32_t min_color) {
            assert(min_color <= num_partitions());
            while (partition_id() < min_color) next_partition_id();
            assert(partition_id() >= min_color);
        }

        void update_partition() {
            m_partition_min_color = m_ptr->m_partition_endpoints[m_curr_partition_id].min_color;

            m_pos_in_partial_color = 0;
            m_curr_partition_it =
                m_ptr->m_partial_color_sets[m_curr_partition_id].color_set(m_curr_relative_color);
            m_curr_partition_size = m_curr_partition_it.size();

            update_curr_val();
        }

        uint64_t size() const {
            uint64_t size = 0;
            auto partition_set_it =  //
                (m_ptr->m_partition_sets).get_iterator_at(m_begin_partition_set);
            auto rel_it =  //
                (m_ptr->m_relative_colors).get_iterator_at(m_begin_rel);
            uint64_t partition_id = 0;
            bits::util::read_delta(partition_set_it);
            for (uint64_t i = 0; i != m_meta_color_set_size; ++i) {
                partition_id += bits::util::read_delta(partition_set_it);
                uint8_t relative_color_size =
                    bits::util::msbll(m_ptr->m_partition_endpoints[partition_id].num_color_sets) +
                    1;
                uint64_t relative_color = rel_it.take(relative_color_size);
                size += m_ptr->m_partial_color_sets[partition_id].color_set(relative_color).size();
            }
            return size;
        }

        uint32_t partial_set_size() const { return m_curr_partition_it.size(); }

        uint32_t partition_id() const { return m_curr_partition_id; }
        uint32_t partition_min_color() const { return m_partition_min_color; }
        uint32_t partition_max_color() const {
            return m_partition_min_color + m_curr_partition_it.num_colors();
        }
        uint32_t meta_color() const { return m_num_color_sets_before + m_curr_relative_color; }
        uint32_t num_colors() const { return m_ptr->num_colors(); }
        uint32_t num_partitions() const { return m_ptr->num_partitions(); }
        uint64_t meta_color_set_size() const { return m_meta_color_set_size; }

        differential::iterator_type partition_it() const { return m_curr_partition_it; }

    private:
        meta_differential const* m_ptr;

        differential::iterator_type m_curr_partition_it;
        bits::bit_vector::iterator m_partition_set_id, m_relative_colors_it;

        uint64_t m_meta_color_set_size;
        uint64_t m_begin_partition_set, m_begin_rel;
        uint64_t m_pos_in_meta_color, m_pos_in_partial_color;
        uint64_t m_curr_relative_color;
        uint64_t m_curr_partition_id, m_curr_partition_size;
        uint64_t m_curr_val;
        uint64_t m_partition_min_color;
        uint64_t m_num_color_sets_before;

        void update_curr_val() { m_curr_val = m_partition_min_color + *m_curr_partition_it; }
    };

    typedef forward_iterator iterator_type;

    forward_iterator color_set(uint64_t color_set_id) const {
        assert(color_set_id < num_color_sets());
        uint64_t begin_partition_set =
            m_partition_sets_offsets.access(m_partition_sets_partitions_rank1_index.rank1(
                m_partition_sets_partitions, color_set_id));
        uint64_t begin_rel = m_relative_colors_offsets.access(color_set_id);
        return forward_iterator(this, begin_partition_set, begin_rel);
    }

    uint32_t num_colors() const { return m_num_colors; }
    uint64_t num_color_sets() const { return m_relative_colors_offsets.size() - 1; }
    uint64_t num_partitions() const { return m_partition_endpoints.size(); }
    uint64_t num_partition_sets() const { return m_num_partition_sets; }

    uint64_t num_bits() const {
        uint64_t num_bits_partial_color_sets = 0;
        for (auto const& c : m_partial_color_sets) num_bits_partial_color_sets += c.num_bits();
        return num_bits_partial_color_sets +  //
               (sizeof(size_t) + sizeof(m_num_colors) + sizeof(m_num_partition_sets) +
                m_relative_colors_offsets.num_bytes() + m_partition_sets_offsets.num_bytes() +
                essentials::vec_bytes(m_partition_endpoints) + m_relative_colors.num_bytes() +
                m_partition_sets.num_bytes() + m_partition_sets_partitions.num_bytes() +
                m_partition_sets_partitions_rank1_index.num_bytes()) *
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
        visitor.visit(t.m_num_partition_sets);
        visitor.visit(t.m_partition_sets_offsets);
        visitor.visit(t.m_relative_colors_offsets);
        visitor.visit(t.m_partition_endpoints);
        visitor.visit(t.m_partial_color_sets);
        visitor.visit(t.m_relative_colors);
        visitor.visit(t.m_partition_sets);
        visitor.visit(t.m_partition_sets_partitions);
        visitor.visit(t.m_partition_sets_partitions_rank1_index);
    }

    uint32_t m_num_colors;
    uint32_t m_num_partition_sets;
    bits::elias_fano<false, false> m_partition_sets_offsets, m_relative_colors_offsets;
    std::vector<partition_endpoint> m_partition_endpoints;
    std::vector<differential> m_partial_color_sets;
    bits::bit_vector m_relative_colors;
    bits::bit_vector m_partition_sets;
    bits::bit_vector m_partition_sets_partitions;
    bits::rank9 m_partition_sets_partitions_rank1_index;
};

}  // namespace fulgor

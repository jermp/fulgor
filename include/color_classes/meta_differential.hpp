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
        builder() : m_prev_docs(0), m_prev_base_id(0) { m_bases_offsets.push_back(0); m_relative_colors_offsets.push_back(0); }

        void init(uint64_t num_docs, uint64_t num_partitions) {
            m_num_docs = num_docs;
            m_partition_endpoints.reserve(num_partitions);
            m_partial_colors.reserve(num_partitions);
        }

        void init_meta_color_bases(uint64_t num_bases){
            m_num_bases = num_bases;
            m_bases_offsets.reserve(num_bases);
        }

        void process_meta_color_base(vector<uint64_t>& base){
            uint64_t size = base.size();
            uint64_t prev_val = base[0];

            util::write_delta(m_bases, size);
            util::write_delta(m_bases, prev_val);
            for(uint64_t i = 1; i < size; i++){
                assert(prev_val < base[i]);
                util::write_delta(m_bases, base[i] - prev_val);
                prev_val = base[i];
            }

            m_bases_offsets.push_back(m_bases.num_bits());
        }

        void process_partition(differential& d){
            m_partial_colors.push_back(d);
            m_partition_endpoints.push_back({m_prev_docs, d.num_color_classes()});
            m_prev_docs += d.num_docs();
        }

        void process_metacolors(uint64_t base_id, vector<uint64_t>& base, vector<uint64_t>& relative_colors){
            assert(base.size() == relative_colors.size());
            if (base_id != m_prev_base_id) {
                m_prev_base_id = base_id;
                m_bases_partitions.set(m_bases_partitions.size() - 1);
            }
            m_bases_partitions.push_back(false);

            uint64_t size = base.size();
            for (uint64_t i = 0; i < size; i++){
                uint64_t partition_id = base[i];
                uint64_t relative_id = relative_colors[i];

                uint64_t base_size = m_partition_endpoints[partition_id].num_lists;
                m_relative_colors.append_bits(relative_id, msb(base_size));
            }
            m_relative_colors_offsets.push_back(m_relative_colors.num_bits());
        }

        void build(meta_differential& m) {
            m.m_num_docs = m_num_docs;
            m.m_num_bases = m_num_bases;
            m.m_partial_colors.swap(m_partial_colors);
            m.m_relative_colors.swap(m_relative_colors.bits());
            m.m_bases.swap(m_bases.bits());

            m.m_bases_partitions.build(&m_bases_partitions);
            m.m_bases_offsets.encode(m_bases_offsets.begin(),
                                     m_bases_offsets.size(),
                                     m_bases_offsets.back());
            m.m_relative_colors_offsets.encode(m_relative_colors_offsets.begin(),
                                               m_relative_colors_offsets.size(),
                                               m_relative_colors_offsets.back());
            m.m_partition_endpoints.swap(m_partition_endpoints);
        }

    private:
        uint8_t msb(uint64_t n){
            return 64 - __builtin_clzll(n);
        } 

        std::vector<differential> m_partial_colors;
        bit_vector_builder m_relative_colors;
        bit_vector_builder m_bases;

        pthash::bit_vector_builder m_bases_partitions;

        std::vector<uint64_t> m_bases_offsets;
        std::vector<uint64_t> m_relative_colors_offsets;

        uint64_t m_num_docs, m_prev_docs;
        uint64_t m_num_bases;

        uint64_t m_prev_base_id;

        std::vector<partition_endpoint> m_partition_endpoints;
    };

    struct forward_iterator {
        forward_iterator(meta_differential const* ptr, uint64_t begin_base, uint64_t begin_rel)
            : m_ptr(ptr), m_begin_base(begin_base), m_begin_rel(begin_rel) {
            rewind();
            assert(m_meta_color_list_size > 0);
        }

        void rewind() {
            init();
            change_partition();
        }

        void init() {
            m_pos_in_meta_color = 0;
            m_pos_in_partial_color = 0;
            m_curr_partition_id = 0;
            m_base_it = bit_vector_iterator((m_ptr->m_bases).data(), (m_ptr->m_bases).size(), m_begin_base);
            m_relative_colors_it = bit_vector_iterator((m_ptr->m_relative_colors).data(), (m_ptr->m_relative_colors).size(),
                                                 m_begin_rel);
            m_meta_color_list_size = util::read_delta(m_base_it);
        }

        uint64_t value() const { return m_curr_val; }
        uint64_t operator*() const { return value(); }

        bool has_next() const { return m_pos_in_partial_color != m_curr_partition_size; }
 
        void next() {
            if (m_pos_in_partial_color == m_curr_partition_size - 1) {
                if (m_pos_in_meta_color == meta_color_list_size() - 1) {  // saturate
                    m_curr_val = num_docs();
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
            assert(lower_bound <= num_docs());
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
            if (m_pos_in_meta_color == meta_color_list_size()){
                m_curr_partition_id = num_partitions();
                return;
            }
            read_partition_id();
        }

        void read_partition_id() {
            m_curr_partition_id += util::read_delta(m_base_it);
            uint8_t relative_color_size = msb(m_ptr->m_partition_endpoints[m_curr_partition_id].num_lists);
            m_curr_relative_color = m_relative_colors_it.take(relative_color_size);
        }

        void next_geq_partition_id(const uint32_t lower_bound) {
            assert(lower_bound <= num_partitions());
            while (partition_id() < lower_bound) next_partition_id();
            assert(partition_id() >= lower_bound);
        }

        void update_partition() {
            m_docid_lower_bound = m_ptr->m_partition_endpoints[m_curr_partition_id].docid_lower_bound;
            
            m_pos_in_partial_color = 0;
            m_curr_partition_it = m_ptr->m_partial_colors[m_curr_partition_id].colors(m_curr_relative_color);
            m_curr_partition_size = m_curr_partition_it.size();
            
            update_curr_val();
        }

        uint64_t size() const { 
            uint64_t size = 0;
            auto base_it = bit_vector_iterator((m_ptr->m_bases).data(), (m_ptr->m_bases).size(), m_begin_base);
            auto rel_it = bit_vector_iterator((m_ptr->m_relative_colors).data(), (m_ptr->m_relative_colors).size(),
                                                 m_begin_rel);
            uint64_t partition_id = 0;
            util::read_delta(base_it); // remove size
            for (uint64_t partial_color_id = 0; partial_color_id < m_meta_color_list_size; partial_color_id++){
                partition_id += util::read_delta(base_it);
                uint8_t relative_color_size = msb(m_ptr->m_partition_endpoints[partition_id].num_lists);
                uint64_t relative_color = rel_it.take(relative_color_size);
                size += m_ptr->m_partial_colors[partition_id].colors(relative_color).size();
            }
            return size; 
        }

        uint32_t partition_id() const { return m_curr_partition_id; }
        uint32_t partition_upper_bound() const { return m_docid_lower_bound + m_curr_partition_it.num_docs(); }
        uint32_t meta_color() const {
            uint32_t meta_color = 0;
            for(uint64_t i = 0; i < m_curr_partition_id; i++){
                meta_color += m_ptr->m_partition_endpoints[i].num_lists;
            }
            return meta_color + m_curr_relative_color;
        }
        uint32_t num_docs() const { return m_ptr->num_docs(); }
        uint32_t num_partitions() const { return m_ptr->num_partitions(); }
        uint64_t meta_color_list_size() const { return m_meta_color_list_size; }

    private:
        meta_differential const* m_ptr;
        differential::iterator_type m_curr_partition_it;
        bit_vector_iterator m_base_it, m_relative_colors_it;
        uint64_t m_meta_color_list_size;
        uint64_t m_begin_base, m_begin_rel;
        uint64_t m_pos_in_meta_color, m_pos_in_partial_color;
        uint64_t m_curr_relative_color;
        uint64_t m_curr_partition_id, m_curr_partition_size;
        uint64_t m_curr_val;
        uint64_t m_docid_lower_bound;

        void update_curr_val(){
            m_curr_val = m_docid_lower_bound + *m_curr_partition_it;
        }

        uint8_t msb(uint64_t n) const{
            return 64 - __builtin_clzll(n);
        } 
    };

    typedef forward_iterator iterator_type;

    forward_iterator colors(uint64_t color_class_id) const {
        assert(color_class_id < num_color_classes());
        uint64_t begin_base = m_bases_offsets.access(m_bases_partitions.rank(color_class_id));
        uint64_t begin_rel = m_relative_colors_offsets.access(color_class_id);
        return forward_iterator(this, begin_base, begin_rel);
    }

    uint32_t num_docs() const { return m_num_docs; }

    /* num. meta color lists */
    uint64_t num_color_classes() const { return m_relative_colors_offsets.size() - 1; }

    /* num. partial color sets */
    uint64_t num_partitions() const { return m_partition_endpoints.size(); }

    uint64_t num_bases() const { return m_num_bases; }

    uint64_t num_bits() const {
        uint64_t partial_colors_size = 0;
        for (auto d: m_partial_colors){
            partial_colors_size += d.num_bits();
        }
        return sizeof(m_num_docs) * 8 +
            sizeof(m_num_bases) * 8 + 
            m_bases_offsets.num_bits() +
            m_relative_colors_offsets.num_bits() + 
            essentials::vec_bytes(m_partition_endpoints) * 8 +
            partial_colors_size +
            essentials::vec_bytes(m_relative_colors) * 8 +
            essentials::vec_bytes(m_bases) * 8 +
            m_bases_partitions.size();
    }

    void print_stats() const {
        std::cout << "Color statistics:\n";
        std::cout << "  Number of partitions: " << num_partitions() << '\n';
        std::cout << "  Number of bases: " << num_bases() << '\n';
        uint64_t colors_bits = 0;

        uint64_t meta_colors_size = (m_relative_colors_offsets.num_bits() + m_bases_offsets.num_bits()) / 8 +
            essentials::vec_bytes(m_relative_colors) + 
            essentials::vec_bytes(m_bases) +
            m_bases_partitions.size();

        for (auto const& c : m_partial_colors) {
            colors_bits += c.num_bits();
        }

        std::cout << "  partial colors: " << colors_bits / 8 << " bytes ("
                  << (colors_bits * 100.0) / num_bits() << "%)\n";
        std::cout << "  meta colors: "
                  << meta_colors_size << " bytes ("
                  << (meta_colors_size * 8 * 100.0) / num_bits()
                  << "%)\n";
        std::cout << "  other: " << essentials::vec_bytes(m_partition_endpoints) << " bytes ("
                  << ((essentials::vec_bytes(m_partition_endpoints) * 8) * 100.0) / num_bits()
                  << "%)\n";
    }

    void dump(std::ofstream& os) const {

    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_num_docs);
        visitor.visit(m_num_bases);
        visitor.visit(m_bases_offsets);
        visitor.visit(m_relative_colors_offsets);
        visitor.visit(m_partition_endpoints);
        visitor.visit(m_partial_colors);
        visitor.visit(m_relative_colors);
        visitor.visit(m_bases);
        visitor.visit(m_bases_partitions);
    }

    
private:
    uint32_t m_num_docs;
    uint32_t m_num_bases; 
    sshash::ef_sequence<false> m_bases_offsets, m_relative_colors_offsets;
    std::vector<partition_endpoint> m_partition_endpoints;
    std::vector<differential> m_partial_colors; 
    std::vector<uint64_t> m_relative_colors; 
    std::vector<uint64_t> m_bases; 
    ranked_bit_vector m_bases_partitions;
};

} // namespace fulgor
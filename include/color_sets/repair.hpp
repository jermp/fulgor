#pragma once

#include "../builders/repair_builder.hpp"

namespace fulgor{

struct repair{
	static const index_t type = index_t::REPAIR;

    struct builder{
        builder() {}
        builder(uint64_t num_colors): m_num_colors(num_colors) {}

        void init_code(std::vector<code_type> C){
            m_C_builder.resize(C.size(), std::ceil(std::log2(m_D.size())));
            m_C_builder.fill(C.begin(), C.size());
        }

        void init_dict(std::vector<std::vector<uint32_t>> D){
            m_D = D;
        }

        void init_sizes(std::vector<uint32_t>& sizes, uint32_t max_size){
            uint64_t n = sizes.size();
            max_size = max_size > m_C_builder.size() ? max_size : m_C_builder.size();
            m_sizes_builder.resize(n, std::ceil(std::log2(max_size)));
            m_sizes_builder.fill(sizes.begin(), n);
            m_sizes = sizes;
        }

        void build(repair& r){
            r.m_num_colors = m_num_colors;
            m_C_builder.build(r.m_C);
            m_sizes_builder.build(r.m_sizes);

            uint32_t n = m_sizes.size();
            vector<uint32_t> offsets;
            offsets.reserve(n);
            uint32_t max_offsets_val = 0;

            uint32_t pos = 0;
            auto it = r.m_C.begin();
            uint32_t offset = 0;
            uint32_t val = *it;

            for(uint32_t i = 0; i < n-1; i++){
                // cout << i << "-" << val << "  " << flush;
                offsets.push_back(pos);
                offsets.push_back(offset);
                max_offsets_val = max(max_offsets_val, pos);
                max_offsets_val = max(max_offsets_val, offset);

                uint32_t target = m_sizes[i] + offset;
                while (m_D[val].size() <= target){
                    target -= m_D[val].size();
                    val = *it;
                    ++it;
                    ++pos;
                }
                offset = target;
            }
            
            m_offsets_builder.resize(2*n, std::ceil(std::log2(max_offsets_val)));
            m_offsets_builder.fill(offsets.begin(), offsets.size());

            m_offsets_builder.build(r.m_offsets);
            r.m_D.swap(m_D);
        }

private:
            uint64_t m_num_colors;
            std::vector<uint32_t> m_sizes;
            pthash::compact_vector::builder m_C_builder, m_sizes_builder, m_offsets_builder;
            std::vector<std::vector<uint32_t>> m_D;
            // std::vector<uint32_t> m_begins;
    };
    
    struct forward_iterator{
        forward_iterator() {}

        forward_iterator(repair const* ptr, uint32_t begin, uint32_t offset, uint32_t size)
            : m_ptr(ptr)
            , m_begin(begin)
            , m_offset(offset)
            , m_size(size) {
            rewind();    
        }

        void rewind(){
            init();
            update_curr_val();
        }

        uint64_t value() const { return m_curr_val; }
        uint64_t operator*() const { return value(); }

        void next() {
            if (m_curr_pos >= m_size) return;

            m_curr_pos++;
            m_curr_offset++;
            uint32_t c = (m_ptr->m_C)[m_curr_symbol];
            if (m_curr_offset >= (m_ptr->m_D)[c].size()){
                m_curr_symbol++;
                m_curr_offset = 0;
            }
            update_curr_val();
        }
        void operator++() { next(); }

        void next_geq(uint32_t lower_bound) {
            while (value() < lower_bound) next();
        }

        uint64_t size() const {
            return m_size;
        }

        uint32_t num_colors() const {
            return m_ptr->num_colors();
        }

    private:
        repair const* m_ptr;
        uint32_t m_begin, m_offset;
        uint32_t m_curr_symbol, m_curr_offset, m_curr_pos;
        uint32_t m_curr_val;
        uint32_t m_size;

        void init(){
            m_curr_pos = 0;
            m_curr_symbol = m_begin;
            m_curr_offset = m_offset;
            m_curr_val = -1;
        }

        void update_curr_val(){
            if (m_curr_pos >= m_size){
                m_curr_val = num_colors();
                return;
            }
            uint32_t c = (m_ptr->m_C)[m_curr_symbol];
            m_curr_val += (m_ptr->m_D)[c][m_curr_offset] + 1;
        }
    };

    typedef forward_iterator iterator_type;

    iterator_type color_set(uint64_t color_set_id) const{
        return forward_iterator(this, 
                m_offsets[color_set_id * 2], 
                m_offsets[color_set_id * 2 + 1], 
                m_sizes[color_set_id]
            );
    }

    uint64_t num_color_sets() const { return m_sizes.size(); } //remove first (universe) and -1
    uint64_t num_colors() const { return m_num_colors; }

    uint64_t num_bits() const {
        return sizeof(m_num_colors)*8 + (m_C.bytes() + m_offsets.bytes() + essentials::vec_bytes(m_D)) * 8; 
    }

    void print_stats() const {
        std::cout << "    C: " << m_C.bytes() << " bytes" << std::endl;
        std::cout << "    D: " << essentials::vec_bytes(m_D) << " bytes" << std::endl;
        std::cout << "    dict size: " << m_D.size() << " items" << std::endl;
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
        visitor.visit(t.m_C);
        visitor.visit(t.m_D);
        visitor.visit(t.m_offsets);
        visitor.visit(t.m_sizes);
    }

    uint64_t m_num_colors;
    pthash::compact_vector m_C;
    std::vector<std::vector<uint32_t>> m_D;
    pthash::compact_vector m_offsets;
    pthash::compact_vector m_sizes;
   //  sshash::ef_sequence<false> m_begins;
};

}


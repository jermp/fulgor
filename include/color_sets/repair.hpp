#pragma once

#include "../builders/repair_builder.hpp"

namespace fulgor{

struct repair{
	static const index_t type = index_t::REPAIR;

    struct builder{
        builder() {}
        builder(uint64_t num_colors): m_num_colors(num_colors) {}

        void set_code(std::vector<code_type> C){
            m_C_builder.resize(C.size(), std::ceil(std::log2(m_D.size())));
            m_C_builder.fill(C.begin(), C.size());
        }

        void set_dict(std::vector<std::vector<uint32_t>> D){
            m_D = D;
        }

        void build(repair& r){
            r.m_num_colors = m_num_colors;
            r.m_D.swap(m_D);
            m_C_builder.build(r.m_C);
        }

        private:
            uint64_t m_num_colors;
            pthash::compact_vector::builder m_C_builder;
            std::vector<std::vector<uint32_t>> m_D;
    };
    
    struct forward_iterator{
        forward_iterator() {}

        uint64_t size(){
            return 0;
        }
    };

    typedef forward_iterator iterator_type;

    iterator_type color_set(uint64_t color_set_id) const{
        return forward_iterator();
    }

    uint64_t num_color_sets() const { return 0; } //remove first (universe) and -1
    uint64_t num_colors() const { return 0; }

    uint64_t num_bits() const {
        return sizeof(m_num_colors)*8 + (m_C.bytes() + essentials::vec_bytes(m_D)) * 8; 
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
    }

    uint64_t m_num_colors;
    pthash::compact_vector m_C;
    std::vector<std::vector<uint32_t>> m_D;
};

}


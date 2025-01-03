#pragma once

#include "../builders/repair_builder.hpp"

namespace fulgor{

struct repair{
	static const index_t type = index_t::REPAIR;

    struct builder{
        builder() {}
        builder(uint64_t num_colors): m_num_colors(num_colors) {}

        void set_code(std::vector<code_type> C){
            m_C = C;
        }

        void set_dict(std::vector<std::vector<uint32_t>> D){
            m_D = D;
        }

        void build(repair& r){
            r.m_num_colors = m_num_colors;
            r.m_D.swap(m_D);
            r.m_C.swap(m_C);
        }

        private:
            uint64_t m_num_colors;
            std::vector<code_type> m_C;
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
        return sizeof(m_num_colors)*8 + (essentials::vec_bytes(m_C) + essentials::vec_bytes(m_D)) * 8; 
    }

    void print_stats() const {
        // TODO: implement
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
    std::vector<code_type> m_C;
    std::vector<std::vector<uint32_t>> m_D;
};

}


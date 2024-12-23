#pragma once

#include "external/s_indexes/include/s_index.hpp"

namespace fulgor{

struct sliced{
	static const index_t type = index_t::SLICED;

    struct builder{
        builder() {}
        builder(uint64_t num_colors): m_num_colors(num_colors) {}

        void push_offsets(std::vector<uint64_t> offsets){
            m_offsets = offsets;
        }

        void push_sequences(std::vector<uint8_t> sequences){
            m_sequences = sequences;
        }

        void build(sliced& s){
            s.m_num_colors = m_num_colors;
            s.m_offsets.swap(m_offsets);
            s.m_sequences.swap(m_sequences);
        }

        private:
            uint64_t m_num_colors;
            std::vector<uint64_t> m_offsets;
            std::vector<uint8_t> m_sequences;
    };
    
    struct forward_iterator{
        forward_iterator() {}

        forward_iterator(sliced const* ptr, uint64_t color_set_id): m_ptr(ptr), m_color_set_id(color_set_id) {
            rewind();
        }

        void rewind() {
            auto sequence = ::sliced::s_sequence((m_ptr->m_sequences).data() + (m_ptr->m_offsets)[m_color_set_id+1]);
            m_size = sequence.size();
            m_it = sequence.begin();
        }

        uint64_t value() const {
           return m_it.id(); 
        }
        uint64_t operator*() const { return value(); }

        void next() { m_it.next(); }
        void operator++() { next(); }

        void next_geq(const uint64_t lower_bound) {
            m_it.advance(lower_bound);
        }

        uint32_t num_colors() const { return m_ptr->m_num_colors; }

        uint64_t size() const {
            return m_size;
        }

        ::sliced::s_sequence sequence() {
            return ::sliced::s_sequence((m_ptr->m_sequences).data() + (m_ptr->m_offsets)[m_color_set_id+1]);
        }

        private:
            sliced const* m_ptr;
            uint64_t m_color_set_id;

            uint64_t m_size;
            ::sliced::s_sequence::iterator m_it;
    };

    typedef forward_iterator iterator_type;

    iterator_type color_set(uint64_t color_set_id) const{
        return forward_iterator(this, color_set_id);
    }

    uint64_t num_color_sets() const { return m_offsets.size() - 2; } //remove first (universe) and -1
    uint64_t num_colors() const { return m_num_colors; }

    uint64_t num_bits() const {
        return sizeof(m_num_colors) * 8 + (essentials::vec_bytes(m_offsets) + essentials::vec_bytes(m_sequences)) * 8;
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
        visitor.visit(t.m_offsets);
        visitor.visit(t.m_sequences);
    }

    uint64_t m_num_colors;
    std::vector<uint64_t> m_offsets;
    std::vector<uint8_t> m_sequences;
};

}

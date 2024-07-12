#pragma once

#include <vector>

#include "external/sshash/external/pthash/external/essentials/include/essentials.hpp"
#include "external/sshash/external/pthash/include/encoders/bit_vector.hpp"

namespace fulgor {

struct ranked_bit_vector : public pthash::bit_vector {
    ranked_bit_vector() : pthash::bit_vector() {}

    void build(pthash::bit_vector_builder* bvb) {
        pthash::bit_vector::build(bvb);
        build_index();
    }

    inline uint64_t num_ones() const { return *(m_block_rank_pairs.end() - 2); }
    inline uint64_t num_zeros() const { return size() - num_ones(); }

    /* return the number of ones in A[0..pos) */
    inline uint64_t rank(uint64_t pos) const {
        assert(pos <= size());
        if (pos == size()) return num_ones();
        uint64_t sub_block = pos / 64;
        uint64_t r = sub_block_rank(sub_block);
        uint64_t sub_left = pos % 64;
        if (sub_left) r += pthash::util::popcount(m_bits[sub_block] << (64 - sub_left));
        return r;
    }

    uint64_t bytes() const {
        return pthash::bit_vector::bytes() + essentials::vec_bytes(m_block_rank_pairs);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        pthash::bit_vector::visit(visitor);
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        pthash::bit_vector::visit(visitor);
        visit_impl(visitor, *this);
    }

protected:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_block_rank_pairs);
    }

    inline uint64_t block_rank(uint64_t block) const { return m_block_rank_pairs[block * 2]; }

    inline uint64_t sub_block_rank(uint64_t sub_block) const {
        uint64_t r = 0;
        uint64_t block = sub_block / block_size;
        r += block_rank(block);
        uint64_t left = sub_block % block_size;
        r += sub_block_ranks(block) >> ((7 - left) * 9) & 0x1FF;
        return r;
    }

    inline uint64_t sub_block_ranks(uint64_t block) const {
        return m_block_rank_pairs[block * 2 + 1];
    }

    void build_index() {
        std::vector<uint64_t> block_rank_pairs;
        uint64_t next_rank = 0;
        uint64_t cur_subrank = 0;
        uint64_t subranks = 0;
        block_rank_pairs.push_back(0);
        for (uint64_t i = 0; i < m_bits.size(); ++i) {
            uint64_t word_pop = pthash::util::popcount(m_bits[i]);
            uint64_t shift = i % block_size;
            if (shift) {
                subranks <<= 9;
                subranks |= cur_subrank;
            }
            next_rank += word_pop;
            cur_subrank += word_pop;

            if (shift == block_size - 1) {
                block_rank_pairs.push_back(subranks);
                block_rank_pairs.push_back(next_rank);
                subranks = 0;
                cur_subrank = 0;
            }
        }
        uint64_t left = block_size - m_bits.size() % block_size;
        for (uint64_t i = 0; i < left; ++i) {
            subranks <<= 9;
            subranks |= cur_subrank;
        }
        block_rank_pairs.push_back(subranks);

        if (m_bits.size() % block_size) {
            block_rank_pairs.push_back(next_rank);
            block_rank_pairs.push_back(0);
        }

        m_block_rank_pairs.swap(block_rank_pairs);
    }

    static const uint64_t block_size = 8;  // in 64bit words
    std::vector<uint64_t> m_block_rank_pairs;
};

}  // namespace fulgor
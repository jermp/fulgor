#pragma once

#include "bit_vector.hpp"
#include "integer_codes.hpp"

namespace fulgor {

struct inverted_index {
    struct builder {
        builder(uint32_t num_docs, std::string const& filename) : m_num_written_bytes(0) {
            m_out.open(filename.c_str(), std::ofstream::binary);
            if (!m_out.is_open()) throw std::runtime_error("cannot open file");
            m_out.write(reinterpret_cast<char const*>(&num_docs), sizeof(num_docs));

            /* reserved space: to be updated with finalize() */
            uint64_t num_ints = 0;
            m_out.write(reinterpret_cast<char const*>(&num_ints), sizeof(num_ints));

            m_num_written_bytes += sizeof(num_docs) + sizeof(num_ints);
        }

        void finalize(uint64_t num_ints) {
            /* pad the file size to a multiple of 8 */
            uint64_t mod = m_num_written_bytes % 8;
            if (mod) {
                uint64_t num_bytes_pad = 8 - mod;
                assert((num_bytes_pad + m_num_written_bytes) % 8 == 0);
                uint64_t pad = 0;
                m_out.write(reinterpret_cast<char const*>(&pad), num_bytes_pad);
            }

            /* after the 4 bytes fo num_docs, over-write num_ints */
            m_out.seekp(4);
            m_out.write(reinterpret_cast<char const*>(&num_ints), sizeof(num_ints));
            m_out.close();
        }

        template <typename Iterator>
        void write_list(seg_id_t seg_id, Iterator list_begin, uint64_t list_size) {
            m_bvb.clear();

            /* reserve a large amounts of bits, corresponding to uncompressed list */
            m_bvb.reserve(list_size * 32);

            /* write seg_id and list_size */
            m_bvb.append_bits(seg_id, sizeof(seg_id_t) * 8);
            util::write_delta(m_bvb, list_size);

            /* encode first value */
            uint32_t prev_val = *list_begin;
            util::write_delta(m_bvb, prev_val);
            ++list_begin;

            /* encode the gaps */
            for (uint64_t i = 1; i != list_size; ++i, ++list_begin) {
                uint32_t val = *list_begin;
                assert(val >= prev_val + 1);
                util::write_delta(m_bvb, val - prev_val - 1);
                prev_val = val;
            }

            /* save the the compressed representation padded to an integral number of bytes*/
            uint64_t num_bytes = (m_bvb.num_bits() + 8 - 1) / 8;
            m_out.write(reinterpret_cast<char const*>(m_bvb.data()), num_bytes);
            m_num_written_bytes += num_bytes;
        }

    private:
        std::ofstream m_out;
        bit_vector_builder m_bvb;
        uint64_t m_num_written_bytes;
    };

    /* forward iterator */
    template <bool with_cache = false>
    struct iterator {
        iterator() : m_num_docs(0), m_num_ints(0), m_num_64bit_words(0) {}

        iterator(uint64_t const* data, uint64_t num_64bit_words) : m_it(data, num_64bit_words) {
            m_num_docs = m_it.take(32);
            m_num_ints = m_it.take(64);
            m_num_64bit_words = num_64bit_words;
            m_curr_list = inv_list(&m_it);
        }

        uint32_t num_docs() const { return m_num_docs; }
        uint64_t num_ints() const { return m_num_ints; }
        bool has_next() const { return m_it.position() != 64 * m_num_64bit_words; }

        struct inv_list {
            inv_list() : m_it(nullptr) {}
            inv_list(bit_vector_iterator* it) : m_it(it) {
                assert(m_it->position() % 8 == 0);
                m_seg_id = m_it->take(sizeof(seg_id_t) * 8);
                m_it->fill_buf();
                m_size = util::read_delta(*m_it);
                if constexpr (with_cache) decode();
            }

            uint64_t seg_id() const { return m_seg_id; }
            uint32_t size() const { return m_size; }
            std::vector<uint32_t> const& cache() const { return m_cache; }

            struct iterator  // forward_iterator
            {
                iterator(bit_vector_iterator* it, uint32_t pos_in_list, uint32_t size)
                    : m_it(it)
                    , m_pos_in_list(pos_in_list)
                    , m_size(size)
                    , m_prev_val(-1)
                    , m_curr_val(0) {
                    if (pos_in_list != m_size) m_curr_val = util::read_delta(*m_it);
                }

                bool operator==(iterator const& rhs) const {
                    return m_pos_in_list == rhs.m_pos_in_list;
                }
                bool operator!=(iterator const& rhs) const { return !(*this == rhs); }
                uint32_t operator*() const { return m_curr_val; }
                void operator++() {
                    m_pos_in_list += 1;
                    if (m_pos_in_list >= m_size) return;
                    m_prev_val = m_curr_val;
                    m_curr_val = util::read_delta(*m_it) + (m_prev_val + 1);
                }

            private:
                bit_vector_iterator* m_it;
                uint32_t m_pos_in_list;
                uint32_t m_size;
                uint32_t m_prev_val;
                uint32_t m_curr_val;
            };

            iterator begin() const { return iterator(m_it, 0, size()); }
            iterator end() const { return iterator(m_it, size(), size()); }

        private:
            bit_vector_iterator* m_it;
            seg_id_t m_seg_id;
            uint32_t m_size;
            std::vector<uint32_t> m_cache;

            void decode() {
                m_cache.reserve(m_size);
                std::copy(begin(), end(), std::back_inserter(m_cache));
            }
        };

        inv_list& list() { return m_curr_list; }

        void advance_to_next_list() {
            uint64_t mod = m_it.position() % 8;
            uint64_t pad = 0;
            if (mod) {
                /* eat pad for list */
                pad = 8 - mod;
                m_it.take(pad);
            }
            /* eat pad for whole index */
            if (64 * m_num_64bit_words - m_it.position() < 64) {
                pad = 64 * m_num_64bit_words - m_it.position();
                m_it.take(pad);
            }
            if (has_next()) m_curr_list = inv_list(&m_it);
        }

    private:
        uint32_t m_num_docs;
        uint64_t m_num_ints;
        uint64_t m_num_64bit_words;
        bit_vector_iterator m_it;
        inv_list m_curr_list;
    };
};

}  // namespace fulgor

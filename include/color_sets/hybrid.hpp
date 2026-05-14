#pragma once

#include "include/util.hpp"
#include <shared_mutex>

namespace fulgor {

struct hybrid {
    static constexpr index_t type = HYBRID;

    struct builder {
        explicit builder() : m_max_RAM_bytes(1<<30), m_verbose(false) {
            init(0, "color_sets.bin");
        }

        explicit builder(const uint64_t num_colors, const build_configuration& build_config)
            : m_max_RAM_bytes(build_config.ram_limit_in_GiB << 30)
            , m_verbose(build_config.ram_limit_in_GiB) {
            init(num_colors, build_config.tmp_dirname + "/color_sets.bin");
        }

        builder(const uint32_t num_colors, const std::string& output_filename,
                         const uint64_t max_RAM_bytes = 1<<30, const bool verbose = false)
            : m_output_filename(output_filename)
            , m_max_RAM_bytes(max_RAM_bytes)
            , m_verbose(verbose)
        {
            init(num_colors, output_filename);
        }

        void init(const uint64_t num_colors, const std::string& output_filename) {
            m_num_colors = num_colors;
            m_sparse_set_threshold_size = 0.25 * m_num_colors;
            m_very_dense_set_threshold_size = 0.75 * m_num_colors;

            m_offsets.clear();
            m_offsets.push_back(0);
            m_num_color_sets = 0;
            m_num_total_integers = 0;

            m_curr_color_set_id = 0;
            m_queue_size = 0;
            m_base_offset = 0;

            m_output_filename = output_filename;
            std::ofstream file(m_output_filename.c_str());
            constexpr uint64_t zero = 0;
            file.write(reinterpret_cast<const char*>(&zero), sizeof(zero));
            file.write(reinterpret_cast<const char*>(&zero), sizeof(zero));
        }

        void set_max_RAM_bytes(const uint64_t max_RAM_bytes) {
            m_max_RAM_bytes = max_RAM_bytes;
        }
        void set_verbose(const bool verbose) {
            m_verbose = verbose;
        }

        void reserve_num_bits(uint64_t num_bits) { m_color_sets_builder.reserve(num_bits); }

        void encode_color_set(std::span<const uint32_t> color_set, const uint64_t color_set_id)  //
        {
            if (size() >= m_max_RAM_bytes) {
                std::unique_lock lock(*m_flush_mutex);
                if (size() >= m_max_RAM_bytes) {
                    flush();
                }
            }

            const uint64_t cs_size = color_set.size();
            bits::bit_vector::builder bvb;
            bits::util::write_delta(bvb, cs_size);
            if (cs_size < m_sparse_set_threshold_size) {
                encode_sparse(color_set, bvb);
            } else if (cs_size < m_very_dense_set_threshold_size) {
                encode_dense(color_set, bvb);
            } else {
                encode_very_dense(color_set, bvb);
            }

            {
                std::shared_lock flush_lock(*m_flush_mutex);
                std::lock_guard queue_lock(*m_queue_mutex);
                m_num_total_integers += cs_size;
                m_num_color_sets += 1;
                if (color_set_id == m_curr_color_set_id) {
                    m_color_sets_builder.append(bvb);
                    m_offsets.push_back(m_base_offset + m_color_sets_builder.num_bits());
                    ++m_curr_color_set_id;
                    while (!m_csb_queue.empty() && m_curr_color_set_id == m_csb_queue.top().first) {
                        auto& [_, csb] = m_csb_queue.top();
                        m_queue_size -= 8 + essentials::vec_bytes(csb.data()) + 16; // first + bvb struct

                        m_color_sets_builder.append(csb);
                        m_csb_queue.pop();
                        m_offsets.push_back(m_base_offset + m_color_sets_builder.num_bits());
                        ++m_curr_color_set_id;
                    }
                    assert((m_csb_queue.empty() && m_queue_size == 0) || (!m_csb_queue.empty() && m_queue_size > 0));
                } else {
                    m_queue_size += 8 + essentials::vec_bytes(bvb.data()) + 16; // first + bvb struct
                    m_csb_queue.emplace(color_set_id, std::move(bvb));
                }
            }

            if (m_verbose && m_num_color_sets % 500000 == 0) {
                std::cout << "  processed " << m_num_color_sets << " color sets" << std::endl;
            }
        }

        uint64_t size() const {
            return essentials::vec_bytes(m_color_sets_builder.data()) +
                essentials::vec_bytes(m_offsets) + m_queue_size;
        }

        void append(hybrid::builder& hb) {
            if (hb.m_num_color_sets == 0) return;
            m_color_sets_builder.append(hb.m_color_sets_builder);
            assert(m_offsets.size() > 0);
            uint64_t delta = m_offsets.back();
            m_offsets.reserve(m_offsets.size() + hb.m_offsets.size());
            for (uint64_t i = 1; i != hb.m_offsets.size(); ++i) {
                m_offsets.push_back(hb.m_offsets[i] + delta);
            }
            m_num_color_sets += hb.m_num_color_sets;
            m_num_total_integers += hb.m_num_total_integers;
            assert(m_num_color_sets == m_offsets.size() - 1);
        }

        void build(hybrid& h) {
            h.m_num_colors = m_num_colors;
            h.m_sparse_set_threshold_size = m_sparse_set_threshold_size;
            h.m_very_dense_set_threshold_size = m_very_dense_set_threshold_size;

            std::cout << "processed " << m_num_color_sets << " color sets" << std::endl;
            std::cout << "m_num_total_integers " << m_num_total_integers << std::endl;
            assert(m_num_color_sets == m_offsets.size() - 1);

            h.m_offsets.encode(m_offsets.begin(), m_offsets.size(), m_offsets.back());
            // m_color_sets_builder.build(h.m_color_sets);
            essentials::load(h.m_color_sets, m_output_filename.c_str());

            std::cout << "  total bits for ints = " << 8 * h.m_color_sets.num_bytes() << std::endl;
            std::cout << "  total bits per offsets = " << 8 * h.m_offsets.num_bytes() << std::endl;
            std::cout << "  total bits = "
                      << 8 * (h.m_color_sets.num_bytes() + h.m_offsets.num_bytes()) << std::endl;
            std::cout << "  offsets: " << (8.0 * h.m_offsets.num_bytes()) / m_num_total_integers
                      << " bits/int" << std::endl;
            std::cout << "  color sets: "
                      << (8.0 * h.m_color_sets.num_bytes()) / m_num_total_integers << " bits/int"
                      << std::endl;
        }

        void clear() {
            m_offsets.clear();
            m_color_sets_builder.clear();
            init(m_num_colors, m_output_filename);
        }

        void flush() {
            if (m_color_sets_builder.data().empty()) return;
            std::fstream file(m_output_filename, std::ios::in | std::ios::out | std::ios::binary);
            if (!file.is_open()) return;

            uint64_t num_bits;
            file.read(reinterpret_cast<char*>(&num_bits), sizeof(uint64_t));

            const uint64_t pos = num_bits >> 6 ;
            const uint64_t final_num_bits = (pos << 6) + m_color_sets_builder.num_bits();
            const uint64_t final_num_words = (final_num_bits + 63) >> 6;
            m_base_offset = final_num_bits & ~63;
            assert(final_num_bits == m_offsets.back());

            file.seekp(sizeof(uint64_t)*2  + (pos * sizeof(uint64_t)));
            file.write(reinterpret_cast<char*>(
                m_color_sets_builder.data().data()),
                m_color_sets_builder.data().size() * sizeof(uint64_t)
            );

            file.seekp(0);
            file.write(reinterpret_cast<const char*>(&final_num_bits), sizeof(uint64_t));
            file.write(reinterpret_cast<const char*>(&final_num_words), sizeof(uint64_t));

            uint64_t last_word = m_color_sets_builder.data().back();
            m_color_sets_builder.clear();
            m_color_sets_builder.reserve(m_max_RAM_bytes << 3); //bits
            m_color_sets_builder.append_bits(last_word, final_num_bits & 63);
        }

    private:
        uint32_t m_num_colors;
        uint32_t m_sparse_set_threshold_size;
        uint32_t m_very_dense_set_threshold_size;
        uint64_t m_num_color_sets;
        uint64_t m_num_total_integers;

        uint64_t m_curr_color_set_id;
        std::unique_ptr<std::mutex> m_queue_mutex = std::make_unique<std::mutex>();
        std::unique_ptr<std::shared_mutex> m_flush_mutex = std::make_unique<std::shared_mutex>();
        bits::bit_vector::builder m_color_sets_builder;
        std::vector<uint64_t> m_offsets;
        uint64_t m_base_offset;
        std::priority_queue<
            std::pair<uint64_t, bits::bit_vector::builder>,
            std::vector<std::pair<uint64_t, bits::bit_vector::builder>>,
            util::compare_first
        > m_csb_queue;
        uint64_t m_queue_size;

        std::string m_output_filename;
        uint64_t m_max_RAM_bytes;
        bool m_verbose;

        void encode_sparse(const std::span<const uint32_t> color_set, bits::bit_vector::builder& bvb) const {
            const uint64_t size = color_set.size();
            uint32_t prev_val = color_set[0];
            bits::util::write_delta(bvb, prev_val);
            for (uint64_t i = 1; i != size; ++i) {
                const uint32_t val = color_set[i];
                assert(val >= prev_val + 1);
                bits::util::write_delta(bvb, val - (prev_val + 1));
                prev_val = val;
            }
        }

        void encode_dense(const std::span<const uint32_t> color_set, bits::bit_vector::builder& bvb) const {
            const uint64_t size = color_set.size();
            bits::bit_vector::builder tmp_bvb;
            tmp_bvb.resize(m_num_colors);
            for (uint64_t i = 0; i != size; ++i) tmp_bvb.set(color_set[i]);
            bvb.append(tmp_bvb);
        }

        void encode_very_dense(const std::span<const uint32_t> color_set, bits::bit_vector::builder& bvb) const {
            const uint64_t size = color_set.size();
            bool first = true;
            uint32_t val = 0;
            uint32_t prev_val = -1;
            uint32_t written = 0;
            for (uint64_t i = 0; i != size; ++i) {
                uint32_t x = color_set[i];
                while (val < x) {
                    if (first) {
                        bits::util::write_delta(bvb, val);
                        first = false;
                    } else {
                        assert(val >= prev_val + 1);
                        bits::util::write_delta(bvb, val - (prev_val + 1));
                    }
                    prev_val = val;
                    ++val;
                    ++written;
                }
                assert(val == x);
                val++;  // skip x
            }
            while (val < m_num_colors) {
                assert(val >= prev_val + 1);
                bits::util::write_delta(bvb, val - (prev_val + 1));
                prev_val = val;
                ++val;
                ++written;
            }
            assert(val == m_num_colors);
            /* complementary_set_size = m_num_colors - size */
            assert(m_num_colors - size <= m_num_colors);
            assert(written == m_num_colors - size);
            (void)written;  // silence "unused" warning
        }
    };

    struct forward_iterator {
        forward_iterator() {}

        forward_iterator(hybrid const* ptr, uint64_t begin)
            : m_ptr(ptr)
            , m_bitmap_begin(begin)
            , m_color_sets_begin(begin)
            , m_num_colors(ptr->m_num_colors) {
            rewind();
        }

        void rewind() {
            m_pos_in_set = 0;
            m_pos_in_comp_set = 0;
            m_comp_set_size = 0;
            m_comp_val = -1;
            m_prev_val = -1;
            m_curr_val = 0;
            m_it = (m_ptr->m_color_sets).get_iterator_at(m_color_sets_begin);
            m_size = bits::util::read_delta(m_it);
            /* set m_encoding_type and read the first value */
            if (m_size < m_ptr->m_sparse_set_threshold_size) {
                m_encoding_type = encoding_t::delta_gaps;
                m_curr_val = bits::util::read_delta(m_it);
            } else if (m_size < m_ptr->m_very_dense_set_threshold_size) {
                m_encoding_type = encoding_t::bitmap;
                m_bitmap_begin = m_it.position();  // after m_size
                m_it.skip_to(m_bitmap_begin);
                uint64_t pos = m_it.next();
                assert(pos >= m_bitmap_begin);
                m_curr_val = pos - m_bitmap_begin;
            } else {
                m_encoding_type = encoding_t::complement_delta_gaps;
                m_comp_set_size = m_num_colors - m_size;
                if (m_comp_set_size > 0) m_comp_val = bits::util::read_delta(m_it);
                next_comp_val();
            }
        }

        /* this is needed to annul the next_comp_val() done in the constructor
           if we want to iterate through the complemented set */
        void reinit_for_complemented_set_iteration() {
            assert(m_encoding_type == encoding_t::complement_delta_gaps);
            m_pos_in_comp_set = 0;
            m_prev_val = -1;
            m_curr_val = 0;
            m_it = (m_ptr->m_color_sets).get_iterator_at(m_color_sets_begin);
            bits::util::read_delta(m_it); /* skip m_size */
            if (m_comp_set_size > 0) {
                m_comp_val = bits::util::read_delta(m_it);
            } else {
                m_comp_val = m_num_colors;
            }
        }

        uint64_t value() const { return m_curr_val; }
        uint64_t comp_value() const { return m_comp_val; }
        uint64_t operator*() const { return value(); }

        void next() {
            if (m_encoding_type == encoding_t::complement_delta_gaps) {
                ++m_curr_val;
                if (m_curr_val >= m_num_colors) {  // saturate
                    m_curr_val = m_num_colors;
                    return;
                }
                next_comp_val();
            } else if (m_encoding_type == encoding_t::delta_gaps) {
                m_pos_in_set += 1;
                if (m_pos_in_set >= m_size) {  // saturate
                    m_curr_val = m_num_colors;
                    return;
                }
                m_prev_val = m_curr_val;
                m_curr_val = bits::util::read_delta(m_it) + (m_prev_val + 1);
            } else {
                assert(m_encoding_type == encoding_t::bitmap);
                m_pos_in_set += 1;
                if (m_pos_in_set >= m_size) {  // saturate
                    m_curr_val = m_num_colors;
                    return;
                }
                uint64_t pos = m_it.next();
                assert(pos >= m_bitmap_begin);
                m_curr_val = pos - m_bitmap_begin;
            }
        }

        void next_comp() {
            ++m_pos_in_comp_set;
            if (m_pos_in_comp_set >= m_comp_set_size) {  // saturate
                m_comp_val = m_num_colors;
                return;
            }
            m_prev_val = m_comp_val;
            m_comp_val = bits::util::read_delta(m_it) + (m_prev_val + 1);
        }

        void operator++() { next(); }

        /* update the state of the iterator to the element
           which is greater-than or equal-to lower_bound */
        void next_geq(const uint64_t lower_bound) {
            assert(lower_bound <= num_colors());
            if (m_encoding_type == encoding_t::complement_delta_gaps) {
                if (value() > lower_bound) return;
                next_geq_comp_val(lower_bound);
                m_curr_val = lower_bound + (m_comp_val == lower_bound);
            } else {
                while (value() < lower_bound) next();
            }
            assert(value() >= lower_bound);
        }

        uint32_t size() const { return m_size; }
        uint32_t num_colors() const { return m_num_colors; }
        int encoding_type() const { return m_encoding_type; }

    private:
        hybrid const* m_ptr;
        uint64_t m_bitmap_begin;
        uint64_t m_color_sets_begin;
        uint32_t m_num_colors;
        int m_encoding_type;

        bits::bit_vector::iterator m_it;
        uint32_t m_pos_in_set;
        uint32_t m_size;

        uint32_t m_pos_in_comp_set;
        uint32_t m_comp_set_size;

        uint32_t m_comp_val;
        uint32_t m_prev_val;
        uint32_t m_curr_val;

        void next_comp_val() {
            while (m_curr_val == m_comp_val) {
                ++m_curr_val;
                ++m_pos_in_comp_set;
                if (m_pos_in_comp_set >= m_comp_set_size) break;
                m_prev_val = m_comp_val;
                m_comp_val = bits::util::read_delta(m_it) + (m_prev_val + 1);
            }
        }

        void next_geq_comp_val(const uint64_t lower_bound) {
            while (m_comp_val < lower_bound) {
                ++m_pos_in_comp_set;
                if (m_pos_in_comp_set >= m_comp_set_size) break;
                m_prev_val = m_comp_val;
                m_comp_val = bits::util::read_delta(m_it) + (m_prev_val + 1);
            }
        }
    };

    typedef forward_iterator iterator_type;

    forward_iterator color_set(uint64_t color_set_id) const {
        assert(color_set_id < num_color_sets());
        uint64_t begin = m_offsets.access(color_set_id);
        return forward_iterator(this, begin);
    }

    uint32_t num_colors() const { return m_num_colors; }
    uint64_t num_color_sets() const { return m_offsets.size() - 1; }

    uint64_t num_bits() const {
        return (sizeof(m_num_colors) + sizeof(m_sparse_set_threshold_size) +
                sizeof(m_very_dense_set_threshold_size) + m_offsets.num_bytes() +
                m_color_sets.num_bytes()) *
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
        visitor.visit(t.m_sparse_set_threshold_size);
        visitor.visit(t.m_very_dense_set_threshold_size);
        visitor.visit(t.m_offsets);
        visitor.visit(t.m_color_sets);
    }

    uint32_t m_num_colors;
    uint32_t m_sparse_set_threshold_size;
    uint32_t m_very_dense_set_threshold_size;

    bits::bit_vector m_color_sets;
    bits::elias_fano<false, false> m_offsets;
};

}  // namespace fulgor

// 13959930901 218123921
//  5639917171  88123706
//  6561605869 102525092


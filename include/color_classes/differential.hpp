#pragma once

namespace fulgor {

struct differential {
    static const bool meta_colored = false;
    static const bool differential_colored = true;

    struct builder {
        builder() : m_prev_cluster_id(0) {
            // TODO: reserve bvb space
            m_list_offsets.push_back(0);
            m_reference_offsets.push_back(0);
        }

        void init_colors_builder(uint64_t num_docs) {
            m_num_docs = num_docs;
            m_num_total_integers = 0;
            m_num_lists = 0;
        }

        void encode_reference(std::vector<uint32_t> const& reference) {
            uint64_t size = reference.size();
            util::write_delta(m_bvb, size);
            m_num_total_integers += size + 1;  // size plus size number
            m_num_lists += 1;

            if (size > 0) {
                uint32_t prev_val = reference[0];
                util::write_delta(m_bvb, prev_val);

                for (uint64_t i = 1; i < size; ++i) {
                    uint32_t val = reference[i];
                    assert(val >= prev_val + 1);
                    util::write_delta(m_bvb, val - (prev_val + 1));
                    prev_val = val;
                }
            }
            m_reference_offsets.push_back(m_bvb.num_bits());
        }

        void encode_list(uint64_t cluster_id, std::vector<uint32_t> const& reference,
                         hybrid::forward_iterator it) {
            /*
                To be fixed later:
                avoid allocating memory here. We can avoid the vector
                std::vector<uint32_t> edit_list entirely: just write
                directly to m_bvb each element of the edit list
                once computed. (We need to keep track of the previously
                written value to take the gap.)
            */

            std::vector<uint32_t> edit_list;
            uint64_t ref_size = reference.size();
            uint64_t it_size = it.size();
            edit_list.reserve(ref_size + it_size);

            if (cluster_id != m_prev_cluster_id) {
                m_prev_cluster_id = cluster_id;
                m_clusters.set(m_clusters.size() - 1);
            }
            m_clusters.push_back(false);

            uint64_t i = 0, j = 0;
            while (i < it_size && j < ref_size) {
                if (*it == reference[j]) {
                    i += 1;
                    ++it;
                    j += 1;
                } else if (*it < reference[j]) {
                    edit_list.push_back(*it);
                    i += 1;
                    ++it;
                } else {
                    edit_list.push_back(reference[j]);
                    j += 1;
                }
            }
            while (i < it_size) {
                edit_list.push_back(*it);
                i += 1;
                ++it;
            }
            while (j < ref_size) {
                edit_list.push_back(reference[j]);
                j += 1;
            }

            uint64_t size = edit_list.size();
            util::write_delta(m_bvb, size);
            util::write_delta(m_bvb, it.size());
            m_num_total_integers += size + 2;  // size plus edit_list size plus original list size
            m_num_lists += 1;

            if (size > 0) {
                uint32_t prev_val = edit_list[0];
                util::write_delta(m_bvb, prev_val);

                for (uint64_t pos = 1; pos < size; ++pos) {
                    uint32_t val = edit_list[pos];
                    assert(val >= prev_val + 1);
                    util::write_delta(m_bvb, val - (prev_val + 1));
                    prev_val = val;
                }
            }

            uint64_t last_offset = m_reference_offsets[m_reference_offsets.size() - 1];
            m_list_offsets.push_back(m_bvb.num_bits() - last_offset);
        }

        void build(differential& d) {
            d.m_num_docs = m_num_docs;
            d.m_colors.swap(m_bvb.bits());
            d.m_clusters.build(&m_clusters);

            d.m_reference_offsets.encode(m_reference_offsets.begin(), m_reference_offsets.size(),
                                         m_reference_offsets.back());
            d.m_list_offsets.encode(m_list_offsets.begin(), m_list_offsets.size(),
                                    m_list_offsets.back());

            std::cout << "processed " << m_num_lists << " lists\n";
            std::cout << "m_num_total_integers " << m_num_total_integers << '\n';

            std::cout << "  total bits for ints = " << d.m_colors.size() * 64 << '\n';
            std::cout << "  total bits per offset = "
                      << d.m_list_offsets.num_bits() + d.m_reference_offsets.num_bits()
                      << " (lists: " << d.m_list_offsets.num_bits()
                      << ", references: " << d.m_reference_offsets.num_bits() << ")\n";
            std::cout << "  offsets: "
                      << static_cast<double>(d.m_list_offsets.num_bits() +
                                             d.m_reference_offsets.num_bits()) /
                             m_num_total_integers
                      << " bits/int\n";
            std::cout << "  lists: "
                      << static_cast<double>(d.m_colors.size() * 64) / m_num_total_integers
                      << " bits/int\n";
        }

    private:
        bit_vector_builder m_bvb;
        pthash::bit_vector_builder m_clusters;
        uint64_t m_num_total_integers, m_num_lists;

        uint64_t m_num_docs;
        uint64_t m_prev_cluster_id;
        std::vector<uint64_t> m_reference_offsets, m_list_offsets;
    };

    struct forward_iterator {
        forward_iterator(differential const* ptr, uint64_t list_begin, uint64_t reference_begin)
            : m_ptr(ptr), m_edit_list_begin(list_begin), m_reference_begin(reference_begin) {
            rewind();
        }

        void rewind() {
            m_edit_list_it = bit_vector_iterator((m_ptr->m_colors).data(), (m_ptr->m_colors).size(),
                                                 m_edit_list_begin);
            m_reference_it = bit_vector_iterator((m_ptr->m_colors).data(), (m_ptr->m_colors).size(),
                                                 m_reference_begin);
            m_edit_list_size = util::read_delta(m_edit_list_it);
            m_reference_size = util::read_delta(m_reference_it);
            m_size = util::read_delta(m_edit_list_it);

            m_curr_edit_val = m_edit_list_size == 0 ? num_docs() : util::read_delta(m_edit_list_it);
            m_prev_edit_val = 0;
            m_curr_reference_val =
                m_reference_size == 0 ? num_docs() : util::read_delta(m_reference_it);
            m_prev_reference_val = 0;

            m_pos_in_edit_list = 0;
            m_pos_in_reference = 0;
            update_curr_val();
        }

        uint32_t size() { return m_size; }

        uint64_t value() const { return m_curr_val; }
        uint64_t operator*() const { return value(); }

        void next() {
            if (m_pos_in_reference >= m_reference_size && m_pos_in_edit_list >= m_edit_list_size) {
                m_curr_val = num_docs();
                return;
            }
            if (m_pos_in_reference >= m_reference_size || m_curr_edit_val < m_curr_reference_val) {
                next_edit_val();
            } else if (m_pos_in_edit_list >= m_edit_list_size ||
                       m_curr_reference_val < m_curr_edit_val) {
                next_reference_val();
            }
            update_curr_val();
        }
        void operator++() { next(); }

        uint32_t num_docs() const { return m_ptr->m_num_docs; }

    private:
        differential const* m_ptr;
        uint64_t m_edit_list_begin, m_reference_begin;
        uint64_t m_reference_size, m_edit_list_size;
        uint64_t m_pos_in_edit_list, m_pos_in_reference;
        uint32_t m_curr_reference_val, m_curr_edit_val;
        uint32_t m_prev_reference_val, m_prev_edit_val;
        uint32_t m_curr_val;
        uint32_t m_size;
        bit_vector_iterator m_reference_it, m_edit_list_it;

        void next_reference_val() {
            m_pos_in_reference += 1;
            m_prev_reference_val = m_curr_reference_val;
            if (m_pos_in_reference < m_reference_size) {
                m_curr_reference_val = m_prev_reference_val + util::read_delta(m_reference_it) + 1;
            } else {
                m_curr_reference_val = num_docs();
            }
        }

        void next_edit_val() {
            m_pos_in_edit_list += 1;
            m_prev_edit_val = m_curr_edit_val;
            if (m_pos_in_edit_list < m_edit_list_size) {
                m_curr_edit_val = m_prev_edit_val + util::read_delta(m_edit_list_it) + 1;
            } else {
                m_curr_edit_val = num_docs();
            }
        }

        void update_curr_val() {
            while (m_curr_reference_val == m_curr_edit_val &&
                   m_pos_in_reference <= m_reference_size &&
                   m_pos_in_edit_list <= m_edit_list_size) {
                next_edit_val();
                next_reference_val();
            }
            m_curr_val = min(m_curr_edit_val, m_curr_reference_val);
        }
    };

    typedef forward_iterator iterator_type;

    forward_iterator colors(uint64_t color_id) const {
        assert(color_id < num_color_classes());
        uint64_t last_reference = m_reference_offsets.access(num_partitions());
        uint64_t list_begin = m_list_offsets.access(color_id) + last_reference;
        uint64_t reference_begin = m_reference_offsets.access(m_clusters.rank(color_id));
        return forward_iterator(this, list_begin, reference_begin);
    }

    uint64_t num_color_classes() const { return m_list_offsets.size() - 1; }
    uint64_t num_partitions() const { return m_clusters.num_ones() + 1; }
    uint64_t num_docs() const { return m_num_docs; }

    uint64_t num_bits() const {
        return sizeof(m_num_docs) * 8 + m_reference_offsets.num_bits() + m_list_offsets.num_bits() +
               essentials::vec_bytes(m_colors) * 8 + m_clusters.size();
    }

    void print_stats() const {
        std::cout << "Color statistics:\n";
        std::cout << "  Number of partitions: " << num_partitions() << std::endl;

        uint64_t num_reference_offsets = m_reference_offsets.num_bits();
        uint64_t num_list_offsets = m_list_offsets.num_bits();
        uint64_t num_colors = essentials::vec_bytes(m_colors) * 8;
        uint64_t num_clusters = m_clusters.size();

        uint64_t num_references = 0;
        uint64_t num_edit_lists = 0;
        uint64_t num_metadata = 0;

        uint64_t num_docs_tenth = num_docs()/10;

        std::vector<uint64_t> distribution(11, 0);
        std::vector<uint64_t> distribution_0_10(num_docs_tenth,0);
        std::vector<uint64_t> distribution_0_10_bits(num_docs_tenth,0);

        for (uint64_t reference_id = 0; reference_id < num_partitions(); reference_id++) {
            uint64_t reference_begin = m_reference_offsets.access(reference_id);
            auto it = bit_vector_iterator(m_colors.data(), m_colors.size(), reference_begin);
            uint64_t prev_position = it.position();

            uint64_t size = util::read_delta(it);
            num_metadata += it.position() - prev_position;
            prev_position = it.position();

            for (uint64_t i = 0; i < size; i++) {
                util::read_delta(it);
                num_references += it.position() - prev_position;
                prev_position = it.position();
            }
        }
        uint64_t last_reference = m_reference_offsets.access(num_partitions());
        for (uint64_t color_id = 0; color_id < num_color_classes(); color_id++) {
            uint64_t list_begin = m_list_offsets.access(color_id) + last_reference;
            auto it = bit_vector_iterator(m_colors.data(), m_colors.size(), list_begin);
            uint64_t prev_position = it.position();

            uint64_t size = util::read_delta(it);
            num_metadata += it.position() - prev_position;
            prev_position = it.position();

            util::read_delta(it); // original list size
            num_metadata += it.position() - prev_position;
            prev_position = it.position();

            uint64_t curr_edit_list_size = 0;
            for (uint64_t i = 0; i < size; i++) {
                util::read_delta(it);
                uint64_t delta_size = it.position() - prev_position;
                num_edit_lists += delta_size;
                curr_edit_list_size += delta_size;

                prev_position = it.position();
            }

            uint64_t q = size/(num_docs_tenth) > 10 ? 10 : size/(num_docs_tenth);
            if (q == 0){
                distribution_0_10[size]++;
                distribution_0_10_bits[size] += curr_edit_list_size;
            }

            distribution[q]++;
        }

        std::cout << "  reference offsets: " << num_reference_offsets / 8 << " bytes ("
                  << (num_reference_offsets * 100.0) / num_bits() << "%)" << std::endl;
        std::cout << "  edit list offsets: " << num_list_offsets / 8 << " bytes ("
                  << (num_list_offsets * 100.0) / num_bits() << "%)" << std::endl;
        std::cout << "  clusters: " << num_clusters / 8 << " bytes ("
                  << (num_clusters * 100.0) / num_bits() << "%)" << std::endl;
        std::cout << "  differential colors: " << num_colors / 8 << " bytes ("
                  << (num_colors * 100.0) / num_bits() << "%)" << std::endl;
        std::cout << "    references: " << num_references / 8 << " bytes ("
                  << (num_references * 100.0) / num_colors << "%)" << std::endl;
        std::cout << "    edit lists: " << num_edit_lists / 8 << " bytes ("
                  << (num_edit_lists * 100.0) / num_colors << "%)" << std::endl;
        std::cout << "    metadata: " << num_metadata / 8 << " bytes ("
                  << (num_metadata * 100.0) / num_colors << "%)" << std::endl;
        std::cout << "  edit lists size distribution:" << std::endl;
        for (uint64_t partition = 0; partition < 11; partition++){
            std::cout << "    range " << partition * num_docs_tenth << " -> " << (partition+1) * num_docs_tenth-1 << ": " << distribution[partition] << std::endl;
        }
        std::cout << "  edit lists size distribution, detail 0% - 10%:" << std::endl;
        for (uint64_t length = 0; length < distribution_0_10.size(); length++){
            std::cout << "    [" << length << "] num_edit_lists: " << distribution_0_10[length] << ", num_bits: "<< distribution_0_10_bits[length] <<
                " (" << distribution_0_10_bits[length]*100.0 / num_edit_lists << "%), bits/int: " << 1.*distribution_0_10_bits[length] / distribution_0_10[length] / length << std::endl;
        }

        std::cout << std::endl;
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_num_docs);
        visitor.visit(m_reference_offsets);
        visitor.visit(m_list_offsets);
        visitor.visit(m_colors);
        visitor.visit(m_clusters);
    }

private:
    uint32_t m_num_docs;

    sshash::ef_sequence<false> m_reference_offsets, m_list_offsets;

    std::vector<uint64_t> m_colors;
    ranked_bit_vector m_clusters;
};

}  // namespace fulgor
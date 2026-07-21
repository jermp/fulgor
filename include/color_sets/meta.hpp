#pragma once
#include <future>

namespace fulgor {

template <typename ColorSets>
struct meta {
    static constexpr index_t type = META;

    struct partition_endpoint {
        template <typename Visitor>
        void visit(Visitor& visitor) {
            visitor.visit(min_color);
            visitor.visit(num_color_sets_before);
        }
        uint32_t min_color;
        uint32_t num_color_sets_before;
    };

    struct builder {
        explicit builder(const uint64_t num_colors, util::external_saver& saver,
                         const std::vector<uint32_t>& partition_starts,
                         const std::string& tmp_dirname, const uint64_t max_RAM_bytes = 8,
                         const bool verbose = false)
            : m_offset(0)
            , m_tmp_dirname(tmp_dirname)
            , m_max_RAM_bytes(max_RAM_bytes)
            , m_verbose(verbose) {
            m_meta_color_sets_offsets.push_back(0);
            init(num_colors, saver, partition_starts);
        }

        void init(const uint64_t num_colors, util::external_saver& saver,
                  std::vector<uint32_t> const& partition_starts) {
            m_num_colors = num_colors;
            m_saver = &saver;

            const uint32_t num_partitions = partition_starts.size() - 1;
            m_partitions_mutex.reserve(num_partitions);
            m_partition_savers.reserve(num_partitions);
            m_hashes.resize(num_partitions);
            m_partition_starts = partition_starts;

            for (uint64_t partition_id = 0; partition_id < num_partitions; ++partition_id) {
                m_partitions_mutex.push_back(std::make_unique<std::shared_mutex>());

                m_partition_savers.emplace_back(partition_filename(partition_id));
                const uint32_t num_colors_in_partition =
                    partition_starts[partition_id + 1] - partition_starts[partition_id];

                m_color_sets_builders.emplace_back(num_colors_in_partition,
                                                   m_partition_savers.back());
                m_color_sets_builders[partition_id].set_verbose(m_verbose);
                m_color_sets_builders[partition_id].set_max_RAM_bytes(m_max_RAM_bytes);
            }
        }

        void init_meta_color_sets_builder(const uint64_t num_integers_in_metacolor_sets)  //
        {
            uint64_t num_partial_color_sets = 0;
            m_partial_sets_starts.push_back(0);
            for (auto& dict : m_hashes) {
                num_partial_color_sets += dict.size();
                m_partial_sets_starts.push_back(num_partial_color_sets);
            }
            m_meta_color_sets_builder.resize(num_integers_in_metacolor_sets,
                                             std::ceil(std::log2(num_partial_color_sets)));
            m_partition_endpoints.reserve(m_hashes.size() + 1);

            assert(m_partition_starts.size() == m_partial_sets_starts.size());
            for (uint32_t partition_id = 0; partition_id != m_partition_starts.size();
                 ++partition_id) {
                m_partition_endpoints.push_back(
                    {m_partition_starts[partition_id], m_partial_sets_starts[partition_id]});
            }
        }

        [[deprecated("External memory construction makes this useless")]]
        void reserve_num_bits(uint64_t partition_id, uint64_t num_bits) {
            assert(partition_id < m_color_sets_builders.size());
            m_color_sets_builders[partition_id].reserve_num_bits(num_bits);
        }

        uint32_t num_partitions() const { return m_color_sets_builders.size(); }

        std::vector<uint32_t> encode(std::vector<uint32_t>& color_set) {
            uint32_t partition_id = 0;
            uint32_t part_end = m_partition_starts[1];
            std::vector<uint32_t> metacolor_set;
            metacolor_set.reserve(num_partitions() * 2);

            const std::span data(color_set);
            uint32_t span_start = 0;
            for (uint64_t i = 0; i < data.size(); ++i) {
                const uint32_t color = data[i];
                while (color >= part_end) {
                    auto partition_view = data.subspan(span_start, i - span_start);
                    if (!partition_view.empty()) {
                        const uint32_t partial_color_set_id =
                            encode_partial_color_set(partition_id, partition_view);
                        metacolor_set.push_back(partition_id);
                        metacolor_set.push_back(partial_color_set_id);
                    }

                    ++partition_id;
                    part_end = m_partition_starts[partition_id + 1];
                    span_start = i;
                }
                assert(color >= m_partition_starts[partition_id]);
            }
            const auto final_view = data.subspan(span_start);
            if (!final_view.empty()) {
                const uint32_t partial_color_set_id =
                    encode_partial_color_set(partition_id, final_view);
                metacolor_set.push_back(partition_id);
                metacolor_set.push_back(partial_color_set_id);
            }
            return metacolor_set;
        }

        void encode_metacolor_set(std::span<std::pair<uint32_t, uint32_t>> metacolor_set) {
            const uint64_t size = metacolor_set.size();
            assert(size < (1ULL << m_meta_color_sets_builder.width()));
            m_meta_color_sets_builder.push_back(size);
            for (auto& [partition_id, partial_cs_id] : metacolor_set) {
                m_meta_color_sets_builder.push_back(partial_cs_id +
                                                    m_partial_sets_starts[partition_id]);
            }
            m_offset += size + 1;
            m_meta_color_sets_offsets.push_back(m_offset);
        }

        void flush() {
            if (m_is_flushing.exchange(true)) {
                return;
            }

            for (auto& builder : m_color_sets_builders) {
                builder.flush();
            }
            m_is_flushing = false;
        }

        uint64_t num_bytes() const {
            uint64_t size = 0;
            for (const auto& builder : m_color_sets_builders) {
                size += builder.size();
            }
            return size;
        }

        [[deprecated("External memory construction requires build()")]]
        void build(meta& m) {
            flush();
            m.m_num_colors = m_num_colors;
            m_meta_color_sets_builder.build(m.m_meta_color_sets);
            m.m_partial_color_sets.resize(m_color_sets_builders.size());
            for (uint64_t partition_id = 0; partition_id < m_color_sets_builders.size();
                 ++partition_id) {
                m_color_sets_builders[partition_id].build();
            }

            m.m_partial_color_sets.reserve(num_partitions());
            for (uint64_t partition_id = 0; partition_id != m_color_sets_builders.size();
                 ++partition_id) {
                essentials::load(m.m_partial_color_sets[partition_id],
                                 partition_filename(partition_id).c_str());
            }

            m.m_meta_color_sets_offsets.encode(m_meta_color_sets_offsets.begin(),
                                               m_meta_color_sets_offsets.size(),
                                               m_meta_color_sets_offsets.back());
            m.m_partition_endpoints.swap(m_partition_endpoints);
        }

        void build() {
            flush();
            assert(m_saver != nullptr);
            m_saver->visit(m_num_colors);

            {
                bits::compact_vector metacolor_sets;
                m_meta_color_sets_builder.build(metacolor_sets);
                m_saver->visit(metacolor_sets);
            }
            {
                bits::elias_fano offsets;
                offsets.encode(m_meta_color_sets_offsets.begin(), m_meta_color_sets_offsets.size(),
                               m_meta_color_sets_offsets.back());
                m_saver->visit(offsets);
            }

            const size_t num_parts = num_partitions();
            m_saver->visit(num_parts);
            for (uint64_t partition_id = 0; partition_id < num_parts; ++partition_id) {
                m_color_sets_builders[partition_id].build();
                m_partition_savers[partition_id].close();

                std::ifstream part_colors(partition_filename(partition_id), std::ios::binary);
                m_saver->append(part_colors);
            }

            m_saver->visit(m_partition_endpoints);
        }

        ~builder() {
            for (uint64_t i = 0; i != m_color_sets_builders.size(); ++i) {
                std::remove(partition_filename(i).c_str());
            }
        }

    private:
        bits::compact_vector::builder m_meta_color_sets_builder;
        std::vector<typename ColorSets::builder> m_color_sets_builders;
        std::vector<util::external_saver> m_partition_savers;
        util::external_saver* m_saver = nullptr;

        std::vector<std::unique_ptr<std::shared_mutex>> m_partitions_mutex;
        using hash_id_map =
            std::unordered_map<__uint128_t, std::shared_future<uint32_t>, util::hasher_uint128_t>;
        std::vector<hash_id_map> m_hashes;  // (hash, id)
        std::vector<uint32_t> m_partition_starts;
        std::vector<uint32_t> m_partial_sets_starts;

        uint32_t m_num_colors;
        uint64_t m_offset;
        std::vector<uint64_t> m_meta_color_sets_offsets;
        std::vector<partition_endpoint> m_partition_endpoints;

        std::atomic<bool> m_is_flushing;

        std::string m_tmp_dirname;
        uint64_t m_max_RAM_bytes;
        bool m_verbose;

        std::string partition_filename(const uint64_t partition_id) const {
            return m_tmp_dirname + "/partial_sets_" + std::to_string(partition_id) + ".bin";
        }

        uint32_t encode_partial_color_set(const uint64_t partition_id,
                                          std::span<uint32_t> partial_color_set) {
            assert(partition_id < m_color_sets_builders.size());
            assert(!partial_color_set.empty());
            uint32_t partial_color_set_id;
            auto partition_begin = m_partition_starts[partition_id];
            auto hash = util::hash128(reinterpret_cast<char const*>(partial_color_set.data()),
                                      partial_color_set.size() * sizeof(uint32_t));
            bool requires_compression = false;
            std::promise<uint32_t> promise;
            std::shared_future future = promise.get_future();

            {
                std::lock_guard lock(*m_partitions_mutex[partition_id]);
                const auto it = m_hashes[partition_id].find(hash);

                if (it == m_hashes[partition_id].cend()) {
                    m_hashes[partition_id].insert({hash, future});
                    requires_compression = true;
                } else {
                    future = it->second;
                }
            }
            if (requires_compression) {
                std::ranges::transform(
                    partial_color_set, partial_color_set.begin(),
                    [partition_begin](const uint32_t n) { return n - partition_begin; });
                partial_color_set_id =
                    m_color_sets_builders[partition_id].encode(partial_color_set);
                promise.set_value(partial_color_set_id);
            } else {
                partial_color_set_id = future.get();
            }

            /*  Note: at this stage, partial_color_set_id is relative
             *  to its partition (is not global yet). */
            if (num_bytes() > m_max_RAM_bytes) {
                flush();
            }
            return partial_color_set_id;
        }
    };

    struct iterator_sentinel {};

    struct forward_iterator {
        using sentinel_type = iterator_sentinel;

        forward_iterator(meta<ColorSets> const* ptr, uint64_t begin)
            : m_ptr(ptr)
            , m_begin(begin)
            , m_meta_color_set_size((m_ptr->m_meta_color_sets)[m_begin]) {
            rewind();
        }

        void rewind() {
            init();
            assert(m_meta_color_set_size > 0);
            change_partition();
        }

        void init() {
            m_pos_in_meta_color_list = 0;
            m_partition_id = 0;
            m_partition_min_color = 0;
        }

        uint64_t value() const { return m_curr_val; }
        uint64_t operator*() const { return value(); }

        bool has_next() const { return m_pos_in_curr_partition != m_curr_partition_size; }
        void next_in_partition() {
            m_pos_in_curr_partition += 1;
            m_curr_partition_it.next();
            update_curr_val();
        }

        void next() {
            if (m_pos_in_curr_partition == m_curr_partition_size - 1) {
                if (m_pos_in_meta_color_list == meta_color_set_size() - 1) {  // saturate
                    m_curr_val = num_colors();
                    return;
                }
                m_pos_in_meta_color_list += 1;
                change_partition();
            } else {
                next_in_partition();
            }
        }
        void operator++() { next(); }

        bool operator==(iterator_sentinel) const { return m_curr_val == num_colors(); }

        /* update the state of the iterator to the element
           which is greater-than or equal-to lower_bound */
        void next_geq(const uint64_t lower_bound) {
            assert(lower_bound <= num_colors());
            while (value() < lower_bound) next();
            assert(value() >= lower_bound);
        }

        /* Warning: this might be slow. */
        uint32_t size() const {
            uint64_t n = 0;
            for (uint32_t i = 0, partition_id = 0; i != meta_color_set_size(); ++i) {
                uint32_t meta_color = (m_ptr->m_meta_color_sets)[m_begin + 1 + i];
                partition_id = update_partition_id(meta_color, partition_id);
                uint32_t num_color_sets_before =
                    (m_ptr->m_partition_endpoints)[partition_id].num_color_sets_before;
                n += (m_ptr->m_partial_color_sets)[partition_id]
                         .color_set(meta_color - num_color_sets_before)
                         .size();
            }
            return n;
        }

        uint32_t partial_set_size() const { return m_curr_partition_it.size(); }

        uint32_t meta_color() const { return m_curr_meta_color; }

        void read_partition_id() {
            m_curr_meta_color = (m_ptr->m_meta_color_sets)[m_begin + 1 + m_pos_in_meta_color_list];
            m_partition_id = update_partition_id(m_curr_meta_color, m_partition_id);
        }

        void next_partition_id() {
            m_pos_in_meta_color_list += 1;
            if (m_pos_in_meta_color_list == meta_color_set_size()) {  // saturate
                m_partition_id = num_partitions();
                return;
            }
            read_partition_id();
        }

        void next_geq_partition_id(const uint32_t lower_bound) {
            assert(lower_bound <= num_partitions());
            while (partition_id() < lower_bound) next_partition_id();
            assert(partition_id() >= lower_bound);
        }

        void update_partition() {
            /* update partition min/max color */
            auto const& endpoints = m_ptr->m_partition_endpoints;
            m_partition_min_color = endpoints[m_partition_id].min_color;
            m_partition_max_color = endpoints[m_partition_id + 1].min_color;

            uint32_t num_color_sets_before = endpoints[m_partition_id].num_color_sets_before;
            m_curr_partition_it = (m_ptr->m_partial_color_sets)[m_partition_id].color_set(
                m_curr_meta_color - num_color_sets_before);
            m_curr_partition_size = m_curr_partition_it.size();
            assert(m_curr_partition_size > 0);
            m_pos_in_curr_partition = 0;

            update_curr_val();
        }

        void change_partition() {
            read_partition_id();
            update_partition();
        }

        uint32_t partition_id() const { return m_partition_id; }
        uint32_t meta_color_set_size() const { return m_meta_color_set_size; }
        uint32_t num_colors() const { return m_ptr->num_colors(); }
        uint32_t num_partitions() const { return m_ptr->num_partitions(); }
        uint32_t partition_min_color() const { return m_partition_min_color; }
        uint32_t partition_max_color() const { return m_partition_max_color; }
        uint32_t num_color_sets_before() const {
            return m_ptr->m_partition_endpoints[m_partition_id].num_color_sets_before;
        }

    private:
        meta<ColorSets> const* m_ptr;
        typename ColorSets::iterator_type m_curr_partition_it;
        uint64_t m_begin;
        uint32_t m_curr_meta_color, m_curr_val;
        uint32_t m_meta_color_set_size, m_pos_in_meta_color_list;
        uint32_t m_curr_partition_size, m_pos_in_curr_partition;
        uint32_t m_partition_id;
        uint32_t m_partition_min_color, m_partition_max_color;

        void update_curr_val() { m_curr_val = m_curr_partition_it.value() + m_partition_min_color; }

        uint32_t update_partition_id(const uint32_t meta_color, uint32_t partition_id) const {
            auto const& endpoints = m_ptr->m_partition_endpoints;
            while (partition_id + 1 < endpoints.size() and
                   meta_color >= endpoints[partition_id + 1].num_color_sets_before) {
                partition_id += 1;
            }
            assert(partition_id < m_ptr->num_partitions());
            return partition_id;
        }
    };

    typedef forward_iterator iterator_type;

    forward_iterator color_set(uint64_t color_set_id) const {
        assert(color_set_id < num_color_sets());
        uint64_t begin = m_meta_color_sets_offsets.access(color_set_id);
        return forward_iterator(this, begin);
    }

    essentials::owning_span<ColorSets> const& partial_colors() const {
        return m_partial_color_sets;
    }

    uint32_t num_colors() const { return m_num_colors; }
    uint64_t num_color_sets() const { return m_meta_color_sets_offsets.size() - 1; }
    uint64_t num_partitions() const { return m_partition_endpoints.size() - 1; }

    uint64_t num_bits() const {
        uint64_t num_bits_colors = sizeof(size_t) * 8;  // for std::vector::size
        for (auto const& c : m_partial_color_sets) num_bits_colors += c.num_bits();
        return num_bits_colors +
               (m_meta_color_sets_offsets.num_bytes() + m_meta_color_sets.num_bytes() +
                essentials::vec_bytes(m_partition_endpoints) + sizeof(m_num_colors)) *
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
        visitor.visit(t.m_meta_color_sets);
        visitor.visit(t.m_meta_color_sets_offsets);
        visitor.visit(t.m_partial_color_sets);
        visitor.visit(t.m_partition_endpoints);
    }

    uint32_t m_num_colors;
    bits::compact_vector m_meta_color_sets;
    bits::elias_fano<false, false> m_meta_color_sets_offsets;
    essentials::owning_span<ColorSets> m_partial_color_sets;
    essentials::owning_span<partition_endpoint> m_partition_endpoints;
};

}  // namespace fulgor
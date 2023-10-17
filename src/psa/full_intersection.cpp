#include "index.hpp"
#include "../external/sshash/include/query/streaming_query_canonical_parsing.hpp"

namespace fulgor {

template <typename Iterator>
void intersect(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors) {
    assert(colors.empty());

    if (iterators.empty()) return;

    bool all_very_dense = true;
    for (auto const& it : iterators) {
        if (it.type() != color_classes::hybrid::list_type::complementary_delta_gaps) {
            all_very_dense = false;
            break;
        }
    }

    if (all_very_dense) {
        /* step 1: take the union of complementary sets */
        std::vector<uint32_t> tmp;
        for (auto& it : iterators) it.reinit_for_complemented_set_iteration();

        uint32_t candidate = (*std::min_element(iterators.begin(), iterators.end(),
                                                [](auto const& x, auto const& y) {
                                                    return x.comp_value() < y.comp_value();
                                                }))
                                 .comp_value();

        const uint32_t num_docs = iterators[0].num_docs();
        tmp.reserve(num_docs);
        while (candidate < num_docs) {
            uint32_t next_candidate = num_docs;
            for (uint64_t i = 0; i != iterators.size(); ++i) {
                if (iterators[i].comp_value() == candidate) iterators[i].next_comp();
                /* compute next minimum */
                if (iterators[i].comp_value() < next_candidate) {
                    next_candidate = iterators[i].comp_value();
                }
            }
            tmp.push_back(candidate);
            assert(next_candidate > candidate);
            candidate = next_candidate;
        }

        /* step 2: compute the intersection by scanning tmp */
        candidate = 0;
        for (uint32_t i = 0; i != tmp.size(); ++i) {
            while (candidate < tmp[i]) {
                colors.push_back(candidate);
                candidate += 1;
            }
            candidate += 1;  // skip the candidate because it is equal to tmp[i]
        }
        while (candidate < num_docs) {
            colors.push_back(candidate);
            candidate += 1;
        }

        return;
    }

    /* traditional intersection code based on next_geq() and next() */

    std::sort(iterators.begin(), iterators.end(),
              [](auto const& x, auto const& y) { return x.size() < y.size(); });

    const uint32_t num_docs = iterators[0].num_docs();
    uint32_t candidate = iterators[0].value();
    uint64_t i = 1;
    while (candidate < num_docs) {
        for (; i != iterators.size(); ++i) {
            iterators[i].next_geq(candidate);
            uint32_t val = iterators[i].value();
            if (val != candidate) {
                candidate = val;
                i = 0;
                break;
            }
        }
        if (i == iterators.size()) {
            colors.push_back(candidate);
            iterators[0].next();
            candidate = iterators[0].value();
            i = 1;
        }
    }
}

template <typename Iterator>
void meta_intersect(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors,
                    std::vector<uint32_t>& partition_ids) {
    assert(colors.empty());
    assert(partition_ids.empty());

    if (iterators.empty()) return;

    std::sort(iterators.begin(), iterators.end(), [](auto const& x, auto const& y) {
        return x.meta_color_list_size() < y.meta_color_list_size();
    });

    /* step 1: determine partitions in common */
    const uint32_t num_partitions = iterators[0].num_partitions();
    partition_ids.reserve(num_partitions);  // at most

    uint32_t candidate = iterators[0].partition_id();
    uint64_t i = 1;
    while (candidate < num_partitions) {
        for (; i != iterators.size(); ++i) {
            iterators[i].next_geq_partition_id(candidate);
            uint32_t val = iterators[i].partition_id();
            if (val != candidate) {
                candidate = val;
                i = 0;
                break;
            }
        }
        if (i == iterators.size()) {
            partition_ids.push_back(candidate);
            iterators[0].next_partition_id();
            candidate = iterators[0].partition_id();
            i = 1;
        }
    }

    /* step 2: intersect partial colors in the same partitions only */
    for (auto& it : iterators) {
        it.init();
        it.change_partition();
    }
    for (auto partition_id : partition_ids) {
        bool same_meta_color = true;
        auto& front_it = iterators.front();
        front_it.next_geq_partition_id(partition_id);
        front_it.update_partition();
        uint32_t meta_color = front_it.meta_color();

        for (uint32_t i = 1; i != iterators.size(); ++i) {
            auto& it = iterators[i];
            it.next_geq_partition_id(partition_id);
            it.update_partition();
            if (it.meta_color() != meta_color) same_meta_color = false;
        }

        if (same_meta_color) {  // do not intersect, just write the whole partial color once
            while (front_it.has_next()) {
                uint32_t val = front_it.value();
                colors.push_back(val);
                front_it.next_in_partition();
            }
        } else {  // intersect partial colors in the partition
            const uint32_t num_docs = iterators[0].partition_upper_bound();
            uint32_t candidate = iterators[0].value();
            uint64_t i = 1;
            while (candidate < num_docs) {
                for (; i != iterators.size(); ++i) {
                    iterators[i].next_geq(candidate);
                    uint32_t val = iterators[i].value();
                    if (val != candidate) {
                        candidate = val;
                        i = 0;
                        break;
                    }
                }
                if (i == iterators.size()) {
                    colors.push_back(candidate);
                    iterators[0].next();
                    candidate = iterators[0].value();
                    i = 1;
                }
            }
        }
    }
}

void stream_through(sshash::dictionary const& k2u, std::string const& sequence,
                    std::vector<uint32_t>& unitig_ids) {
    sshash::streaming_query_canonical_parsing query(&k2u);
    query.start();
    const uint64_t num_kmers = sequence.length() - k2u.k() + 1;
    for (uint64_t i = 0, prev_unitig_id = -1; i != num_kmers; ++i) {
        char const* kmer = sequence.data() + i;
        auto answer = query.lookup_advanced(kmer);
        if (answer.kmer_id != sshash::constants::invalid_uint64) {  // kmer is positive
            if (answer.contig_id != prev_unitig_id) {
                unitig_ids.push_back(answer.contig_id);
                prev_unitig_id = answer.contig_id;
            }
        }
    }
}

template <typename ColorClasses>
void index<ColorClasses>::pseudoalign_full_intersection(std::string const& sequence,
                                                        std::vector<uint32_t>& colors) const {
    if (sequence.length() < m_k2u.k()) return;
    colors.clear();
    std::vector<uint32_t> unitig_ids;
    stream_through(m_k2u, sequence, unitig_ids);
    intersect_unitigs(unitig_ids, colors);
}

template <typename ColorClasses>
void index<ColorClasses>::intersect_unitigs(std::vector<uint32_t>& unitig_ids,
                                            std::vector<uint32_t>& colors) const {
    /* here we use it to hold the color class ids;
       in meta_intersect we use it to hold the partition ids */
    std::vector<uint32_t> tmp;
    std::vector<typename ColorClasses::iterator_type> iterators;

    /* deduplicate unitig_ids */
    std::sort(unitig_ids.begin(), unitig_ids.end());
    auto end = std::unique(unitig_ids.begin(), unitig_ids.end());
    tmp.reserve(end - unitig_ids.begin());
    for (auto it = unitig_ids.begin(); it != end; ++it) {
        uint32_t unitig_id = *it;
        uint32_t color_class_id = u2c(unitig_id);
        tmp.push_back(color_class_id);
    }

    /* deduplicate color class ids */
    std::sort(tmp.begin(), tmp.end());
    end = std::unique(tmp.begin(), tmp.end());
    iterators.reserve(end - tmp.begin());
    for (auto it = tmp.begin(); it != end; ++it) {
        uint64_t color_class_id = *it;
        auto fwd_it = m_ccs.colors(color_class_id);
        iterators.push_back(fwd_it);
    }

    tmp.clear();  // don't need color class ids anymore
    if constexpr (ColorClasses::meta_colored) {
        meta_intersect(iterators, colors, tmp);
    } else {
        intersect(iterators, colors);
    }
}

}  // namespace fulgor

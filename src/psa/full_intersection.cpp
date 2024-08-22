#include "include/index.hpp"
#include "external/sshash/include/query/streaming_query_canonical_parsing.hpp"

namespace fulgor {

template <typename Iterator>
void intersect(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors,
               std::vector<uint32_t>& complement_set) {
    assert(colors.empty());
    assert(complement_set.empty());

    if (iterators.empty()) return;

    bool all_very_dense = true;
    for (auto const& it : iterators) {
        if (it.type() != list_type::complement_delta_gaps) {
            all_very_dense = false;
            break;
        }
    }

    if (all_very_dense) {
        /* step 1: take the union of complementary sets */
        for (auto& it : iterators) it.reinit_for_complemented_set_iteration();

        uint32_t candidate = (*std::min_element(iterators.begin(), iterators.end(),
                                                [](auto const& x, auto const& y) {
                                                    return x.comp_value() < y.comp_value();
                                                }))
                                 .comp_value();

        const uint32_t num_colors = iterators[0].num_colors();
        complement_set.reserve(num_colors);
        while (candidate < num_colors) {
            uint32_t next_candidate = num_colors;
            for (uint64_t i = 0; i != iterators.size(); ++i) {
                if (iterators[i].comp_value() == candidate) iterators[i].next_comp();
                /* compute next minimum */
                if (iterators[i].comp_value() < next_candidate) {
                    next_candidate = iterators[i].comp_value();
                }
            }
            complement_set.push_back(candidate);
            assert(next_candidate > candidate);
            candidate = next_candidate;
        }

        /* step 2: compute the intersection by scanning complement_set */
        candidate = 0;
        for (uint32_t i = 0; i != complement_set.size(); ++i) {
            while (candidate < complement_set[i]) {
                colors.push_back(candidate);
                candidate += 1;
            }
            candidate += 1;  // skip the candidate because it is equal to complement_set[i]
        }
        while (candidate < num_colors) {
            colors.push_back(candidate);
            candidate += 1;
        }

        return;
    }

    /* traditional intersection code based on next_geq() and next() */

    std::sort(iterators.begin(), iterators.end(),
              [](auto const& x, auto const& y) { return x.size() < y.size(); });

    const uint32_t num_colors = iterators[0].num_colors();
    uint32_t candidate = iterators[0].value();
    uint64_t i = 1;
    while (candidate < num_colors) {
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
void diff_intersect(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors) {
    assert(colors.empty());

    if (iterators.empty()) return;

    std::sort(iterators.begin(), iterators.end(),
              [](auto const& x, auto const& y) { return x.size() < y.size(); });

    const uint32_t num_colors = iterators[0].num_colors();
    uint32_t candidate = iterators[0].value();
    uint64_t i = 1;
    while (candidate < num_colors) {
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
            const uint32_t num_colors = iterators[0].partition_upper_bound();
            uint32_t candidate = iterators[0].value();
            uint64_t i = 1;
            while (candidate < num_colors) {
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
                    std::vector<uint64_t>& unitig_ids) {
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

template <typename ColorSets>
void index<ColorSets>::pseudoalign_full_intersection(std::string const& sequence,
                                                     std::vector<uint32_t>& colors) const {
    if (sequence.length() < m_k2u.k()) return;
    colors.clear();
    std::vector<uint64_t> unitig_ids;
    stream_through(m_k2u, sequence, unitig_ids);
    intersect_unitigs(unitig_ids, colors);
}

template <typename ColorSets>
void index<ColorSets>::intersect_unitigs(std::vector<uint64_t>& unitig_ids,
                                         std::vector<uint32_t>& colors) const {
    /* here we use it to hold the color class ids;
       in meta_intersect we use it to hold the partition ids */
    std::vector<uint32_t> tmp;
    std::vector<typename ColorSets::iterator_type> iterators;

    /* deduplicate unitig_ids */
    std::sort(unitig_ids.begin(), unitig_ids.end());
    auto end_unitigs = std::unique(unitig_ids.begin(), unitig_ids.end());
    tmp.reserve(end_unitigs - unitig_ids.begin());
    for (auto it = unitig_ids.begin(); it != end_unitigs; ++it) {
        uint32_t unitig_id = *it;
        uint32_t color_set_id = u2c(unitig_id);
        tmp.push_back(color_set_id);
    }

    /* deduplicate color class ids */
    std::sort(tmp.begin(), tmp.end());
    auto end_tmp = std::unique(tmp.begin(), tmp.end());
    iterators.reserve(end_tmp - tmp.begin());
    for (auto it = tmp.begin(); it != end_tmp; ++it) {
        uint64_t color_set_id = *it;
        auto fwd_it = m_color_sets.color_set(color_set_id);
        iterators.push_back(fwd_it);
    }

    tmp.clear();  // don't need color class ids anymore
    if constexpr (ColorSets::meta_colored) {
        meta_intersect(iterators, colors, tmp);
    } else if constexpr (ColorSets::differential_colored) {
        diff_intersect(iterators, colors);
    } else {
        intersect(iterators, colors, tmp);
    }

    assert(util::check_intersection(iterators, colors));
}

}  // namespace fulgor

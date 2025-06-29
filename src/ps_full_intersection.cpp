#include "include/index.hpp"
#include "external/sshash/include/streaming_query.hpp"

namespace fulgor {

template <typename Iterator>
void next_geq_intersect(Iterator const& begin, Iterator const& end,                //
                        std::vector<uint32_t>& colors, const uint32_t num_colors)  //
{
    uint32_t candidate = begin->value();
    uint64_t i = 1;
    const uint64_t size = end - begin;
    while (candidate < num_colors) {
        for (; i != size; ++i) {
            (begin + i)->next_geq(candidate);
            uint32_t val = (begin + i)->value();
            if (val != candidate) {
                candidate = val;
                i = 0;
                break;
            }
        }
        if (i == size) {
            colors.push_back(candidate);
            begin->next();
            candidate = begin->value();
            i = 1;
        }
    }
}

template <typename Iterator>
void intersect(std::vector<Iterator>& iterators,       //
               std::vector<uint32_t>& colors,          //
               std::vector<uint32_t>& complement_set)  //
{
    assert(colors.empty());
    assert(complement_set.empty());

    if (iterators.empty()) return;
    std::sort(iterators.begin(), iterators.end(),
              [](auto const& x, auto const& y) { return x.size() < y.size(); });

    const uint32_t num_colors = iterators[0].num_colors();
    uint32_t num_sparse = 0;
    while (num_sparse != iterators.size() &&
           iterators[num_sparse].encoding_type() != encoding_t::complement_delta_gaps) {
        ++num_sparse;
    }

    if (num_sparse == 0) {
        /* step 1: take the union of complementary sets */
        for (auto& it : iterators) it.reinit_for_complemented_set_iteration();

        uint32_t candidate = (*std::min_element(iterators.begin(), iterators.end(),
                                                [](auto const& x, auto const& y) {
                                                    return x.comp_value() < y.comp_value();
                                                }))
                                 .comp_value();

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
        for (uint32_t i : complement_set) {
            while (candidate < i) {
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

    std::vector<bool> complement_union(num_colors, true);
    for (uint32_t i = num_sparse; i < iterators.size(); ++i) {
        auto it = iterators[i];
        it.reinit_for_complemented_set_iteration();
        while (it.comp_value() < num_colors) {
            complement_union[it.comp_value()] = false;
            it.next_comp();
        }
    }

    /* traditional intersection code based on next_geq() and next() */

    assert(iterators[0].encoding_type() != encoding_t::complement_delta_gaps);
    uint32_t candidate = iterators[0].value();
    uint64_t i = 1;
    uint64_t size = num_sparse;

    while (candidate < num_colors) {
        for (; i != size; ++i) {
            iterators[i].next_geq(candidate);
            uint32_t val = iterators[i].value();
            if (val != candidate) {
                candidate = val;
                i = 0;
                break;
            }
        }
        if (i == size) {
            if (complement_union[candidate]) { colors.push_back(candidate); }
            iterators[0].next();
            candidate = iterators[0].value();
            i = 1;
        }
    }
}

template <typename Iterator>
void diff_intersect(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors,
                    uint32_t lower_bound = 0)  //
{
    if (iterators.empty()) return;
    const uint32_t num_colors = iterators[0].num_colors();

    std::sort(iterators.begin(), iterators.end(), [](const Iterator& a, const Iterator& b) {
        return a.representative_begin() < b.representative_begin();
    });

    const uint32_t num_iterators = iterators.size();
    uint32_t num_partitions = 1;
    {
        uint32_t prev_partition = iterators[0].representative_begin();
        for (const auto& it : iterators) {
            uint32_t partition_id = it.representative_begin();
            if (partition_id != prev_partition) {
                prev_partition = partition_id;
                ++num_partitions;
            }
        }
    }

    std::vector<std::vector<uint32_t>> partitions(num_partitions);

    {
        std::vector<uint32_t> counts(num_colors, 0);
        uint32_t partition_id = 0;
        uint32_t partition_size = 0;

        for (uint32_t iterator_id = 0; iterator_id < num_iterators; iterator_id++) {
            Iterator it = iterators[iterator_id];
            partition_size++;

            bool is_last_in_partition =
                iterator_id + 1 == num_iterators ||
                iterators[iterator_id + 1].representative_begin() != it.representative_begin();

            if (partition_size == 1 && is_last_in_partition) {
                // if one element in partition, decode the color set
                for (uint32_t i = 0; i < it.size(); ++i, ++it) {
                    partitions[partition_id].push_back(*it);
                }
                partition_id++;
                partition_size = 0;
                continue;
            }

            it.full_rewind();

            uint32_t val = it.differential_val();
            while (val != num_colors) {
                ++counts[val];
                it.next_differential_val();
                val = it.differential_val();
            }

            if (is_last_in_partition) {
                it.full_rewind();
                val = it.representative_val();
                for (uint32_t color = 0; color < num_colors; color++) {
                    if (val < color) {
                        it.next_representative_val();
                        val = it.representative_val();
                    }
                    if ((counts[color] == partition_size && val != color) ||
                        (counts[color] == 0 && val == color)) {
                        partitions[partition_id].push_back(color);
                    }
                }
                partition_id++;
                partition_size = 0;
                fill(counts.begin(), counts.end(), 0);
            }
        }
    }

    std::sort(partitions.begin(), partitions.end(),
              [](auto const& x, auto const& y) { return x.size() < y.size(); });

    std::vector<std::vector<uint32_t>::iterator> its(num_partitions);
    for (uint32_t i = 0; i < num_partitions; i++) {
        if (partitions[i].empty()) return;
        its[i] = partitions[i].begin();
    }

    uint32_t candidate = *its[0];
    uint64_t i = 1;
    while (candidate < num_colors) {
        for (; i != its.size(); ++i) {
            while (its[i] != partitions[i].end() && *its[i] < candidate) ++its[i];
            if (its[i] == partitions[i].end()) {
                candidate = num_colors;
                break;
            }
            uint32_t val = *its[i];
            if (val != candidate) {
                candidate = val;
                i = 0;
                break;
            }
        }
        if (i == its.size()) {
            colors.push_back(candidate + lower_bound);
            ++its[0];
            if (its[0] == partitions[0].end()) break;
            candidate = *its[0];
            i = 1;
        }
    }
}

template <typename Iterator, bool is_differential>
void meta_intersect(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors,
                    std::vector<uint32_t>& partition_ids) {
    assert(colors.empty());
    assert(partition_ids.empty());

    if (iterators.empty()) return;

    for (auto it : iterators) {
        while (it.partition_id() != it.num_partitions()) it.next_partition_id();
        it.init();
        it.change_partition();
    }

    std::sort(iterators.begin(), iterators.end(), [](auto const& x, auto const& y) {
        return x.meta_color_set_size() < y.meta_color_set_size();
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
            std::sort(iterators.begin(), iterators.end(), [](const Iterator& a, const Iterator& b) {
                return a.partial_set_size() < b.partial_set_size() ||
                       (a.partial_set_size() == b.partial_set_size() &&
                        a.meta_color() < b.meta_color());
            });

            uint64_t back_pos = 0;
            for (uint64_t curr_pos = 1; curr_pos < iterators.size(); curr_pos++) {
                if (iterators[curr_pos].meta_color() != iterators[back_pos].meta_color()) {
                    std::swap(iterators[++back_pos], iterators[curr_pos]);
                }
            }
            auto end_it = iterators.begin() + back_pos + 1;

            if constexpr (is_differential) {
                std::vector<differential::iterator_type> diff_iterators;
                diff_iterators.reserve(end_it - iterators.begin());
                std::transform(iterators.begin(), end_it, back_inserter(diff_iterators),
                               [](Iterator a) { return a.partition_it(); });
                uint32_t lower_bound =
                    iterators[0].partition_max_color() - diff_iterators[0].num_colors();
                diff_intersect(diff_iterators, colors, lower_bound);
            } else {
                const uint32_t num_colors = iterators[0].partition_max_color();
                next_geq_intersect(iterators.begin(), end_it, colors, num_colors);
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

    { /* stream through */
        sshash::streaming_query<kmer_type, true> query(&m_k2u);
        query.reset();
        const uint64_t num_kmers = sequence.length() - m_k2u.k() + 1;
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
    if constexpr (ColorSets::type == index_t::META) {
        meta_intersect<typename ColorSets::iterator_type, false>(iterators, colors, tmp);
    } else if constexpr (ColorSets::type == index_t::META_DIFF) {
        meta_intersect<typename ColorSets::iterator_type, true>(iterators, colors, tmp);
    } else if constexpr (ColorSets::type == index_t::DIFF) {
        diff_intersect(iterators, colors);
    } else if constexpr (ColorSets::type == index_t::HYBRID) {
        intersect(iterators, colors, tmp);
    }

    assert(util::check_intersection(iterators, colors));
}

}  // namespace fulgor

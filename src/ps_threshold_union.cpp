#include <numeric>  // for std::accumulate

#include "include/index.hpp"
#include "external/sshash/include/query/streaming_query_canonical_parsing.hpp"

namespace fulgor {

template <typename T>
struct scored {
    T item;
    uint32_t score;
};

typedef scored<uint64_t> scored_id;

template <typename Iterator>
void merge(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors, int64_t min_score) {
    if (iterators.empty()) return;

    uint32_t num_colors = iterators[0].item.num_colors();
    std::vector<int32_t> scores(num_colors, 0);
    for (auto& it : iterators) {
        if (it.item.type() == list_type::complement_delta_gaps) {
            it.item.reinit_for_complemented_set_iteration();
            min_score -= it.score;
            while (it.item.comp_value() < num_colors) {
                scores[it.item.comp_value()] -= it.score;
                it.item.next_comp();
            }
        } else {
            uint32_t size = it.item.size();
            for (uint32_t i = 0; i < size; ++i, it.item.next()) {
                scores[it.item.value()] += it.score;
            }
        }
    }
    for (uint32_t color = 0; color < num_colors; color++) {
        if (scores[color] >= min_score) colors.push_back(color);
    }
}

template <typename Iterator>
void merge_meta(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors,
                const uint64_t min_score) {
    if (iterators.empty()) return;

    const uint32_t num_partitions = iterators[0].item.num_partitions();
    const uint32_t num_colors = iterators[0].item.num_colors();
    std::vector<uint32_t> partition_ids;
    partition_ids.reserve(num_partitions);

    // the number of partitions is relatively small, so this does not impact efficiency
    uint32_t candidate_partition =
        (*std::min_element(iterators.begin(), iterators.end(), [](auto const& x, auto const& y) {
            return x.item.partition_id() < y.item.partition_id();
        })).item.partition_id();

    while (candidate_partition < num_partitions) {
        uint32_t next_partition = num_partitions;
        uint32_t score = 0;
        for (uint64_t i = 0; i != iterators.size(); ++i) {
            if (iterators[i].item.partition_id() == candidate_partition) {
                score += iterators[i].score;
                iterators[i].item.next_partition_id();
            }
            /* compute next minimum */
            if (iterators[i].item.partition_id() < next_partition) {
                next_partition = iterators[i].item.partition_id();
            }
        }
        if (score >= min_score) partition_ids.push_back(candidate_partition);
        assert(next_partition > candidate_partition);
        candidate_partition = next_partition;
    }

    std::vector<uint32_t> scores(num_colors, 0);
    for (auto& it : iterators) {
        it.item.init();
        it.item.change_partition();
    }
    for (auto partition_id : partition_ids) {
        uint32_t upper_bound = 0;
        for (auto& it : iterators) {
            it.item.next_geq_partition_id(partition_id);
            if (it.item.partition_id() == partition_id) {
                it.item.update_partition();
                upper_bound = it.item.partition_upper_bound();
            }
        }

        std::sort(iterators.begin(), iterators.end(), [](const Iterator& a, const Iterator& b) {
            return a.item.partition_id() < b.item.partition_id() ||
                   (a.item.partition_id() == b.item.partition_id() &&
                    a.item.meta_color() < b.item.meta_color());
        });

        uint64_t meta_score = iterators.front().score;
        auto process_meta = [&](Iterator& it) {
            while (it.item.value() < upper_bound) {
                scores[it.item.value()] += meta_score;
                it.item.next();
            }
        };
        uint64_t i = 1;
        for (; i < iterators.size(); ++i) {
            auto& it = iterators[i];
            if (it.item.partition_id() != partition_id) break;
            if (it.item.meta_color() != iterators[i - 1].item.meta_color()) {
                process_meta(iterators[i - 1]);
                meta_score = 0;
            }
            meta_score += it.score;
        }
        process_meta(iterators[i - 1]);
    }

    for (uint32_t color = 0; color < num_colors; color++) {
        if (scores[color] >= min_score) colors.push_back(color);
    }
}

template <typename Iterator>
void merge_diff(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors,
                const uint64_t min_score) {
    if (iterators.empty()) return;
    const uint32_t num_colors = iterators[0].item.num_colors();
    const uint32_t num_iterators = iterators.size();

    std::sort(iterators.begin(), iterators.end(), [](const Iterator& a, const Iterator& b) {
        return a.item.representative_begin() < b.item.representative_begin();
    });

    std::vector<uint32_t> partition_scores(num_colors, 0);
    std::vector<uint32_t> scores(num_colors, 0);
    uint32_t score = 0;
    uint32_t partition_size = 0;
    for (uint32_t iterator_id = 0; iterator_id < num_iterators; iterator_id++) {
        Iterator it = iterators[iterator_id];
        partition_size++;
        score += it.score;

        bool is_last_in_partition = iterator_id + 1 == num_iterators ||
                                    iterators[iterator_id + 1].item.representative_begin() !=
                                        it.item.representative_begin();

        if (partition_size == 1 && is_last_in_partition) {
            for (uint32_t i = 0; i < it.item.size(); ++i, it.item.next()) {
                scores[it.item.value()] += it.score;
            }
            score = 0;
            partition_size = 0;
            continue;
        }

        it.item.full_rewind();

        uint32_t val = it.item.differential_val();
        while (val != num_colors) {
            partition_scores[val] += it.score;
            it.item.next_differential_val();
            val = it.item.differential_val();
        }

        if (is_last_in_partition) {
            it.item.full_rewind();
            val = it.item.representative_val();
            for (uint32_t color = 0; color < num_colors; color++) {
                if (val == color) {
                    scores[color] += score - partition_scores[color];
                    it.item.next_representative_val();
                    val = it.item.representative_val();
                } else {
                    scores[color] += partition_scores[color];
                }
            }
            score = 0;
            partition_size = 0;
            fill(partition_scores.begin(), partition_scores.end(), 0);
        }
    }

    for (uint32_t color = 0; color < num_colors; color++) {
        if (scores[color] >= min_score) colors.push_back(color);
    }
}

template <typename Iterator>
void merge_metadiff(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors,
                    const uint64_t min_score) {
    if (iterators.empty()) return;

    const uint32_t num_partitions = iterators[0].item.num_partitions();
    const uint32_t num_colors = iterators[0].item.num_colors();
    const uint32_t num_iterators = iterators.size();
    std::vector<uint32_t> partition_ids;
    partition_ids.reserve(num_partitions);

    // the number of partitions is relatively small, so this does not impact efficiency
    uint32_t candidate_partition =
        std::min_element(iterators.begin(), iterators.end(), [](auto const& x, auto const& y) {
            return x.item.partition_id() < y.item.partition_id();
        })->item.partition_id();

    while (candidate_partition < num_partitions) {
        uint32_t next_partition = num_partitions;
        uint32_t score = 0;
        for (uint64_t i = 0; i != iterators.size(); ++i) {
            if (iterators[i].item.partition_id() == candidate_partition) {
                score += iterators[i].score;
                iterators[i].item.next_partition_id();
            }
            /* compute next minimum */
            if (iterators[i].item.partition_id() < next_partition) {
                next_partition = iterators[i].item.partition_id();
            }
        }
        if (score >= min_score) partition_ids.push_back(candidate_partition);
        assert(next_partition > candidate_partition);
        candidate_partition = next_partition;
    }

    std::vector<uint32_t> scores(num_colors, 0);
    std::vector<uint32_t> partition_scores(num_colors, 0);
    for (auto& it : iterators) {
        it.item.init();
        it.item.change_partition();
    }
    for (auto partition_id : partition_ids) {
        uint32_t num_partition_colors = 0;
        uint32_t lower_bound = 0;
        uint32_t num_sets = 0;
        for (auto& it : iterators) {
            it.item.next_geq_partition_id(partition_id);
            if (it.item.partition_id() == partition_id) {
                it.item.update_partition();
                num_sets++;
            }
        }

        std::sort(iterators.begin(), iterators.end(), [&](const Iterator& a, const Iterator& b) {
            uint32_t a_part = a.item.partition_id();
            uint32_t b_part = b.item.partition_id();
            if (a_part == partition_id && b_part == partition_id) {
                uint32_t a_repr = a.item.partition_it().representative_begin();
                uint32_t a_meta = a.item.meta_color();
                uint32_t b_repr = b.item.partition_it().representative_begin();
                uint32_t b_meta = b.item.meta_color();
                return a_meta < b_meta || (a_meta == b_meta && a_repr < b_repr);
            }
            return a_part < b_part;
        });

        lower_bound = iterators.front().item.partition_lower_bound();
        num_partition_colors = iterators.front().item.partition_it().num_colors();

        uint32_t partition_score = 0;
        uint32_t partition_size = 0;
        uint32_t meta_score = 0;
        for (uint32_t iterator_id = 0; iterator_id < num_iterators; iterator_id++) {
            Iterator it = iterators[iterator_id];
            if (it.item.partition_id() != partition_id) break;
            meta_score += it.score;
            num_sets--;
            partition_size++;
            if (num_sets != 0 &&
                iterators[iterator_id + 1].item.meta_color() == it.item.meta_color())
                continue;

            auto diff_it = it.item.partition_it();
            partition_score += meta_score;

            bool is_last_in_partition =
                num_sets == 0 ||
                iterators[iterator_id + 1].item.partition_it().representative_begin() !=
                    diff_it.representative_begin();

            if (is_last_in_partition && partition_size == 1) {
                for (uint32_t i = 0; i < diff_it.size(); ++i, ++diff_it) {
                    scores[lower_bound + *diff_it] += meta_score;
                }
                partition_score = 0;
                partition_size = 0;
                meta_score = 0;
                continue;
            }

            diff_it.full_rewind();

            uint32_t val = diff_it.differential_val();
            while (val != num_partition_colors) {
                partition_scores[val] += meta_score;
                diff_it.next_differential_val();
                val = diff_it.differential_val();
            }
            meta_score = 0;

            if (is_last_in_partition) {
                diff_it.full_rewind();
                val = diff_it.representative_val();
                for (uint32_t color = 0; color < num_partition_colors; color++) {
                    if (val == color) {
                        scores[lower_bound + color] += partition_score - partition_scores[color];
                        diff_it.next_representative_val();
                        val = diff_it.representative_val();
                    } else {
                        scores[lower_bound + color] += partition_scores[color];
                    }
                }
                partition_score = 0;
                partition_size = 0;
                fill(partition_scores.begin(), partition_scores.begin() + num_partition_colors, 0);
            }
        }
    }
    for (uint32_t color = 0; color < num_colors; color++) {
        if (scores[color] >= min_score) colors.push_back(color);
    }
}

template <typename ColorSets>
void index<ColorSets>::pseudoalign_threshold_union(std::string const& sequence,
                                                   std::vector<uint32_t>& colors,
                                                   const double threshold) const {
    if (sequence.length() < m_k2u.k()) return;
    colors.clear();

    std::vector<scored_id> unitig_ids;
    uint64_t num_positive_kmers_in_sequence = 0;
    { /* stream through with multiplicities */
        sshash::streaming_query_canonical_parsing<kmer_type> query(&m_k2u);
        query.start();
        const uint64_t num_kmers = sequence.length() - m_k2u.k() + 1;
        for (uint64_t i = 0, prev_unitig_id = -1; i != num_kmers; ++i) {
            char const* kmer = sequence.data() + i;
            auto answer = query.lookup_advanced(kmer);
            if (answer.kmer_id != sshash::constants::invalid_uint64) {  // kmer is positive
                num_positive_kmers_in_sequence += 1;
                if (answer.contig_id != prev_unitig_id) {
                    unitig_ids.push_back({answer.contig_id, 1});
                    prev_unitig_id = answer.contig_id;
                } else {
                    assert(!unitig_ids.empty());
                    unitig_ids.back().score += 1;
                }
            }
        }
    }

    /* num_positive_kmers_in_sequence must be equal to the sum of the scores  */
    assert(num_positive_kmers_in_sequence ==
           std::accumulate(unitig_ids.begin(), unitig_ids.end(), uint64_t(0),
                           [](uint64_t curr_sum, auto const& u) { return curr_sum + u.score; }));

    std::vector<scored_id> color_set_ids;
    std::vector<scored<typename ColorSets::iterator_type>> iterators;

    /* deduplicate unitig_ids */
    std::sort(unitig_ids.begin(), unitig_ids.end(),
              [](auto const& x, auto const& y) { return x.item < y.item; });
    uint32_t prev_unitig_id = -1;
    for (uint64_t i = 0; i != unitig_ids.size(); ++i) {
        uint32_t unitig_id = unitig_ids[i].item;
        if (unitig_id != prev_unitig_id) {
            uint32_t color_set_id = u2c(unitig_id);
            color_set_ids.push_back({color_set_id, unitig_ids[i].score});
            prev_unitig_id = unitig_id;
        } else {
            assert(!color_set_ids.empty());
            color_set_ids.back().score += unitig_ids[i].score;
        }
    }

    /* deduplicate color_set_ids */
    std::sort(color_set_ids.begin(), color_set_ids.end(),
              [](auto const& x, auto const& y) { return x.item < y.item; });
    uint32_t prev_color_set_id = -1;
    for (uint64_t i = 0; i != color_set_ids.size(); ++i) {
        uint64_t color_set_id = color_set_ids[i].item;
        if (color_set_id != prev_color_set_id) {
            auto fwd_it = m_color_sets.color_set(color_set_id);
            iterators.push_back({fwd_it, color_set_ids[i].score});
            prev_color_set_id = color_set_id;
        } else {
            assert(!iterators.empty());
            iterators.back().score += color_set_ids[i].score;
        }
    }

    const uint64_t min_score = static_cast<double>(num_positive_kmers_in_sequence) * threshold;

    if constexpr (ColorSets::type == index_t::META) {
        merge_meta(iterators, colors, min_score);
    } else if constexpr (ColorSets::type == index_t::DIFF) {
        merge_diff(iterators, colors, min_score);
    } else if constexpr (ColorSets::type == index_t::META_DIFF) {
        merge_metadiff(iterators, colors, min_score);
    } else if constexpr (ColorSets::type == index_t::HYBRID) {
        merge(iterators, colors, min_score);
    }

    assert(util::check_union(iterators, colors, min_score));
}

}  // namespace fulgor

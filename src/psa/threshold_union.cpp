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
void merge(std::vector<Iterator>& iterators, std::vector<uint32_t>& colors,
           const uint64_t min_score) {
    if (iterators.empty()) return;

    uint32_t candidate =
        (*std::min_element(iterators.begin(), iterators.end(), [](auto const& x, auto const& y) {
            return x.item.value() < y.item.value();
        })).item.value();

    const uint32_t num_colors = iterators[0].item.num_colors();
    while (candidate < num_colors) {
        uint32_t next_candidate = num_colors;
        uint32_t score = 0;
        for (uint64_t i = 0; i != iterators.size(); ++i) {
            if (iterators[i].item.value() == candidate) {
                score += iterators[i].score;
                iterators[i].item.next();
            }
            /* compute next minimum */
            if (iterators[i].item.value() < next_candidate) {
                next_candidate = iterators[i].item.value();
            }
        }
        if (score >= min_score) colors.push_back(candidate);
        assert(next_candidate > candidate);
        candidate = next_candidate;
    }
}

uint64_t stream_through_with_multiplicities(sshash::dictionary const& k2u,
                                            std::string const& sequence,
                                            std::vector<scored_id>& unitig_ids) {
    sshash::streaming_query_canonical_parsing query(&k2u);
    query.start();
    const uint64_t num_kmers = sequence.length() - k2u.k() + 1;
    uint64_t num_positive_kmers_in_sequence = 0;
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
    return num_positive_kmers_in_sequence;
}

template <typename ColorSets>
void index<ColorSets>::pseudoalign_threshold_union(std::string const& sequence,
                                                   std::vector<uint32_t>& colors,
                                                   const double threshold) const {
    if (sequence.length() < m_k2u.k()) return;
    colors.clear();

    std::vector<scored_id> unitig_ids;
    uint64_t num_positive_kmers_in_sequence =
        stream_through_with_multiplicities(m_k2u, sequence, unitig_ids);

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

    /* as Themisto does */
    uint64_t min_score = static_cast<double>(num_positive_kmers_in_sequence) * threshold;

    /* as Bifrost and Metagraph do */
    // uint64_t num_kmers_in_sequence = sequence.length() - m_k2u.k() + 1;
    // uint64_t min_score = static_cast<double>(num_kmers_in_sequence) * threshold;

    merge(iterators, colors, min_score);
}

}  // namespace fulgor

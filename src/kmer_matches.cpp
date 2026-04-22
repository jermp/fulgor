#include "include/index.hpp"
#include "external/sshash/include/streaming_query.hpp"

namespace fulgor {

template <typename ColorSets>
void index<ColorSets>::kmer_matches(std::string const& sequence,
                                    bits::bit_vector::builder& positive_kmers_in_sequence,
                                    std::vector<count_type>& counts) const  //
{
    if (sequence.length() < m_k2u.k()) return;

    const uint64_t num_kmers = sequence.length() - m_k2u.k() + 1;
    positive_kmers_in_sequence.resize(num_kmers, 0);
    std::fill(counts.begin(), counts.end(), 0);
    sshash::streaming_query<sshash_type, true> query(&m_k2u);
    query.reset();

    for (uint64_t i = 0; i != num_kmers; ++i) {
        char const* kmer = sequence.data() + i;
        auto answer = query.lookup(kmer);
        if (answer.kmer_id != sshash::constants::invalid_uint64) {  // kmer is positive
            positive_kmers_in_sequence.set(i);
            uint64_t color_set_id = u2c(answer.string_id);
            auto it = color_set(color_set_id);
            uint64_t color_set_size = it.size();
            for (uint64_t i = 0; i != color_set_size; ++i, ++it) counts[*it] += 1;
        }
    }
}

}  // namespace fulgor

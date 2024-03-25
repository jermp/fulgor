#pragma once

#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "projected_hits.hpp"
#include "external/sshash/include/query/streaming_query_canonical_parsing.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>

namespace piscem_psa {

template <typename FulgorIndex>
class hit_searcher {
    enum class ExpansionTerminationType : uint8_t { MISMATCH = 0, CONTIG_END, READ_END };

public:
    explicit hit_searcher(FulgorIndex const* pfi) : pfi_(pfi) {
        k = static_cast<size_t>(pfi_->k());
    }

    bool get_raw_hits_sketch(std::string const& read, sshash::streaming_query_canonical_parsing& qc,
                             bool isLeft = false, bool verbose = false);

    void clear() {
        left_rawHits.clear();
        right_rawHits.clear();
    }

    void setAltSkip(uint32_t as) { altSkip = as; }
    inline std::vector<std::pair<int, projected_hits>>& get_left_hits() { return left_rawHits; }
    inline std::vector<std::pair<int, projected_hits>>& get_right_hits() { return right_rawHits; }

private:
    FulgorIndex const* pfi_;
    size_t k;
    uint32_t altSkip{5};
    std::vector<std::pair<int, projected_hits>> left_rawHits;
    std::vector<std::pair<int, projected_hits>> right_rawHits;
};

}  // namespace piscem_psa

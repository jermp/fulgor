#pragma once

#include <vector>
#include <utility>
#include <string_view>

#include "piscem_psa/CanonicalKmer.hpp"
#include "piscem_psa/CanonicalKmerIterator.hpp"
#include "piscem_psa/projected_hits.hpp"

template <typename FulgorIndex>
void match(const char* s, int l, FulgorIndex const* idx,
           std::vector<std::pair<projected_hits, int>>& v);

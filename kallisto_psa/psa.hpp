#pragma once

#include <vector>
#include <utility>
#include <string_view>

#include "../piscem_psa/CanonicalKmer.hpp"
#include "../piscem_psa/CanonicalKmerIterator.hpp"
#include "../piscem_psa/projected_hits.hpp"
#include "../include/index_types.hpp"

void match(const char* s, int l, fulgor::index_type const* idx,
           std::vector<std::pair<projected_hits, int>>& v);

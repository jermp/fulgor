#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>
#include <chrono>
#include <algorithm>  // for std::set_intersection

#include "external/smhasher/src/City.h"
#include "external/smhasher/src/City.cpp"

namespace fulgor {

enum index_t { HYBRID, DIFF, META, META_DIFF };
enum encoding_t { delta_gaps, bitmap, complement_delta_gaps, symmetric_difference };

namespace constants {

constexpr double invalid_threshold = -1.0;
constexpr uint64_t default_ram_limit_in_GiB = 8;
static const std::string default_tmp_dirname(".");
static const std::string fulgor_filename_extension("fur");
static const std::string meta_colored_fulgor_filename_extension("mfur");
static const std::string diff_colored_fulgor_filename_extension("dfur");
static const std::string meta_diff_colored_fulgor_filename_extension("mdfur");

namespace current_version_number {
constexpr uint8_t x = 4;
constexpr uint8_t y = 0;
constexpr uint8_t z = 0;
}  // namespace current_version_number

}  // namespace constants

struct build_configuration {
    build_configuration()
        : k(31)
        , m(20)
        , num_threads(1)
        , ram_limit_in_GiB(constants::default_ram_limit_in_GiB)
        , num_colors(0)
        , tmp_dirname(constants::default_tmp_dirname)
        //
        , verbose(false)
        , canonical_parsing(true)
        , check(false)
        //
        , meta_colored(false)
        , diff_colored(false)  //
    {}

    uint32_t k;            // kmer length
    uint32_t m;            // minimizer length
    uint32_t num_threads;  // for building and checking correctness
    uint32_t ram_limit_in_GiB;
    uint64_t num_colors;

    std::string tmp_dirname;
    std::string file_base_name;
    std::string filenames_list;

    std::string index_filename_to_partition;

    bool verbose;
    bool canonical_parsing;
    bool check;

    bool meta_colored;
    bool diff_colored;
};

namespace util {

void print_cmd(int argc, char** argv) {
    for (int i = 0; i != argc; ++i) std::cout << argv[i] << ' ';
    std::cout << std::endl;
}

std::string filename(std::string const& path) { return path.substr(path.find_last_of("/\\") + 1); }

void check_version_number(essentials::version_number const& vnum) {
    if (vnum.x != constants::current_version_number::x) {
        throw std::runtime_error("MAJOR index version mismatch: Fulgor index needs rebuilding");
    }
}

template <typename ForwardIterator>
bool check_intersection(std::vector<ForwardIterator>& iterators,
                        std::vector<uint32_t> const& got)  //
{
    if (iterators.empty()) return true;

    /* re-init iterators */
    for (auto& it : iterators) it.rewind();

    /* decompress the color sets */
    const uint32_t num_colors = iterators[0].num_colors();
    std::vector<std::vector<uint32_t>> sets(iterators.size());
    for (uint64_t i = 0; i != iterators.size(); ++i) {
        auto& it = iterators[i];
        uint32_t val = it.value();
        while (val < num_colors) {
            sets[i].push_back(val);
            it.next();
            val = it.value();
        }
    }

    /* compute intersectiom using std::set_intersection */
    std::vector<uint32_t> expected;
    if (iterators.size() > 1) {
        std::vector<uint32_t> l = sets[0];
        for (uint64_t i = 1; i != sets.size(); ++i) {
            auto r = sets[i];
            expected.clear();
            std::set_intersection(l.begin(), l.end(), r.begin(), r.end(),
                                  std::back_inserter(expected));
            if (i != sets.size() - 1) l.swap(expected);
        }
    } else {
        expected.swap(sets[0]);
    }

    /* compare the results */
    if (expected.size() != got.size()) {
        std::cerr << "expected intersection size " << expected.size() << " but got " << got.size()
                  << std::endl;
        return false;
    }
    for (uint64_t i = 0; i != got.size(); ++i) {
        if (expected[i] != got[i]) {
            std::cerr << "error at " << i << "/" << got.size() << ": expected " << expected[i]
                      << " but got " << got[i] << std::endl;
            return false;
        }
    }

    return true;
}

template <typename ForwardIterator>
bool check_union(std::vector<ForwardIterator>& iterators,                     //
                 std::vector<uint32_t> const& got, const uint64_t min_score)  //
{
    if (iterators.empty()) return true;

    /* re-init iterators */
    for (auto& p : iterators) p.item.rewind();

    /* compute the num. occs of each color */
    const uint32_t num_colors = iterators[0].item.num_colors();
    std::vector<uint32_t> scores(num_colors, 0);
    for (auto& [it, score] : iterators) {
        uint32_t val = it.value();
        while (val < num_colors) {
            scores[val] += score;
            it.next();
            val = it.value();
        }
    }

    /* compare the results */
    uint64_t expected_size = 0;
    auto it = got.begin();
    for (uint64_t i = 0; i != num_colors; ++i) {
        if (scores[i] >= min_score) {
            if (it == got.end()) {
                std::cerr << "error: more elements than expected in thershold-union result"
                          << std::endl;
                return false;
            }
            if (i != *it) {
                std::cerr << "error at " << expected_size << "/" << got.size() << ": expected " << i
                          << " but got " << *it << std::endl;
                return false;
            }
            ++expected_size;
            ++it;
        }
    }

    if (expected_size != got.size()) {
        std::cerr << "expected thershold-union size " << expected_size << " but got " << got.size()
                  << std::endl;
        return false;
    }

    return true;
}

__uint128_t hash128(char const* bytes, uint64_t num_bytes, const uint64_t seed = 1234567890) {
    auto ret = CityHash128WithSeed(bytes, num_bytes, {seed, seed});
    __uint128_t out = 0;
    out += __uint128_t(ret.first);
    out += __uint128_t(ret.second) << 64;
    return out;
}

struct hasher_uint128_t {
    uint64_t operator()(const __uint128_t x) const { return static_cast<uint64_t>(x) ^ (x >> 64); }
};

}  // namespace util
}  // namespace fulgor

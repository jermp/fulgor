#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>
#include <chrono>

#include "GGCAT.hpp"

namespace fulgor {

namespace constants {
constexpr double invalid_threshold = -1.0;
constexpr uint64_t default_ram_limit_in_GiB = 8;
static const std::string default_tmp_dirname(".");
}  // namespace constants

struct build_configuration {
    build_configuration()
        : k(31)
        , m(20)
        , num_threads(1)
        , ram_limit_in_GiB(constants::default_ram_limit_in_GiB)
        , num_docs(0)
        , tmp_dirname(constants::default_tmp_dirname)
        , verbose(false)
        , canonical_parsing(true)
        , check(false) {}

    uint32_t k;            // kmer length
    uint32_t m;            // minimizer length
    uint64_t num_threads;  // for building and checking correctness
    uint64_t ram_limit_in_GiB;
    uint64_t num_docs;
    std::string tmp_dirname;
    std::string file_base_name;
    bool verbose;
    bool canonical_parsing;
    bool check;
};

struct ccdbg_builder {
    ccdbg_builder() : ggcat(new GGCAT()) {}
    ~ccdbg_builder() { delete ggcat; }

    void build_ccdbg(std::string const& filenames_list) {
        ggcat->build(filenames_list, config.ram_limit_in_GiB, config.k, config.num_threads,
                     config.tmp_dirname, config.file_base_name);
        config.num_docs = ggcat->num_docs();
    }

    build_configuration config;
    GGCAT* ggcat;
};

namespace util {

static void print_cmd(int argc, char** argv) {
    for (int i = 0; i != argc; ++i) std::cout << argv[i] << ' ';
    std::cout << std::endl;
}

/* return the number of 64-bit words for num_bits */
static uint64_t num_64bit_words_for(uint64_t num_bits) { return (num_bits + 64 - 1) / 64; }

/*
    Good reference for built-in functions:
    http://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html
*/

/* position of the most significant bit (msb) */
static uint32_t msbll(uint64_t x) {
    assert(x > 0);                   // if x is 0, the result is undefined
    return 63 - __builtin_clzll(x);  // count leading zeros (clz)
}

/* position of the least significant bit (lsb) */
static uint64_t lsbll(uint64_t x) {
    assert(x > 0);              // if x is 0, the result is undefined
    return __builtin_ctzll(x);  // count trailing zeros (ctz)
}

inline uint8_t lsb(uint64_t x, unsigned long& ret) {
    if (x) {
        ret = (unsigned long)__builtin_ctzll(x);
        return true;
    }
    return false;
}

static uint64_t binary_bitsize(uint64_t x) { return util::msbll(x) + 1; }
static uint64_t unary_bitsize(uint64_t x) { return x; }
static uint64_t gamma_bitsize(uint64_t x) {
    uint64_t b = binary_bitsize(x + 1);
    return unary_bitsize(b) + b - 1;
}
static uint64_t delta_bitsize(uint64_t x) {
    uint64_t b = binary_bitsize(x + 1);
    return gamma_bitsize(b - 1) + b - 1;
}

}  // namespace util
}  // namespace fulgor

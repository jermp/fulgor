#pragma once

#include <algorithm>  // for std::set_intersection
#include <array>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <mutex>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

#include "external/sshash/external/cityhash/cityhash.hpp"

#include "external/FQFeeder/include/blockingconcurrentqueue.h"

#include "external/sshash/external/pthash/external/bits/external/essentials/include/essentials.hpp"

#include "include/profiling.hpp"

namespace fulgor {

enum index_t { HYBRID, DIFF, META, META_DIFF };
enum encoding_t { delta_gaps, bitmap, complement_delta_gaps, symmetric_difference };

namespace constants {

constexpr double invalid_threshold = -1.0;
constexpr uint64_t default_ram_limit_in_GiB = 8;
static const std::string default_tmp_dirname(".");
static const std::string hfur_filename_extension(".fur");
static const std::string mfur_filename_extension(".mfur");
static const std::string dfur_filename_extension(".dfur");
static const std::string mdfur_filename_extension(".mdfur");

namespace current_version_number {
constexpr uint8_t major = 4;
constexpr uint8_t minor = 2;
constexpr uint8_t patch = 0;
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
        //
        , meta_colored(false)
        , diff_colored(false)  //
        , fresh_cdbg(false)   //
    {}

    uint32_t k;            // kmer length
    uint32_t m;            // minimizer length
    uint32_t num_threads;  // for building and checking correctness
    uint64_t ram_limit_in_GiB;
    uint64_t num_colors;

    std::filesystem::path tmp_dirname;
    std::filesystem::path output_filename;
    std::filesystem::path filenames_list;

    std::filesystem::path base_filename;
    std::filesystem::path index_filename_to_partition;

    bool verbose;

    bool meta_colored;
    bool diff_colored;
    bool fresh_cdbg;  // ignore a cached cdbg build and reconstruct it

    std::filesystem::path tmp_filename(const std::string& filename) const {
        return tmp_dirname / filename;
    }
};

struct kmer_conservation_triple {
    uint32_t start_pos_in_query;
    uint32_t num_kmers;
    uint32_t color_set_id;
};

typedef uint32_t count_type;

namespace util {
inline void print_cmd(int argc, char** argv) {
    for (int i = 0; i != argc; ++i) std::cout << argv[i] << ' ';
    std::cout << std::endl;
}

inline std::string filename(std::string const& path) {
    return path.substr(path.find_last_of("/\\") + 1);
}

inline void check_version_number(essentials::version_number const& vnum) {
    if (vnum.x != constants::current_version_number::major) {
        throw std::runtime_error("MAJOR index version mismatch: Fulgor index needs rebuilding");
    }
}

/* The index stores counts, ids and offsets as 32-bit values; validate them once at the
   boundary where they enter the system (after the cdbg build, or when loading an index). */
inline void check_index_limits(uint64_t num_colors, uint64_t num_unitigs, uint64_t num_color_sets) {
    std::string exceeded;
    auto check = [&exceeded](uint64_t value, char const* name) {
        if (value > UINT32_MAX) {
            if (!exceeded.empty()) exceeded += ", ";
            exceeded += std::string(name) + "=" + std::to_string(value);
        }
    };
    check(num_colors, "num_colors");
    check(num_unitigs, "num_unitigs");
    check(num_color_sets, "num_color_sets");
    if (!exceeded.empty()) {
        throw std::runtime_error("32-bit index limit exceeded: " + exceeded + " (max " +
                                 std::to_string(UINT32_MAX) + ")");
    }
}

/* The builder appends the version number (3 bytes) at the end of an index file. */
inline essentials::version_number read_version_number(std::string const& index_filename) {
    std::ifstream in(index_filename, std::ios::binary | std::ios::ate);
    if (!in) throw std::runtime_error("cannot open index file: " + index_filename);

    const auto size = static_cast<std::streamoff>(in.tellg());
    if (size < 3) throw std::runtime_error("file too small: " + index_filename);

    in.seekg(size - 3);
    uint8_t b[3];
    in.read(reinterpret_cast<char*>(b), 3);
    if (!in) throw std::runtime_error("cannot read the version number: " + index_filename);

    return essentials::version_number(b[0], b[1], b[2]);
}

inline void check_index_version(std::string const& index_filename) {
    check_version_number(read_version_number(index_filename));
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

inline __uint128_t hash128(char const* bytes, uint64_t num_bytes,
                           const uint64_t seed = 1234567890) {
    auto ret = cityhash::CityHash128WithSeed(bytes, num_bytes, {seed, seed});
    __uint128_t out = 0;
    out += __uint128_t(ret.first);
    out += __uint128_t(ret.second) << 64;
    return out;
}

struct hasher_uint128_t {
    uint64_t operator()(const __uint128_t x) const { return static_cast<uint64_t>(x) ^ (x >> 64); }
};

inline int num_digits(const uint32_t n) {
    if (n >= 10000) {
        if (n >= 10000000) {
            if (n >= 100000000) {
                if (n >= 1000000000) return 10;
                return 9;
            }
            return 8;
        }
        if (n >= 100000) {
            if (n >= 1000000) return 7;
            return 6;
        }
        return 5;
    }
    if (n >= 100) {
        if (n >= 1000) return 4;
        return 3;
    }
    if (n >= 10) return 2;
    return 1;
}

template <typename T>
void vec_to_tsv(std::vector<T> const& vec, std::string& s) {
    s.clear();
    if (vec.empty()) return;
    s.reserve(vec.size() * (std::numeric_limits<T>::digits10 + 2));
    char buffer[std::numeric_limits<T>::digits10 + 3];
    for (T x : vec) {
        auto [ptr, ec] = std::to_chars(buffer, buffer + sizeof(buffer), x);
        if (ec != std::errc()) throw std::runtime_error("vec_to_tsv: conversion failed");
        s.append(buffer, ptr);
        s.push_back('\t');
    }
    s.pop_back();
}

struct compare_first {
    template <typename T>
    bool operator()(const std::pair<uint64_t, T>& a, const std::pair<uint64_t, T>& b) const {
        return a.first > b.first;
    }
};

template <typename It, typename Sentinel = It::sentinel_type>
struct range_view {
    explicit range_view(It it) : _it(it) {}

    It _it;
    It begin() const { return _it; }
    Sentinel end() const { return {}; }
    auto size() const { return _it.size(); }
};

template <typename T, typename Compare = std::less<T>>
class bounded_priority_queue {
public:
    explicit bounded_priority_queue(const size_t capacity, Compare c = Compare())
        : comp(c), max_capacity(capacity) {}

    void set_gauge(prof::gauge* g) { m_gauge = g; }

    void push(T item) {
        std::unique_lock lock(mtx);
        cv_push.wait(lock, [this, &item] {
            return heap.size() < max_capacity || (!heap.empty() && comp(heap.front(), item));
        });

        const uint64_t extra = entry_bytes(item);
        heap.push_back(std::move(item));
        std::push_heap(heap.begin(), heap.end(), comp);
        if (m_gauge) m_gauge->add(static_cast<int64_t>(extra));
    }

    template <typename Predicate>
    bool try_pop_if(T& out_item, Predicate condition) {
        std::lock_guard lock(mtx);
        if (heap.empty() || !condition(heap.front())) {
            return false;
        }
        std::pop_heap(heap.begin(), heap.end(), comp);

        const uint64_t extra = entry_bytes(heap.back());
        out_item = std::move(heap.back());
        heap.pop_back();
        if (m_gauge) m_gauge->add(-static_cast<int64_t>(extra));

        cv_push.notify_all();
        return true;
    }

private:
    static uint64_t entry_bytes(T const& item) {
        uint64_t bytes = sizeof(T);
        if constexpr (requires { item.second.capacity(); }) {
            bytes += item.second.capacity() *
                     sizeof(typename std::remove_reference_t<decltype(item.second)>::value_type);
        }
        return bytes;
    }

    std::vector<T> heap;
    Compare comp;
    std::mutex mtx;
    std::condition_variable cv_push;
    size_t max_capacity;
    prof::gauge* m_gauge = nullptr;
};

/* Removes `path` when the guard goes out of scope: normal return, early return or exception. */
struct file_guard {
    std::filesystem::path path;

    explicit file_guard(std::filesystem::path filename) : path(std::move(filename)) {}

    ~file_guard() {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    file_guard(file_guard const&) = delete;
    file_guard& operator=(file_guard const&) = delete;
    file_guard(file_guard&&) = delete;
    file_guard& operator=(file_guard&&) = delete;
};

class external_saver {
public:
    external_saver() = default;

    explicit external_saver(const std::string& output_filename)
        : output_stream(output_filename, std::ios::binary | std::ios::trunc)
        , saver(output_stream) {}

    explicit external_saver(std::ofstream existing_stream)
        : output_stream(std::move(existing_stream)), saver(output_stream) {}

    void set_stream(std::ofstream out_stream) {
        output_stream = std::move(out_stream);
        saver.emplace(output_stream);
    }

    template <typename T>
    void visit(T const& item) {
        saver->visit(item);
    }

    template <typename T>
    void write(const T& value) {
        static_assert(std::is_trivially_copyable_v<T>,
                      "Type must be trivially copyable for raw binary writes.");

        output_stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    template <typename T>
    void write_vec_data(const std::vector<T>& vec) {
        if (!vec.empty()) {
            static_assert(std::is_trivially_copyable_v<T>,
                          "Vector elements must be trivially copyable for bulk binary writes.");

            output_stream.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(T));
        }
    }

    void seek(const std::streampos& pos) { output_stream.seekp(pos); }
    void seek_end() { output_stream.seekp(0, std::ios::end); }

    std::streampos tell() { return output_stream.tellp(); }

    void append(const std::ifstream& input_stream) {
        if (input_stream.is_open()) {
            output_stream << input_stream.rdbuf();
        }
    }

    void append(std::ifstream& input_stream, const std::streamsize bytes_to_copy) {
        if (!input_stream.is_open() || bytes_to_copy <= 0) {
            return;
        }

        char buffer[4096];
        std::streamsize total_bytes_written = 0;

        while (total_bytes_written < bytes_to_copy && input_stream) {
            const std::streamsize bytes_to_read = std::min(
                static_cast<std::streamsize>(sizeof(buffer)), bytes_to_copy - total_bytes_written);
            input_stream.read(buffer, bytes_to_read);
            const std::streamsize bytes_read = input_stream.gcount();

            if (bytes_read > 0) {
                output_stream.write(buffer, bytes_read);
                total_bytes_written += bytes_read;
            }

            if (bytes_read < bytes_to_read) {
                break;
            }
        }
    }

    void close() { output_stream.close(); }

private:
    std::ofstream output_stream;
    std::optional<essentials::generic_saver> saver;
};

struct timed_phase {
    explicit timed_phase(const std::string& name) : name(name) {
        essentials::logger(name);
        idx = breakdown.size();
        breakdown.emplace_back(name, 0);
        timer.start();
    }

    ~timed_phase() {
        timer.stop();
        breakdown[idx].second = timer.elapsed();
        std::cout << "** " << name << " took " << timer.elapsed() << " seconds / "
                  << timer.elapsed() / 60 << " minutes" << std::endl;
    }

    static void print_breakdown() {
        if (breakdown.empty()) return;

        // Find the maximum length among all phase names
        std::size_t max_len = 0;
        for (const auto& name : breakdown | std::views::keys) {
            max_len = std::max(max_len, name.length());
        }

        std::cout << "--- Execution Breakdown ---\n";
        for (const auto& [name, seconds] : breakdown) {
            double minutes = static_cast<double>(seconds) / 60.0;

            std::cout << std::format("{:<{}} : {:>8} s | {:>8.2f} m\n", name, max_len, seconds,
                                     minutes);
        }
    }

private:
    std::string name;
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    uint64_t idx;

    inline static std::vector<std::pair<std::string, std::uint64_t>> breakdown;
};

constexpr auto make_dna_lut() {
    std::array<bool, 256> lut{};
    for (unsigned char c : {'A', 'C', 'G', 'T', 'U', 'a', 'c', 'g', 't', 'u'}) {
        lut[c] = true;
    }
    return lut;
}
inline constexpr auto dna_lut = make_dna_lut();

inline bool is_dna(std::string_view sv) noexcept {
    return std::ranges::all_of(sv, [](char c) { return dna_lut[static_cast<unsigned char>(c)]; });
}

template <typename Index>
void load_index(Index& index, const std::string& index_filename, bool mmap, const bool verbose) {
    check_index_version(index_filename);

    if (mmap && Index::color_sets_type::type == META_DIFF) {
        std::cerr << "warning: --mmap is not supported for meta-differential indexes; "
                     "loading into RAM instead"
                  << std::endl;
        mmap = false;
    }

    if (mmap) {
        if (verbose) essentials::logger("*** START: memory mapping the index");
        if (essentials::mmap(index, index_filename.c_str()) == 0) {
            throw std::runtime_error("cannot mmap index file: " + index_filename);
        }
        if (verbose) essentials::logger("*** DONE: memory mapping the index");
    } else {
        if (verbose) essentials::logger("*** START: loading the index");
        essentials::load(index, index_filename.c_str());
        if (verbose) essentials::logger("*** DONE: loading the index");
    }

    check_index_limits(index.num_colors(), index.num_unitigs(), index.num_color_sets());
}

}  // namespace util
}  // namespace fulgor

#pragma once

#include <vector>
#include <algorithm>  // for std::set_intersection
#include <condition_variable>

#include "external/smhasher/src/City.h"

#include "external/FQFeeder/include/blockingconcurrentqueue.h"

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
        , check(false)
        //
        , meta_colored(false)
        , diff_colored(false)  //
    {}

    uint32_t k;            // kmer length
    uint32_t m;            // minimizer length
    uint32_t num_threads;  // for building and checking correctness
    uint64_t ram_limit_in_GiB;
    uint64_t num_colors;

    std::filesystem::path tmp_dirname;
    std::filesystem::path output_filename;
    std::filesystem::path filenames_list;

    std::filesystem::path index_filename_to_partition;

    bool verbose;
    bool check;

    bool meta_colored;
    bool diff_colored;

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

void print_cmd(int argc, char** argv) {
    for (int i = 0; i != argc; ++i) std::cout << argv[i] << ' ';
    std::cout << std::endl;
}

std::string filename(std::string const& path) { return path.substr(path.find_last_of("/\\") + 1); }

void check_version_number(essentials::version_number const& vnum) {
    if (vnum.x != constants::current_version_number::major) {
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

inline void vec_to_tsv(std::vector<uint32_t> const& vec, std::string& s) {
    s.clear();
    s.reserve(vec.size() * 12);
    char buffer[32];
    buffer[31] = '\t';
    uint32_t tmp;
    for (uint32_t x : vec) {
        int len = 0;
        do {
            tmp = x / 10;
            buffer[30 - len++] = '0' + (x - tmp * 10);
            x = tmp;
        } while (x > 0);
        s.append(buffer + 31 - len, len + 1);
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

    void push(T item) {
        std::unique_lock lock(mtx);
        cv_push.wait(lock, [this, &item] {
            return heap.size() < max_capacity || (!heap.empty() && comp(heap.front(), item));
        });

        heap.push_back(std::move(item));
        std::push_heap(heap.begin(), heap.end(), comp);
    }

    template <typename Predicate>
    bool try_pop_if(T& out_item, Predicate condition) {
        std::lock_guard lock(mtx);
        if (heap.empty() || !condition(heap.front())) {
            return false;
        }
        std::pop_heap(heap.begin(), heap.end(), comp);

        out_item = std::move(heap.back());
        heap.pop_back();

        cv_push.notify_all();
        return true;
    }

private:
    std::vector<T> heap;
    Compare comp;
    std::mutex mtx;
    std::condition_variable cv_push;
    size_t max_capacity;
};

class external_saver {
public:
    explicit external_saver(const std::string& output_filename)
        : output_stream(output_filename, std::ios::binary | std::ios::trunc)
        , saver(output_stream) {}

    template <typename T>
    void visit(T const& item) {
        saver.visit(item);
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
    essentials::generic_saver saver;
};

}  // namespace util
}  // namespace fulgor

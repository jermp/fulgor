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

#include "external/FQFeeder/include/blockingconcurrentqueue.h"

namespace fulgor {

enum index_t { HYBRID, DIFF, META, META_DIFF };
enum encoding_t { delta_gaps, bitmap, complement_delta_gaps, symmetric_difference };

namespace constants {

constexpr double invalid_threshold = -1.0;
constexpr uint64_t default_ram_limit_in_GiB = 8;
static const std::string default_tmp_dirname(".");
static const std::string hfur_filename_extension("fur");
static const std::string mfur_filename_extension("mfur");
static const std::string dfur_filename_extension("dfur");
static const std::string mdfur_filename_extension("mdfur");

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
    bool check;

    bool meta_colored;
    bool diff_colored;
};

struct kmer_conservation_triple {
    uint32_t start_pos_in_query;
    uint32_t num_kmers;
    uint32_t color_set_id;
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

/*
    Return the largest power of 2 that is <= n.
    Note: could use bitwise tricks for more efficiency.
*/
uint64_t largest_power_of_2(const uint64_t n) {
    if (n == 0) return 0;
    return 1 << static_cast<uint64_t>(std::floor(std::log2(n)));
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

inline int num_digits(const uint32_t n) {
    if (n >= 10000) {
        if (n >= 10000000) {
            if (n >= 100000000) {
                if (n >= 1000000000)
                    return 10;
                return 9;
            }
            return 8;
        }
        if (n >= 100000) {
            if (n >= 1000000)
                return 7;
            return 6;
        }
        return 5;
    }
    if (n >= 100) {
        if (n >= 1000)
            return 4;
        return 3;
    }
    if (n >= 10)
        return 2;
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

template <typename Formatter>
struct formatter_buffer {
    explicit formatter_buffer(Formatter* ptr) : m_formatter(ptr), m_num_bytes(0) {}

    void write(const uint32_t query_id, std::vector<uint32_t> const& vec) {
        m_num_bytes += m_formatter->format(m_buffer, query_id, vec);

        if (m_num_bytes > (1 << 14)) {
            m_formatter->flush(m_buffer, m_num_bytes);
            m_num_bytes = 0;
        }
    }

    ~formatter_buffer() {
        m_formatter->flush(m_buffer, m_num_bytes);
    }

private:
    Formatter* m_formatter;
    typename Formatter::buffer_t m_buffer;
    uint32_t m_num_bytes;
};

struct psa_ascii_formatter {
    typedef std::stringstream buffer_t;

    explicit psa_ascii_formatter(const std::string &output_filename)
        : m_file(output_filename){}

    formatter_buffer<psa_ascii_formatter> buffer() {
        return formatter_buffer(this);
    }

    uint32_t format(buffer_t& ss, uint32_t query_id, std::vector<uint32_t> const& colors) {
        uint32_t num_bytes = 0;
        ss << query_id << "\t" << colors.size();
        num_bytes += num_digits(query_id) + 1 + num_digits(colors.size()); // query_id + tab + size

        if (!colors.empty()) {
            std::string tmpstr;
            vec_to_tsv(colors, tmpstr);
            ss << "\t" << tmpstr;
            num_bytes += 1 + tmpstr.size(); // tab + size of string
        }

        ss << "\n";
        num_bytes += 1; // newline
        return num_bytes;
    }

    void flush(buffer_t& buffer, const uint32_t num_bytes) {
        m_mut.lock();
        m_file.write(buffer.str().data(), num_bytes);
        m_mut.unlock();
        buffer.str("");
    }

protected:
    std::ofstream m_file;
    std::mutex m_mut;
};

struct psa_slow_formatter {
    typedef std::stringstream buffer_t;

    explicit psa_slow_formatter(const std::string &output_filename)
        : m_file(output_filename){}

    formatter_buffer<psa_slow_formatter> buffer() {
        return formatter_buffer(this);
    }

    uint32_t format(buffer_t& ss, uint32_t query_id, std::vector<uint32_t> const& colors) {
        uint32_t num_bytes = 0;
        ss << query_id << "\t" << colors.size();
        num_bytes += num_digits(query_id) + 1 + num_digits(colors.size()); // query_id + tab + size

        for(auto c : colors) {
            ss << '\t' << c;
            num_bytes += 1 + num_digits(c); // tab + color
        }
        ss << "\n";
        num_bytes += 1; // newline

        return num_bytes;
    }

    void flush(buffer_t& buffer, const uint32_t num_bytes) {
        m_mut.lock();
        m_file.write(buffer.str().data(), num_bytes);
        m_mut.unlock();
        buffer.str("");
    }

protected:
    std::ofstream m_file;
    std::mutex m_mut;
};

struct psa_binary_formatter {
    typedef std::stringstream buffer_t;

    explicit psa_binary_formatter(const std::string &output_filename)
        : m_file(output_filename){}

    formatter_buffer<psa_binary_formatter> buffer() {
        return formatter_buffer(this);
    }

    uint32_t format(buffer_t& ss, uint32_t query_id, std::vector<uint32_t> const& colors) {
        uint32_t colors_size = colors.size();
        ss.write(reinterpret_cast<char*>(&query_id), sizeof(uint32_t));
        ss.write(reinterpret_cast<char*>(&colors_size), sizeof(uint32_t));
        if (!colors.empty())
            ss.write(reinterpret_cast<const char*>(colors.data()), colors.size() * sizeof(uint32_t));
        return (2 + colors_size) * sizeof(uint32_t); // (query id + size + colors) * 4 Bytes
    }

    void flush(buffer_t& buffer, const uint32_t num_bytes) {
        m_mut.lock();
        m_file.write(buffer.str().data(), num_bytes);
        m_mut.unlock();
        buffer.str("");
    }

protected:
    std::ofstream m_file;
    std::mutex m_mut;
};

struct psa_compressed_formatter {
    typedef bits::bit_vector::builder buffer_t;

    explicit psa_compressed_formatter(const std::string &output_filename)
        : m_file(output_filename), m_num_colors(0), m_sparse_set_threshold_size(0), m_very_dense_set_threshold_size(0) {}

    void set_num_colors(uint32_t num_colors) {
        assert(m_num_colors == 0);
        m_num_colors = num_colors;
        m_file.write(reinterpret_cast<char*>(&m_num_colors), sizeof(uint64_t));
        m_sparse_set_threshold_size = 0.25 * m_num_colors;
        m_very_dense_set_threshold_size = 0.75 * m_num_colors;
    }

    formatter_buffer<psa_compressed_formatter> buffer() {
        return formatter_buffer(this);
    }

    uint32_t format(buffer_t& bvb, uint32_t query_id, std::vector<uint32_t> const& colors) {
        assert(m_num_colors != 0);
        const uint32_t num_bytes_start = bvb.data().size() * sizeof(uint64_t);
        const uint32_t size = colors.size();
        bits::util::write_delta(bvb, query_id);
        bits::util::write_delta(bvb, size); /* encode size */
        if (size == 0) {

        } else if (size < m_sparse_set_threshold_size) {
            uint32_t prev_val = colors[0];
            bits::util::write_delta(bvb, prev_val);
            for (uint64_t i = 1; i != size; ++i) {
                uint32_t val = colors[i];
                assert(val >= prev_val + 1);
                bits::util::write_delta(bvb, val - (prev_val + 1));
                prev_val = val;
            }
        } else if (size < m_very_dense_set_threshold_size) {
            bits::bit_vector::builder tmp;
            tmp.resize(m_num_colors);
            for (uint64_t i = 0; i != size; ++i) tmp.set(colors[i]);
            bvb.append(tmp);
        } else {
            bool first = true;
            uint32_t val = 0;
            uint32_t prev_val = -1;
            uint32_t written = 0;
            for (uint64_t i = 0; i != size; ++i) {
                uint32_t x = colors[i];
                while (val < x) {
                    if (first) {
                        bits::util::write_delta(bvb, val);
                        first = false;
                        ++written;
                    } else {
                        assert(val >= prev_val + 1);
                        bits::util::write_delta(bvb, val - (prev_val + 1));
                        ++written;
                    }
                    prev_val = val;
                    ++val;
                }
                assert(val == x);
                val++;  // skip x
            }
            while (val < m_num_colors) {
                assert(val >= prev_val + 1);
                bits::util::write_delta(bvb, val - (prev_val + 1));
                prev_val = val;
                ++val;
                ++written;
            }
            assert(val == m_num_colors);
            /* complementary_set_size = m_num_colors - size */
            assert(m_num_colors - size <= m_num_colors);
            assert(written == m_num_colors - size);
        }
        // ss.write(reinterpret_cast<const char*>(bvb.data().data()), bvb.data().size() * sizeof(uint64_t));
        return bvb.data().size() * sizeof(uint64_t) - num_bytes_start;
    }

    void flush(buffer_t& buffer, const uint32_t num_bytes) {
        m_mut.lock();
        const uint64_t num_bits = buffer.num_bits();
        m_file.write(reinterpret_cast<const char*>(&num_bits), sizeof(uint64_t));
        m_file.write(reinterpret_cast<const char*>(buffer.data().data()), num_bytes);
        m_mut.unlock();
        buffer.clear();
    }

protected:
    std::ofstream m_file;
    std::mutex m_mut;
    uint32_t m_num_colors, m_sparse_set_threshold_size, m_very_dense_set_threshold_size;
};

template <typename FulgorIndex>
struct fastq_query_reader {
    fastq_query_reader(std::string& query_filename, uint64_t num_threads, FulgorIndex& index)
        : rparser({query_filename}, num_threads, num_threads-1)
        , index(index){
        rparser.start();
    }

    struct query_t {
        query_t() : id(-1) {}
        query_t(uint32_t id, uint32_t size) :
            id(id) {
            cids.resize(size);
        }

        query_t(uint32_t id, std::vector<uint32_t>&& cids_) :
            id(id), cids(std::move(cids_)) {}

        uint32_t id;
        std::vector<uint32_t> cids;
        std::string seq;
    };

    struct query_group {
        query_group(fastq_query_reader* qb_)
            : qb(qb_), rg(qb->rparser.getReadGroup()) {}

        bool has_next() {
            return curr_record != rg.end();
        }

        void next() {
            ++curr_record;
            ++curr_read_id;
        }
        void operator++() { next(); }

        void value(query_t& query) {
            query.id = curr_read_id;
            query.cids.clear();
            qb->index.fetch_color_set_ids(curr_record->seq, query.cids);
            query.seq = curr_record->seq;
        }

        bool refill() {
            const bool result = qb->rparser.refill(rg);
            if (result) {
                curr_record = rg.begin();
                curr_read_id = rg.chunk_frag_offset().frag_idx;
            }
            return result;
        }

    private:
        fastq_query_reader* qb;
        fastx_parser::ReadGroup<klibpp::KSeq> rg;
        vector<klibpp::KSeq>::iterator curr_record;
        uint32_t curr_read_id;
    };

    query_group get_query_group() {
        return query_group(this);
    }

    ~fastq_query_reader() {
        rparser.stop();
    }

private:
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser;
    FulgorIndex& index;
};

struct preprocessed_query_reader {
    struct query_t {
        query_t(): ids(0), cids(0) {}
        query_t(std::vector<uint32_t>&& ids_, std::vector<uint32_t>&& cids_) :
            ids(std::move(ids_)), cids(std::move(cids_)) {}

        void clear() {
            ids.clear();
            cids.clear();
        }

        std::vector<uint32_t> ids;
        std::vector<uint32_t> cids;
        std::string seq;
    };

    preprocessed_query_reader(std::string const& query_filename, uint64_t num_threads)
        : m_query_file(query_filename), m_queue(3*num_threads)
    {
        const uint64_t batch_thresh = 10000;

        m_producer = std::thread([this] () -> void {
            auto* batch = new std::vector<query_t>();
            batch->reserve(batch_thresh);
            uint32_t list_len = 0, query_id;
            size_t nbatches = 0;

            query_t query;
            this->m_query_file.read(reinterpret_cast<char*>(&list_len), sizeof(list_len));
            this->m_query_file.read(reinterpret_cast<char*>(&query_id), sizeof(query_id));
            query.ids.push_back(query_id);
            query.cids.resize(list_len-1);
            this->m_query_file.read(reinterpret_cast<char*>(query.cids.data()),
                               (list_len - 1) * sizeof(list_len));

            while (this->m_query_file.read(reinterpret_cast<char*>(&list_len), sizeof(list_len))) {
                assert(list_len > 0);
                this->m_query_file.read(reinterpret_cast<char*>(&query_id), sizeof(query_id));

                if (list_len == 1) {
                    query.ids.push_back(query_id);
                    continue;
                }
                batch->push_back(query);

                if (batch->size() >= 10000) {
                    while (!this->m_queue.try_enqueue(batch)) {}
                    batch = new std::vector<query_t>();
                    this->in_flight += 1;
                    ++nbatches;
                }

                query.clear();
                query.ids.push_back(query_id);
                query.cids.resize(list_len-1);
                this->m_query_file.read(reinterpret_cast<char*>(query.cids.data()),
                           (list_len - 1) * sizeof(uint32_t));
            }
            batch->push_back(query);


            if (!batch->empty()) {
                this->m_queue.enqueue(batch);
                this->in_flight += 1;
            }
            this->done = true;
        });
    }

    struct query_group {
        query_group(preprocessed_query_reader* qb_)
            : reader(qb_), sbatch(nullptr) {}

        bool has_next() {
            return sbatch != nullptr && curr_query != sbatch->end();
        }

        void next() {
            ++curr_query;
        }
        void operator++() { next(); }

        void value(query_t& query) {
            query.ids = std::move(curr_query->ids);
            query.cids = std::move(curr_query->cids);
        }

        bool refill() {
            std::lock_guard lock(reader->mut);
            while (!reader->done or (reader->in_flight > 0)) {
                if (reader->m_queue.try_dequeue(sbatch)) break;
            }
            if (sbatch != nullptr) {
                curr_query = sbatch->begin();
            }
            if (reader->in_flight > 0) {
                reader->in_flight -= 1;
                return true;
            }
            return !reader->done;
        }

    private:
        preprocessed_query_reader* reader;
        std::vector<query_t>* sbatch;
        std::vector<query_t>::iterator curr_query;
    };

    query_group get_query_group() {
        return query_group(this);
    }

    ~preprocessed_query_reader() {
        m_producer.join();
    }

private:
    std::ifstream m_query_file;
    moodycamel::BlockingConcurrentQueue<std::vector<query_t>*> m_queue;
    std::thread m_producer;
    std::atomic<bool> done{false};
    std::atomic<uint32_t> in_flight{0};
    std::mutex mut;


};

}  // namespace util
}  // namespace fulgor

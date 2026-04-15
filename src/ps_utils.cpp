#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "include/index.hpp"
#include "external/FQFeeder/include/FastxParser.hpp"
#include "external/FQFeeder/include/blockingconcurrentqueue.h"

namespace fulgor {

enum class pseudoalignment_algorithm : uint8_t { FULL_INTERSECTION, THRESHOLD_UNION };

std::string to_string(const pseudoalignment_algorithm algo, const double threshold) {
    std::string o;
    switch (algo) {
        case pseudoalignment_algorithm::FULL_INTERSECTION:
            o = "full-intersection";
            break;
        case pseudoalignment_algorithm::THRESHOLD_UNION:
            o = "threshold-union (threshold = " + std::to_string(threshold) + ")";
            break;
    }
    return o;
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
        num_bytes += util::num_digits(query_id) + 1 + util::num_digits(colors.size()); // query_id + tab + size

        if (!colors.empty()) {
            std::string tmpstr;
            util::vec_to_tsv(colors, tmpstr);
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
        num_bytes += util::num_digits(query_id) + 1 + util::num_digits(colors.size()); // query_id + tab + size

        for(auto c : colors) {
            ss << '\t' << c;
            num_bytes += 1 + util::num_digits(c); // tab + color
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
        explicit query_group(fastq_query_reader* qb_)
            : qb(qb_), rg(qb->rparser.getReadGroup()), curr_read_id(0) {}

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
        std::vector<klibpp::KSeq>::iterator curr_record;
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
         constexpr uint64_t batch_thresh = 10000;

        m_producer = std::thread([this] () -> void {
            auto* batch = new std::vector<query_t>();
            batch->reserve(batch_thresh);
            uint32_t list_len = 0, query_id;

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
        explicit query_group(preprocessed_query_reader* qb_)
            : reader(qb_), sbatch(nullptr) {}

        bool has_next() const {
            return sbatch != nullptr && curr_query != sbatch->end();
        }

        void next() {
            ++curr_query;
        }
        void operator++() { next(); }

        void value(query_t& query) const {
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

struct ps_options {
    explicit ps_options(const pseudoalignment_algorithm algo, const bool verbose, const uint64_t num_threads)
        : verbose(verbose)
        , algo(algo)
        , num_threads(num_threads)
        , num_reads(0)
        , num_mapped_reads(0) {}

    void increment_processed_reads(const int val = 1) {
        uint64_t prev = num_reads.fetch_add(val);
        if (verbose && prev >> batch_size != (prev + val) >> batch_size) {
            io_mut.lock();
            std::cout << "processed " << num_reads << " reads" << std::endl;
            io_mut.unlock();
        }
    }

    void increment_mapped_reads(const int val = 1) {
        num_mapped_reads += val;
    }

    const bool verbose;
    const pseudoalignment_algorithm algo;
    const uint64_t num_threads;
    std::mutex io_mut;
    std::atomic<uint64_t> num_reads;
    std::atomic<uint64_t> num_mapped_reads;

private:
    const uint64_t batch_size = 20;
};

}
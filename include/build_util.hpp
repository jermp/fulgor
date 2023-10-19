#pragma once

#include "sketch/hll.h"
#include "kmeans.hpp"
#include "../external/FQFeeder/include/FastxParser.hpp"
#include "../external/FQFeeder/src/FastxParser.cpp"
#include "../external/sshash/include/query/streaming_query_canonical_parsing.hpp"

namespace fulgor {

void exe(sshash::dictionary const& k2u, std::vector<sketch::hll_t>& sketches, uint64_t i,
         std::string filename) {
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser({filename}, 1);
    sshash::streaming_query_canonical_parsing query(&k2u);
    parser.start();
    auto chunk = parser.getReadGroup();
    while (parser.refill(chunk)) {
        for (auto const& record : chunk) {
            query.start();
            const uint64_t num_kmers = record.seq.length() - k2u.k() + 1;
            for (uint64_t j = 0, prev_unitig_id = -1; j != num_kmers; ++j) {
                char const* kmer = record.seq.data() + j;
                auto answer = query.lookup_advanced(kmer);
                // kmer should always be found actually...
                if (answer.kmer_id != sshash::constants::invalid_uint64) {
                    if (answer.contig_id != prev_unitig_id) {
                        uint32_t unitig_id = answer.contig_id;
                        sketches[i].addh(unitig_id);
                        prev_unitig_id = unitig_id;
                    }
                }
            }
        }
    }
    parser.stop();
}

void build_reference_sketches(index_type const& index,
                              uint64_t p,  // use 2^p bytes per HLL sketch
                              uint64_t num_threads,
                              std::string output_filename  // where the sketches will be serialized
) {
    uint64_t num_docs = index.num_docs();
    std::vector<sketch::hll_t> sketches(num_docs, sketch::hll_t(p));
    typename sketch::hll_t::HashType hasher;

    sshash::dictionary const& k2u = index.get_dict();

    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    // auto exe = [&](uint64_t i, std::string filename) {
    //     fastx_parser::FastxParser<fastx_parser::ReadSeq> parser({filename}, 1);
    //     sshash::streaming_query_canonical_parsing query(&k2u);
    //     parser.start();
    //     auto chunk = parser.getReadGroup();
    //     while (parser.refill(chunk)) {
    //         for (auto const& record : chunk) {
    //             query.start();
    //             const uint64_t num_kmers = record.seq.length() - k2u.k() + 1;
    //             for (uint64_t j = 0, prev_unitig_id = -1; j != num_kmers; ++j) {
    //                 char const* kmer = record.seq.data() + j;
    //                 auto answer = query.lookup_advanced(kmer);
    //                 // kmer should always be found actually...
    //                 if (answer.kmer_id != sshash::constants::invalid_uint64) {
    //                     if (answer.contig_id != prev_unitig_id) {
    //                         uint32_t unitig_id = answer.contig_id;
    //                         sketches[i].addh(unitig_id);
    //                         prev_unitig_id = unitig_id;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     parser.stop();
    // };

    for (uint64_t i = 0; i != num_docs;) {
        for (uint64_t j = 0; j != num_threads and i != num_docs; ++j, ++i) {
            std::string filename(index.filename(i));
            threads.push_back(
                std::thread([&k2u, &sketches, i, filename]() { exe(k2u, sketches, i, filename); }));
        }
        for (auto& t : threads) {
            if (t.joinable()) t.join();
        }
    }

    std::ofstream out(output_filename, std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("cannot open file");
    uint64_t num_bytes = 1ULL << p;
    out.write(reinterpret_cast<char const*>(&num_bytes), 8);
    out.write(reinterpret_cast<char const*>(&num_docs), 8);
    for (auto const& x : sketches) {
        assert(x.m() == num_bytes);
        assert(x.m() == x.core().size());
        uint8_t const* data = x.data();
        out.write(reinterpret_cast<char const*>(data), num_bytes);
    }
    out.close();
}

// void build_reference_sketches(index_type const& index,
//                               uint64_t p,                  // use 2^p bytes per HLL sketch
//                               std::string output_filename  // where the sketches will be
//                               serialized
// ) {
//     uint64_t num_docs = index.num_docs();
//     std::vector<sketch::hll_t> sketches(num_docs, sketch::hll_t(p));
//     typename sketch::hll_t::HashType hasher;

//     auto const& u2c = index.get_u2c();
//     pthash::bit_vector::unary_iterator unary_it(u2c);
//     auto const& ccs = index.color_classes();
//     uint64_t num_color_classes = ccs.num_color_classes();
//     uint64_t num_ones = u2c.num_ones();
//     uint64_t pop_count = 0;
//     uint64_t prev_pos = 0;
//     std::vector<uint64_t> hashes;
//     for (uint64_t color_id = 0; color_id != num_color_classes; ++color_id) {
//         uint64_t curr_pos = pop_count != num_ones ? unary_it.next() : (u2c.size() - 1);
//         auto it = ccs.colors(color_id);
//         uint64_t size = it.size();
//         hashes.reserve(curr_pos - prev_pos + 1);
//         for (uint64_t unitig_id = prev_pos; unitig_id <= curr_pos; ++unitig_id) {
//             assert(index.u2c(unitig_id) == color_id);
//             hashes.push_back(hasher.hash(unitig_id));
//         }
//         for (uint64_t i = 0; i != size; ++i, ++it) {
//             uint32_t ref_id = *it;
//             assert(ref_id < num_docs);
//             for (auto hash : hashes) sketches[ref_id].add(hash);
//         }
//         pop_count += 1;
//         prev_pos = curr_pos + 1;
//         hashes.clear();
//     }
//     assert(prev_pos == index.get_dict().num_contigs());

//     std::ofstream out(output_filename, std::ios::binary);
//     if (!out.is_open()) throw std::runtime_error("cannot open file");
//     uint64_t num_bytes = 1ULL << p;
//     out.write(reinterpret_cast<char const*>(&num_bytes), 8);
//     out.write(reinterpret_cast<char const*>(&num_docs), 8);
//     for (auto const& x : sketches) {
//         assert(x.m() == num_bytes);
//         assert(x.m() == x.core().size());
//         uint8_t const* data = x.data();
//         out.write(reinterpret_cast<char const*>(data), num_bytes);
//     }
//     out.close();
// }

}  // namespace fulgor
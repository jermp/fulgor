#include "sketch/hll.h"

using namespace fulgor;

struct sketch_data {
    sketch_data(uint64_t p) : set_size(0), sketch(p) {}
    uint32_t set_size;
    sketch::hll_t sketch;
};

int sketch_references(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The fulgor index filename.", "-i", true);
    parser.add("output_filename", "Output filename for the sketches.", "-o", true);
    parser.add("p", "Use 2^p bytes for each HLL sketc.", "-p", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    auto index_filename = parser.get<std::string>("index_filename");
    index_type index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");

    uint64_t num_docs = index.num_docs();
    uint64_t p = parser.get<uint64_t>("p");  // use 2^p bytes per HLL sketch
    std::vector<sketch_data> sketches(num_docs, sketch_data(p));
    typename sketch::hll_t::HashType hasher;

    essentials::logger("building sketches...");

    essentials::timer_type timer;
    timer.start();

    auto const& u2c = index.get_u2c();
    pthash::bit_vector::unary_iterator unary_it(u2c);
    auto const& ccs = index.color_classes();
    uint64_t num_color_classes = ccs.num_color_classes();
    uint64_t num_ones = u2c.num_ones();
    uint64_t pop_count = 0;
    uint64_t prev_pos = 0;
    std::vector<uint64_t> hashes;
    for (uint64_t color_id = 0; color_id != num_color_classes; ++color_id) {
        uint64_t curr_pos = pop_count != num_ones ? unary_it.next() : (u2c.size() - 1);
        auto it = ccs.colors(color_id);
        uint64_t size = it.size();
        hashes.reserve(curr_pos - prev_pos + 1);
        for (uint64_t unitig_id = prev_pos; unitig_id <= curr_pos; ++unitig_id) {
            assert(index.u2c(unitig_id) == color_id);
            hashes.push_back(hasher.hash(unitig_id));
        }
        for (uint64_t i = 0; i != size; ++i, ++it) {
            uint32_t ref_id = *it;
            assert(ref_id < num_docs);
            sketches[ref_id].set_size += hashes.size();
            for (auto hash : hashes) sketches[ref_id].sketch.add(hash);
        }
        pop_count += 1;
        prev_pos = curr_pos + 1;
        hashes.clear();
    }
    assert(prev_pos == index.get_dict().num_contigs());

    timer.stop();
    std::cout << "computing all sketches took: " << timer.elapsed() / 1000000 << " [sec] ("
              << (timer.elapsed() / 1000) / num_docs << " ms/sketch)" << std::endl;

    essentials::logger("writing sketches to disk...");
    auto output_filename = parser.get<std::string>("output_filename");
    std::ofstream out(output_filename, std::ios::binary);
    if (!out.is_open()) return 1;
    uint64_t num_bytes = 1ULL << p;
    out.write(reinterpret_cast<char const*>(&num_bytes), 4);
    out.write(reinterpret_cast<char const*>(&num_docs), 4);
    for (auto const& x : sketches) {
        assert(x.sketch.m() == num_bytes);
        assert(x.sketch.m() == x.sketch.core().size());
        uint8_t const* data = x.sketch.data();
        out.write(reinterpret_cast<char const*>(&x.set_size), sizeof(x.set_size));
        out.write(reinterpret_cast<char const*>(data), num_bytes);
    }
    out.close();
    essentials::logger("DONE");

    return 0;
}
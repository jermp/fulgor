#include "build_util.hpp"

using namespace fulgor;

int sketch_references(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("output_filename", "Output filename for the sketches.", "-o", true);
    parser.add("p", "Use 2^p bytes for each HLL sketch.", "-p", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    auto index_filename = parser.get<std::string>("index_filename");
    index_type index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");

    essentials::logger("building and writing sketches to disk...");

    essentials::timer_type timer;
    timer.start();

    build_reference_sketches(index, parser.get<uint64_t>("p"),
                             parser.get<std::string>("output_filename"));

    timer.stop();
    essentials::logger("DONE");
    std::cout << timer.elapsed() / 1000000 << " [sec] ("
              << (timer.elapsed() / 1000) / index.num_docs() << " ms/sketch)" << std::endl;

    return 0;
}

int sketch_colors(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("output_filename", "Output filename for the sketches.", "-o", true);
    parser.add("p", "Use 2^p bytes for each HLL sketch.", "-p", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    auto index_filename = parser.get<std::string>("index_filename");
    index_type index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");

    uint64_t num_color_classes = index.num_color_classes();
    uint64_t p = parser.get<uint64_t>("p");  // use 2^p bytes per HLL sketch
    std::vector<sketch::hll_t> sketches(num_color_classes, sketch::hll_t(p));
    typename sketch::hll_t::HashType hasher;

    essentials::logger("building sketches...");

    essentials::timer_type timer;
    timer.start();

    auto const& ccs = index.color_classes();

#pragma omp parallel for
    for (uint64_t color_id = 0; color_id != num_color_classes; ++color_id) {
        auto& sketch = sketches[color_id];
        auto it = ccs.colors(color_id);
        uint64_t size = it.size();
        for (uint64_t i = 0; i != size; ++i, ++it) {
            uint32_t ref_id = *it;
            assert(ref_id < index.num_docs());
            uint64_t hash = hasher.hash(ref_id);
            sketch.add(hash);
        }
    }
    std::cout << "processed " << num_color_classes << " colors" << std::endl;

    timer.stop();
    std::cout << "computing all sketches took: " << timer.elapsed() / 1000000 << " [sec] ("
              << (timer.elapsed() / 1000) / num_color_classes << " ms/sketch)" << std::endl;

    essentials::logger("writing sketches to disk...");
    auto output_filename = parser.get<std::string>("output_filename");
    std::ofstream out(output_filename, std::ios::binary);
    if (!out.is_open()) return 1;
    uint64_t num_bytes = 1ULL << p;
    out.write(reinterpret_cast<char const*>(&num_bytes), 8);
    out.write(reinterpret_cast<char const*>(&num_color_classes), 8);
    for (auto const& x : sketches) {
        assert(x.m() == num_bytes);
        assert(x.m() == x.core().size());
        uint8_t const* data = x.data();
        out.write(reinterpret_cast<char const*>(data), num_bytes);
    }
    out.close();
    essentials::logger("DONE");

    return 0;
}
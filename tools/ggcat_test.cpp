#include <iostream>
#include <string.h>
#include <fstream>

#include "../external/ggcat/crates/capi/ggcat-cpp-api/include/ggcat.hh"

using namespace ggcat;

int main(int argc, char const* argv[]) {
    size_t k = 31;

    // AAAAACACACATATACAGTGTGTGAGTAGTATGATGTGCAACCGAT
    //                 AGTGTGTGAGTAGTATGATGTGCAACCGATTTTCAG
    const char* sequences[] = {"AAAAACACACATATACAGTGTGTGAGTAGTATGATGTGCAACCGAT",
                               "AGTGTGTGAGTAGTATGATGTGCAACCGATTTTCAG"};

    std::vector<std::string> input_files{"ref1.fa", "ref2.fa"};

    std::ofstream out0((input_files[0]).c_str());
    out0 << ">\n";
    out0.write(sequences[0], strlen(sequences[0]));
    out0 << '\n';
    out0.close();
    std::ofstream out1((input_files[1]).c_str());
    out1 << ">\n";
    out1.write(sequences[1], strlen(sequences[1]));
    out1 << '\n';
    out1.close();

    std::vector<std::string> color_names{"ref1", "ref2"};

    GGCATConfig config;

    config.use_temp_dir = true;
    config.temp_dir = "tmp_dir";
    config.memory = 2.0;
    config.prefer_memory = true;
    config.total_threads_count = 1;
    config.intermediate_compression_level = -1;
    config.use_stats_file = false;
    config.stats_file = "";

    GGCATInstance* instance = GGCATInstance::create(config);

    std::string graph_file("ccdgb.fa");

    std::string output_file = instance->build_graph_from_files(
        Slice<std::string>(input_files.data(), input_files.size()), graph_file, k, 1, false, 1,
        ExtraElaborationStep_UnitigLinks, true,  // output_colors
        Slice<std::string>(color_names.data(), color_names.size()));

    std::vector<std::string> file_color_names =
        GGCATInstance::dump_colors(GGCATInstance::get_colormap_file(graph_file));

    instance->dump_unitigs(
        graph_file, k, 1, true,
        [&](Slice<char> unitig, Slice<uint32_t> color, bool same_color) {
            std::cout << "unitig: '" << std::string(unitig.data, unitig.data + unitig.size) << "'"
                      << std::endl;
            std::cout << "color: [ ";
            for (uint64_t i = 0; i != color.size; ++i) { std::cout << color.data[i] << ' '; }
            std::cout << "]" << std::endl;
            std::cout << "refs: [ ";
            for (uint64_t i = 0; i != color.size; ++i) {
                std::cout << file_color_names[color.data[i]] << ' ';
            }
            std::cout << "]" << std::endl;
        },
        true  // output_colors
    );

    std::remove("ref1.fa");
    std::remove("ref2.fa");
    std::remove("ccdgb.fa");
    std::remove("ccdgb.colors.dat");

    return 0;
}

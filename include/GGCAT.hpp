#pragma once

#include <vector>

#include "../external/ggcat/crates/capi/ggcat-cpp-api/include/ggcat.hh"

namespace fulgor {

struct GGCAT {
    GGCAT(std::vector<std::string>& input_filenames, uint64_t mem_gigas, uint64_t k,
          uint64_t n_threads, std::string const& tmp_dirname)
    // : k(k)
    {
        ggcat::GGCATConfig config;
        config.use_temp_dir = true;
        config.temp_dir = tmp_dirname;
        config.memory = mem_gigas;
        config.prefer_memory = true;
        config.total_threads_count = n_threads;
        config.intermediate_compression_level = -1;

        config.use_stats_file = false;
        config.stats_file = "";

        instance = ggcat::GGCATInstance::create(config);

        graph_file = tmp_dirname + "/ccdbg.ggcat.fa";

        std::vector<std::string> color_names;
        color_names.reserve(filenames.size());
        for (uint64_t i = 0; i != filenames.size(); ++i) {
            color_names.push_back(std::to_string(i));
        }

        constexpr bool forward_only = false;
        constexpr bool colors = true;
        constexpr size_t min_multiplicity = 1;
        // std::string output_file =
        instance->build_graph_from_files(
            ggcat::Slice<std::string>(input_filenames.data(), input_filenames.size()), graph_file,
            k, n_threads, forward_only, min_multiplicity, ggcat::ExtraElaborationStep_UnitigLinks,
            colors, ggcat::Slice<std::string>(color_names.data(), color_names.size()));

        // ideally they should be the same as those in input
        filenames =
            ggcat::GGCATInstance::dump_colors(ggcat::GGCATInstance::get_colormap_file(graph_file));
    }

    // // The callback takes a unitig, the color set, and the is_same flag
    // void iterate(std::function<void(const std::string&, const vector<int64_t>&, bool)> callback)
    // {
    //     vector<int64_t> prev_colors;
    //     std::mutex callback_mutex;

    //     auto outer_callback = [&](Slice<char> read, Slice<uint32_t> colors, bool same_colors) {
    //         // Calls in callback provided by the caller of iterate.
    //         // WARNING: this function is called asynchronously from multiple threads, so it must
    //         be
    //         // thread-safe. Also the same_colors boolean is referred to the previous call of this
    //         // function from the current thread.
    //         std::lock_guard<std::mutex> _lock(callback_mutex);
    //         try {
    //             string unitig = string(read.data, read.data + read.size);
    //             if (same_colors) {
    //                 callback(unitig, prev_colors, true);
    //             } else {
    //                 prev_colors.clear();
    //                 for (size_t i = 0; i < colors.size; i++) {
    //                     prev_colors.push_back(colors.data[i]);
    //                 }
    //                 callback(unitig, prev_colors, false);
    //             }
    //         } catch (const std::exception& e) {
    //             std::cerr << "Caught Error: " << e.what() << '\n';
    //             exit(1);
    //         }
    //     };

    //     this->instance->dump_unitigs(graph_file, k, 1, true, outer_callback, true, -1);
    // }

    std::string get_unitig_filename() const { return graph_file; }

private:
    ggcat::GGCATInstance* instance;
    std::string graph_file;
    std::vector<std::string> filenames;
    // uint64_t k;
};

}  // namespace fulgor
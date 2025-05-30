#pragma once

#include <functional>
#include <vector>

#include "external/ggcat/crates/capi/ggcat-cpp-api/include/ggcat.hh"
#include "util.hpp"

namespace fulgor {

struct GGCAT {
    GGCAT() : m_k(0), m_instance(nullptr) {}

    ~GGCAT() {
        try {
            /* remove GGCAT's tmp files */
            std::remove((m_graph_file).c_str());
            std::remove((m_color_sets_file).c_str());
        } catch (std::exception const& e) { std::cerr << e.what() << std::endl; }
    }

    void build(build_configuration const& build_config) {
        {
            std::ifstream in(build_config.filenames_list);
            if (!in.is_open()) throw std::runtime_error("error in opening file");
            std::string filename;
            while (in >> filename) m_filenames.push_back(filename);
            std::cout << "about to process " << m_filenames.size() << " files..." << std::endl;
            in.close();
        }

        std::string basename =
            build_config.tmp_dirname + "/" + util::filename(build_config.file_base_name);
        m_color_sets_file = basename + ".ggcat.color_sets.dat";
        m_graph_file = basename + ".ggcat.fa";
        m_k = build_config.k;

        ggcat::GGCATConfig config;
        config.use_temp_dir = true;
        config.temp_dir = build_config.tmp_dirname;
        config.memory = build_config.ram_limit_in_GiB;
        config.prefer_memory = true;
        config.total_threads_count = build_config.num_threads;
        config.intermediate_compression_level = -1;
        config.use_stats_file = false;
        config.stats_file = "";

        // GGCAT bug:
        // This leaks memory (not much) but it can't be easily fixed because
        // this memory is allocated in the Rust API of GGCAT and freed only at the
        // end of the program.
        m_instance = ggcat::GGCATInstance::create(config);

        std::vector<std::string> color_names;
        color_names.reserve(m_filenames.size());
        for (uint64_t i = 0; i != m_filenames.size(); ++i) {
            color_names.push_back(std::to_string(i));
        }

        if (std::filesystem::exists(m_graph_file)){
            cout << "GGCAT files found, skipped GGCAT construction" << endl;
            return;
        }

        constexpr bool forward_only = false;
        constexpr bool output_color_sets = true;
        constexpr size_t min_multiplicity = 1;
        m_instance->build_graph_from_files(
            ggcat::Slice<std::string>(m_filenames.data(), m_filenames.size()), m_graph_file, m_k,
            build_config.num_threads, forward_only, min_multiplicity,
            ggcat::ExtraElaborationStep_UnitigLinks, output_color_sets,
            ggcat::Slice<std::string>(color_names.data(), color_names.size()));
    }

    void loop_through_unitigs(std::function<void(ggcat::Slice<char> const /* unitig */,         //
                                                 ggcat::Slice<uint32_t> const /* color_set */,  //
                                                 bool /* same_color_set */)>
                                  callback,
                              uint64_t num_threads = 1) const  //
    {
        if (m_k == 0) throw std::runtime_error("graph must be built first");
        m_instance->dump_unitigs(m_graph_file, m_k, num_threads, num_threads == 1, callback, true);
    }

    uint64_t num_colors() const { return m_filenames.size(); }
    std::vector<std::string> const& filenames() const { return m_filenames; }

private:
    uint64_t m_k;
    ggcat::GGCATInstance* m_instance;
    std::vector<std::string> m_filenames;
    std::string m_graph_file;
    std::string m_color_sets_file;
};

}  // namespace fulgor

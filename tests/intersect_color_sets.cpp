#include <iostream>
#include <fstream>
#include <sstream>

#include "external/CLI11.hpp"
#include "external/sshash/include/gz/zip_stream.hpp"

#include "tests/utils.cpp"

using namespace fulgor;

template <typename FulgorIndex>
int intersect_color_sets(std::string const& index_filename, uint64_t size, std::string const& algo,
                         const bool quiet) {
    FulgorIndex index;
    if (!quiet) essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    if (!quiet) essentials::logger("DONE");

    uint32_t num_intersections = 1000000;
    double density_threshold = 0.25;
    uint32_t n = size;
    uint32_t num_colors = index.num_colors();
    srand(42);

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> t;

    std::vector<std::vector<uint32_t>> color_sets_by_density(std::ceil(1/density_threshold));

    for(uint32_t color_set_id = 0; color_set_id < index.num_color_sets(); ++color_set_id){
        auto it = index.color_set(color_set_id);
        double density = 1.0 * it.size() / num_colors;
        color_sets_by_density[std::ceil(density/density_threshold)-1].push_back(color_set_id);  // std::floor throws segfault if density == 1
    }

    for (uint32_t density_group = 0; density_group < color_sets_by_density.size(); density_group++){
        std::cout << color_sets_by_density[density_group].size() << " ";
    }
    std::cout << endl;

    for (uint32_t density_group = 0; density_group < color_sets_by_density.size(); density_group++){
        uint16_t group_size = color_sets_by_density[density_group].size();
        t.start();
        for (uint32_t i = 0; i < num_intersections; i++) {
            vector<typename FulgorIndex::color_sets_type::iterator_type> iterators;
            for (uint32_t j = 0; j < n; j++) {
                uint32_t id = color_sets_by_density[density_group][rand() % group_size];
                iterators.push_back(index.color_set(id));
            }
            vector<uint32_t> colors;

            if (algo == "geq") {
                next_geq_intersect(iterators, colors, index.num_colors());
            } else {
                counting_intersect(iterators, colors, index.num_colors());
            }
        }
        t.stop();

        if (!quiet) {
            std::cout << "intersected " << num_intersections * size << " color_sets with density [" 
                      << density_threshold*density_group << ", " << density_threshold*(density_group + 1) 
                      << "]" << std::endl;
            std::cout << "elapsed = " << t.elapsed() << " millisec / "
                      << t.elapsed() / 1000 << " sec / "
                      << t.elapsed() / 1000 / 60 << " min / "
                      << (t.elapsed() * 1000) / num_intersections << " musec/intersection" << std::endl;
        }
        t.reset();
    }


    if (!quiet) essentials::logger("DONE");

    

    return 0;
}

int intersect_color_sets(int argc, char** argv) {
    /* params */
    std::string index_filename;
    std::string algo;
    uint64_t size = 2;
    bool quiet = false;

    CLI::App app{"Perform (color-only) pseudoalignment to a Fulgor index."};
    app.add_option("-i,--index", index_filename, "The Fulgor index filename,")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-s,--size", size, "Number of color sets to intersect.")->required();
    app.add_option("-a,--algorithm", algo, "Intersection algorithm ['geq' or 'count'].")
        ->required();
    app.add_flag("--quiet", quiet, "Quiet mode: do not print status messages to stdout.");
    CLI11_PARSE(app, argc, argv);

    if (!quiet) util::print_cmd(argc, argv);

    if (algo != "geq" && algo != "count") {
        std::cerr << "Wrong intersection algorithm" << std::endl;
        return 1;
    }

    if (sshash::util::ends_with(index_filename,
                                constants::meta_diff_colored_fulgor_filename_extension)) {
        return intersect_color_sets<meta_differential_index_type>(index_filename, size, algo,
                                                                  quiet);
    } else if (sshash::util::ends_with(index_filename,
                                       constants::meta_colored_fulgor_filename_extension)) {
        return intersect_color_sets<meta_index_type>(index_filename, size, algo, quiet);
    } else if (sshash::util::ends_with(index_filename,
                                       constants::diff_colored_fulgor_filename_extension)) {
        return intersect_color_sets<differential_index_type>(index_filename, size, algo, quiet);
    } else if (sshash::util::ends_with(index_filename, constants::fulgor_filename_extension)) {
        return intersect_color_sets<index_type>(index_filename, size, algo, quiet);
    }

    std::cerr << "Wrong filename supplied." << std::endl;

    return 1;
}

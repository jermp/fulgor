#include <iostream>

#include "../include/index_types.hpp"
#include "../src/index.cpp"
#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"

#include "build.cpp"
#include "pseudoalign.cpp"
#include "sketch.cpp"

using namespace fulgor;

template <typename FulgorIndex>
void print_stats(std::string const& index_filename) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    index.print_stats();
}

template <typename FulgorIndex>
void print_filenames(std::string const& index_filename) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    for (uint64_t i = 0; i != index.num_docs(); ++i) {
        std::cout << i << '\t' << index.filename(i) << '\n';
    }
}

int stats(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("meta", "Specify if the Fulgor index is meta-colored.", "--meta", false, true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    auto meta = parser.get<bool>("meta");
    if (meta) {
        print_stats<meta_index_type>(index_filename);
    } else {
        print_stats<index_type>(index_filename);
    }
    return 0;
}

int print_filenames(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("meta", "Specify if the Fulgor index is meta-colored.", "--meta", false, true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    auto meta = parser.get<bool>("meta");
    essentials::logger("loading index from disk...");
    if (meta) {
        print_filenames<meta_index_type>(index_filename);
    } else {
        print_filenames<index_type>(index_filename);
    }
    return 0;
}

int permute_filenames(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("cluster_labels", "Cluster labels file.", "-c", true);
    parser.add("num_clusters", "Number of clusters (should match those in the cluster labels file)",
               "-n", true);
    parser.add("output_filename", "Cluster size filename", "-o", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    index_type index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    uint64_t num_clusters = parser.get<uint64_t>("num_clusters");

    std::vector<std::vector<std::string_view>> clusters(num_clusters,
                                                        std::vector<std::string_view>());
    auto cluster_labels = parser.get<std::string>("cluster_labels");
    std::ifstream in(cluster_labels);  // we assume there are exactly [num_docs] labels
    for (uint64_t i = 0; i != index.num_docs(); ++i) {
        uint32_t label = 0;
        in >> label;
        assert(label < num_clusters);
        clusters[label].push_back(index.filename(i));
    }
    in.close();

    /* sort references lexicographically within each cluster */
    for (auto& c : clusters) {
        std::sort(c.begin(), c.end(), [](auto const& fn1, auto const& fn2) { return fn1 < fn2; });
    }

    /* sort by non-increasing cluster size */
    std::sort(clusters.begin(), clusters.end(),
              [](auto const& c1, auto const& c2) { return c1.size() > c2.size(); });

    auto output_filename = parser.get<std::string>("output_filename");
    std::ofstream out(output_filename);
    uint64_t begin = 0;
    for (auto const& cluster : clusters) {
        for (auto const& fn : cluster) std::cerr << fn << '\n';
        uint64_t end = begin + cluster.size();  // one-past the end
        out << begin << " ";
        begin = end;
    }
    out << begin << "\n";
    out.close();

    return 0;
}

int partition(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename to partition.", "-i", true);
    parser.add("partitions_filename", "Partition file.", "-p", true);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("check", "Check correctness after index construction (it might take some time).",
               "--check", false, true);

    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    build_configuration build_config;
    build_config.index_filename_to_partition = parser.get<std::string>("index_filename");
    build_config.partitions_filename = parser.get<std::string>("partitions_filename");

    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    build_config.check = parser.get<bool>("check");

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();

    meta_index_type index;
    typename meta_index_type::meta_builder builder(build_config);
    builder.build(index);
    index.print_stats();

    timer.stop();
    essentials::logger("DONE");
    std::cout << "** building the index took " << timer.elapsed() << " seconds / "
              << timer.elapsed() / 60 << " minutes" << std::endl;

    std::string output_filename = build_config.index_filename_to_partition + ".partitioned_index";
    essentials::logger("saving index to disk...");
    essentials::save(index, output_filename.c_str());
    essentials::logger("DONE");

    return 0;
}

int help(char* arg0) {
    std::cout << "== Fulgor: a colored compacted de Bruijn graph index =========================="
              << std::endl
              << std::endl;
    std::cout << "Usage: " << arg0 << " <tool> ...\n\n"
              << "Available tools:\n"
              << "  build             \t build a Fulgor index\n"
              << "  pseudoalign       \t pseudoalign reads to references using a Fulgor index\n"
              << "  stats             \t print index statistics\n"
              << "  print-filenames   \t print all reference filenames\n"
              << "  sketch            \t build reference sketches\n"
              << "  permute-filenames \t permute filenames according to clusters\n"
              << "  partition         \t partition a single Fulgor index" << std::endl;
    return 1;
}

int main(int argc, char** argv) {
    if (argc < 2) return help(argv[0]);
    auto tool = std::string(argv[1]);
    if (tool == "build") {
        return build(argc - 1, argv + 1);
    } else if (tool == "pseudoalign") {
        return pseudoalign(argc - 1, argv + 1);
    } else if (tool == "stats") {
        return stats(argc - 1, argv + 1);
    } else if (tool == "print-filenames") {
        return print_filenames(argc - 1, argv + 1);
    } else if (tool == "sketch") {
        return sketch_references(argc - 1, argv + 1);
    } else if (tool == "permute-filenames") {
        return permute_filenames(argc - 1, argv + 1);
    } else if (tool == "partition") {
        return partition(argc - 1, argv + 1);
    }
    std::cout << "Unsupported tool '" << tool << "'.\n" << std::endl;
    return help(argv[0]);
}

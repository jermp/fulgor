using namespace fulgor;

int print_lists_in_cluster(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("labels_filename", "Cluster labels.", "-l", true);
    parser.add("target_label", "Cluster label whose lists we want to print.", "-c", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    auto index_filename = parser.get<std::string>("index_filename");
    index_type index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");

    uint64_t num_color_classes = index.num_color_classes();
    auto const& ccs = index.color_classes();

    std::ifstream in(parser.get<std::string>("labels_filename"));
    if (!in.is_open()) throw std::runtime_error("error in opening file");
    uint64_t label;
    uint64_t target_label = parser.get<uint64_t>("target_label");
    for (uint64_t color_id = 0; color_id != num_color_classes; ++color_id) {
        in >> label;
        if (label == target_label) {  // print list
            auto it = ccs.colors(color_id);
            uint64_t size = it.size();
            std::cout << "color_id " << color_id << ":\n";
            for (uint64_t i = 0; i != size; ++i, ++it) {
                uint32_t ref_id = *it;
                std::cout << ref_id << ' ';
                assert(ref_id < index.num_docs());
            }
            std::cout << std::endl;
        }
    }

    return 0;
}
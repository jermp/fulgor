using namespace fulgor;

bool is_meta(std::string index_filename) {
    return sshash::util::ends_with(index_filename,
                                   constants::meta_colored_fulgor_filename_extension);
}

bool is_meta_diff(std::string index_filename) {
    return sshash::util::ends_with(index_filename,
                                   constants::meta_diff_colored_fulgor_filename_extension);
}

bool is_diff(std::string index_filename) {
    return sshash::util::ends_with(index_filename,
                                   constants::diff_colored_fulgor_filename_extension);
}

bool is_hybrid(std::string index_filename) {
    return sshash::util::ends_with(index_filename, constants::fulgor_filename_extension);
}

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
    for (uint64_t i = 0; i != index.num_colors(); ++i) {
        std::cout << i << '\t' << index.filename(i) << '\n';
    }
}

template <typename FulgorIndex>
void dump(std::string const& index_filename, std::string const& basename) {
    FulgorIndex index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");
    index.dump(basename);
}

int stats(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    if (is_meta(index_filename)) {
        print_stats<meta_index_type>(index_filename);
    } else if (is_meta_diff(index_filename)) {
        print_stats<meta_differential_index_type>(index_filename);
    } else if (is_diff(index_filename)) {
        print_stats<differential_index_type>(index_filename);
    } else if (is_hybrid(index_filename)) {
        print_stats<index_type>(index_filename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int print_filenames(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    if (is_meta_diff(index_filename)) {
        print_filenames<meta_differential_index_type>(index_filename);
    } else if (is_meta(index_filename)) {
        print_filenames<meta_index_type>(index_filename);
    } else if (is_diff(index_filename)) {
        print_filenames<differential_index_type>(index_filename);
    } else if (is_hybrid(index_filename)) {
        print_filenames<index_type>(index_filename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int dump(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    if (is_meta_diff(index_filename)) {
        std::string basename{index_filename.data(),
                             index_filename.length() -
                                 constants::meta_diff_colored_fulgor_filename_extension.length() -
                                 1};
        dump<meta_differential_index_type>(index_filename, basename);
    } else if (is_meta(index_filename)) {
        std::string basename{index_filename.data(),
                             index_filename.length() -
                                 constants::meta_colored_fulgor_filename_extension.length() - 1};
        dump<meta_index_type>(index_filename, basename);
    } else if (is_diff(index_filename)) {
        std::string basename{index_filename.data(),
                             index_filename.length() -
                                 constants::diff_colored_fulgor_filename_extension.length() - 1};
        dump<differential_index_type>(index_filename, basename);
    } else if (is_hybrid(index_filename)) {
        std::string basename{
            index_filename.data(),
            index_filename.length() - constants::fulgor_filename_extension.length() - 1};
        dump<index_type>(index_filename, basename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

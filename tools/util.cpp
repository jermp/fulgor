using namespace fulgor;

bool is_meta(std::string const& index_filename) {
    return sshash::util::ends_with(index_filename,
                                   constants::mfur_filename_extension);
}

bool is_meta_diff(std::string const& index_filename) {
    return sshash::util::ends_with(index_filename,
                                   constants::mdfur_filename_extension);
}

bool is_diff(std::string const& index_filename) {
    return sshash::util::ends_with(index_filename,
                                   constants::dfur_filename_extension);
}

bool is_hybrid(std::string const& index_filename) {
    return sshash::util::ends_with(index_filename, constants::hfur_filename_extension);
}

template <typename FulgorIndex>
void verify(std::string const& index_filename) {
    FulgorIndex index;
    essentials::version_number vnum(constants::current_version_number::x,   //
                                    constants::current_version_number::y,   //
                                    constants::current_version_number::z);  //
    essentials::loader l(index_filename.c_str());
    l.visit(vnum);
    std::cout << "read version number = " << vnum.to_string() << std::endl;
    util::check_version_number(vnum);
    essentials::logger("OK: Fulgor index is compatible with current library version.");
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

int verify(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    if (is_meta(index_filename)) {
        verify<mfur_index_t>(index_filename);
    } else if (is_meta_diff(index_filename)) {
        verify<mdfur_index_t>(index_filename);
    } else if (is_diff(index_filename)) {
        verify<dfur_index_t>(index_filename);
    } else if (is_hybrid(index_filename)) {
        verify<hfur_index_t>(index_filename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int stats(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);
    auto index_filename = parser.get<std::string>("index_filename");
    if (is_meta(index_filename)) {
        print_stats<mfur_index_t>(index_filename);
    } else if (is_meta_diff(index_filename)) {
        print_stats<mdfur_index_t>(index_filename);
    } else if (is_diff(index_filename)) {
        print_stats<dfur_index_t>(index_filename);
    } else if (is_hybrid(index_filename)) {
        print_stats<hfur_index_t>(index_filename);
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
        print_filenames<mdfur_index_t>(index_filename);
    } else if (is_meta(index_filename)) {
        print_filenames<mfur_index_t>(index_filename);
    } else if (is_diff(index_filename)) {
        print_filenames<dfur_index_t>(index_filename);
    } else if (is_hybrid(index_filename)) {
        print_filenames<hfur_index_t>(index_filename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

int dump(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("output_basename",
               "Output basename. If the supplied basename is F, the output will consist in three "
               "files: F.unitigs.fa, F.color_sets.txt, and F.metadata.txt. (If this is not "
               "supplied, the basename of the index is used instead.)",
               "-o", false);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    auto index_filename = parser.get<std::string>("index_filename");
    std::string output_basename("");
    if (parser.parsed("output_basename")) {
        output_basename = parser.get<std::string>("output_basename");
        assert(output_basename.length() != 0);
    }

    if (is_meta_diff(index_filename)) {
        std::string basename{index_filename.data(),
                             index_filename.length() -
                                 constants::mdfur_filename_extension.length() -
                                 1};
        dump<mdfur_index_t>(
            index_filename, output_basename.length() == 0 ? basename : output_basename);
    } else if (is_meta(index_filename)) {
        std::string basename{index_filename.data(),
                             index_filename.length() -
                                 constants::mfur_filename_extension.length() - 1};
        dump<mfur_index_t>(index_filename,
                              output_basename.length() == 0 ? basename : output_basename);
    } else if (is_diff(index_filename)) {
        std::string basename{index_filename.data(),
                             index_filename.length() -
                                 constants::dfur_filename_extension.length() - 1};
        dump<dfur_index_t>(index_filename,
                                      output_basename.length() == 0 ? basename : output_basename);
    } else if (is_hybrid(index_filename)) {
        std::string basename{
            index_filename.data(),
            index_filename.length() - constants::hfur_filename_extension.length() - 1};
        dump<hfur_index_t>(index_filename,
                         output_basename.length() == 0 ? basename : output_basename);
    } else {
        std::cerr << "Wrong filename supplied." << std::endl;
        return 1;
    }
    return 0;
}

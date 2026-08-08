#include <iostream>
#include <filesystem>
#include <map>

#include "external/sshash/external/gz/zip_stream.hpp"
#include "external/sshash/external/gz/zip_stream.cpp"
#include "external/sshash/src/builder/build.cpp"
#include "external/sshash/src/dictionary.cpp"
#include "external/sshash/src/info.cpp"
#include "external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "external/FQFeeder/include/FastxParser.hpp"
#include "external/FQFeeder/src/FastxParser.cpp"

#include "include/index_types.hpp"
#include "src/index.cpp"
#include "src/color_sets.cpp"

#include "util.cpp"
#include "build.cpp"
#include "permute.cpp"
#include "pseudoalign.cpp"
#include "kmer_conservation.cpp"
#include "kmer_matches.cpp"

void help(char* arg0) {
    auto v_num = essentials::version_number(constants::current_version_number::major,
                                            constants::current_version_number::minor,
                                            constants::current_version_number::patch)
                     .to_string();
    std::cout << std::format(
        "== Fulgor: a colored de Bruijn graph index (v{}) "
        "=======================================\n\n",
        v_num);

    std::cout << "Usage: " << arg0 << " <tool> ...\n\n";
}

struct tool_set {
    using tool_function = int (*)(int, char**);

    struct tool {
        tool(const std::string& name, const tool_function function, const std::string& description)
            : name(name), function(function), description(description) {}

        const std::string name;
        const tool_function function;
        const std::string description;
    };

    void add(std::string const& name, const tool_function function,
             std::string const& description) {
        tools.emplace_back(name, function, description);
        sections.back().second++;
        longest_name = std::max(longest_name, name.size());
    }

    void add_section(std::string&& section) { sections.emplace_back(section, 0); }

    int run(std::string const& name, int argc, char** argv) const {
        if (name == "help") {
            help(argv[0]);
            print();
            return 0;
        }
        for (auto& tool : tools) {
            if (tool.name != name) continue;
            return tool.function(argc - 1, argv + 1);
        }
        std::cout << "Unsupported tool '" << name << "'." << std::endl;

        return 1;
    }

    void print() const {
        auto it = tools.begin();
        for (auto& [sec, num] : sections) {
            std::cout << sec << std::endl;
            for (uint64_t i = 0; i < num; ++i, ++it) {
                auto tool = *it;
                std::cout << std::format("  {}{}{}\n", tool.name,
                                         std::string(longest_name + 5 - tool.name.size(), ' '),
                                         tool.description);
            }
            std::cout << std::endl;
        }
    }

private:
    std::string::size_type longest_name = 0;
    std::vector<tool> tools;
    std::vector<std::pair<std::string, uint64_t>> sections;
};

int main(int argc, char** argv) {
    tool_set tools;
    tools.add_section("Construction");
    tools.add("build", build, "build an index");
    tools.add("color", color, "build a meta- or a diff- or a meta-diff- index");
    tools.add("permute", permute, "permute the reference names of an index");

    tools.add_section("Queries");
    tools.add("pseudoalign", pseudoalign, "perform pseudoalignment to an index");
    tools.add("kmer-conservation", kmer_conservation,
              "print color set info for each positive kmer in query");
    tools.add("kmer-matches", kmer_matches,
              "print positive kmers per query and number of kmer matches per color");

    tools.add_section("Debug");
    tools.add("check", probabilistic_check,
              "perform an in-depth check to verify that an index was built correctly");
    tools.add("verify", verify, "verify that index works correctly with current library version");
    tools.add("stats", stats, "print index statistics");
    tools.add("print-filenames", print_filenames, "print all reference filenames");
    tools.add("dump", dump, "write unitigs and color sets of an index in text format");
    tools.add("load", load, "build an index from dump output");

    tools.add_section("Other");
    tools.add("help", nullptr, "print this helper and exit gracefully");

    if (argc < 2) {
        tools.run("help", argc, argv);
        return 1;
    }
    const auto tool = std::string(argv[1]);

    return tools.run(tool, argc, argv);
}

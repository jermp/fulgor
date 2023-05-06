#include <iostream>
#include <numeric>
#include <vector>

#include "../include/inverted_index.hpp"
#include "../include/index.hpp"
#include "../src/index.cpp"
#include "../include/index_types.hpp"
#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../external/sshash/external/pthash/external/mm_file/include/mm_file/mm_file.hpp"

using namespace fulgor;

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("file_base_name", "Cuttlefish input file_base_name.", "-i", true);
    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    auto file_base_name = parser.get<std::string>("file_base_name");
    uint64_t num_docs = 0;
    uint64_t num_total_integers = 0;

    /* build the permutation map */
    std::unordered_map<seg_id_t, seg_id_t> permutation;
    {
        mm::file_source<seg_id_t> in(file_base_name + ".map", mm::advice::sequential);
        uint64_t num_sequences = in.bytes() / (sizeof(seg_id_t) + sizeof(uint64_t));
        seg_id_t const* data = in.data();
        assert(num_sequences < (uint64_t(1) << 32));
        for (uint64_t i = 0; i != num_sequences; ++i) {
            seg_id_t seg_id = *data;
            permutation[seg_id] = i;
            data += 1;  // skip seg_id
            /* skip dictionary_entry */
            if constexpr (sizeof(seg_id_t) == 4) {
                data += 2;
            } else {
                assert(sizeof(seg_id_t) == 8);
                data += 1;
            }
        }
    }

    {
        mm::file_source<uint64_t> mm_index_file(file_base_name + ".inv_idx",
                                                mm::advice::sequential);
        inverted_index::iterator it(mm_index_file.data(), mm_index_file.size());
        num_docs = it.num_docs();
        std::cout << "num_docs: " << num_docs << std::endl;
        std::cout << "num_ints: " << it.num_ints() << std::endl;

        index_type index;
        std::string index_filename =
            file_base_name + "." + index_type::color_classes_type::type() + ".index";
        essentials::logger("loading index from disk...");
        essentials::load(index, index_filename.c_str());
        essentials::logger("DONE");
        index.print_stats();

        auto const& ccs = index.color_classes();

        uint64_t num_lists = 0;
        while (it.has_next()) {
            auto list_exp = it.list();
            auto seg_id = list_exp.seg_id();
            assert(permutation.find(seg_id) != permutation.cend());
            uint64_t unitig_id = permutation[seg_id];
            uint64_t color_class_id = index.u2c(unitig_id);
            auto fwd_it = ccs.colors(color_class_id);
            if (fwd_it.size() != list_exp.size()) {
                std::cerr << "error: expected list of size " << list_exp.size() << " but got "
                          << fwd_it.size() << std::endl;
            }
            uint64_t size = fwd_it.size();
            auto list_exp_it = list_exp.begin();
            for (uint64_t i = 0; i != size; ++i, ++fwd_it, ++list_exp_it) {
                uint64_t got_color = *fwd_it;
                uint64_t exp_color = *list_exp_it;
                if (got_color != exp_color) {
                    std::cerr << "error: expected color " << exp_color << " but got " << got_color
                              << std::endl;
                }
                // else {
                //     std::cout << "OK: expected color " << exp_color << " got " << got_color
                //               << std::endl;
                // }
            }

            num_total_integers += size;
            it.advance_to_next_list();
            ++num_lists;

            if (num_lists % 10000 == 0) {
                std::cerr << "checked " << num_lists << " lists" << std::endl;
            }
        }

        std::cout << "checked " << num_lists << " lists" << std::endl;

        mm_index_file.close();
        assert(num_total_integers == it.num_ints());
    }

    std::cout << "num integers: " << num_total_integers << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;

    return 0;
}
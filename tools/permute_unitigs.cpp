#include <unordered_map>
#include <chrono>

#include "../include/util.hpp"
#include "../external/sshash/external/pthash/external/essentials/include/essentials.hpp"
#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../external/sshash/external/pthash/external/mm_file/include/mm_file/mm_file.hpp"

using namespace fulgor;

struct lines_iterator {
    lines_iterator(char const* begin, char const* end) : m_begin(begin), m_end(end) {
        advance_to_next();
    }

    void advance_to_next() {
        read_line('\n');
        m_seg_id = seg_id_t(-1);
        if (has_next()) m_seg_id = std::stoull(m_sequence);
        read_line('\n');
    }

    seg_id_t seg_id() const { return m_seg_id; }
    std::string const& sequence() const { return m_sequence; }
    bool has_next() const { return !m_sequence.empty(); }

private:
    char const* m_begin;
    char const* m_end;
    seg_id_t m_seg_id;
    std::string m_sequence;

    void read_line(char delim) {
        char const* begin = m_begin;
        while (m_begin != m_end and *m_begin++ != delim)
            ;
        if (begin == m_begin) {
            m_sequence.assign("");
        } else {
            m_sequence.assign(reinterpret_cast<const char*>(begin), m_begin - begin - 1);
        }
    }
};

void permute_and_write(build_configuration const& build_config) {
    /* build the permutation map */
    std::unordered_map<seg_id_t, seg_id_t> permutation;
    {
        mm::file_source<seg_id_t> in(build_config.file_base_name + ".map", mm::advice::sequential);
        uint64_t num_sequences = in.bytes() / (sizeof(seg_id_t) + sizeof(uint64_t));
        seg_id_t const* data = in.data();
        assert(num_sequences < (uint64_t(1) << 32));
        for (uint64_t i = 0; i != num_sequences; ++i) {
            seg_id_t seg_id = *data;
            permutation[seg_id] = i;
            data += 1;  // skip seg_id
            // std::cerr << seg_id << " --> " << i << "; color: ";
            /* skip dictionary_entry */
            if constexpr (sizeof(seg_id_t) == 4) {
                // std::cerr << *reinterpret_cast<uint64_t const*>(data) << '\n';
                data += 2;
            } else {
                assert(sizeof(seg_id_t) == 8);
                // std::cerr << *data << '\n';
                data += 1;
            }
        }
    }

    const uint64_t limit = build_config.ram_limit_in_GiB * essentials::GiB;
    std::vector<std::pair<seg_id_t, std::string>> buffer;  // (header, DNA sequence)

    std::string run_identifier =
        std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    uint64_t num_sequences = permutation.size();
    uint64_t num_bases = 0;
    uint64_t bytes = 0;
    uint64_t num_files_to_merge = 0;

    auto get_tmp_output_filename = [&](uint64_t id) {
        return build_config.tmp_dirname + "/perm.tmp.run" + run_identifier + "." +
               std::to_string(id);
    };

    auto sort_and_flush = [&]() {
        if (buffer.empty()) return;

        std::cout << "sorting buffer..." << std::endl;
        std::sort(buffer.begin(), buffer.end(), [&](auto const& p_x, auto const& p_y) {
            /* must be found */
            assert(permutation.find(p_x.first) != permutation.cend());
            assert(permutation.find(p_y.first) != permutation.cend());
            return permutation[p_x.first] < permutation[p_y.first];
        });

        auto tmp_output_filename = get_tmp_output_filename(num_files_to_merge);
        std::cout << "saving to file '" << tmp_output_filename << "'..." << std::endl;
        std::ofstream out(tmp_output_filename.c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");
        for (auto const& p : buffer) out << p.first << '\n' << p.second << '\n';
        out.close();

        buffer.clear();
        bytes = 0;
        num_files_to_merge += 1;
    };

    { /* create sorted runs */
        std::ifstream is(build_config.file_base_name + ".fa");
        if (!is.is_open()) throw std::runtime_error("cannot open '.fa' file");
        std::string line;
        for (uint64_t i = 0; i != num_sequences; ++i) {
            std::getline(is, line, '>');
            std::getline(is, line, '\n');
            seg_id_t seg_id = std::stoull(line);
            assert(permutation.find(seg_id) != permutation.cend());
            std::getline(is, line, '\n');
            uint64_t seq_bytes = line.length() + sizeof(seg_id_t) + 8;
            if (bytes + seq_bytes > limit) sort_and_flush();
            bytes += seq_bytes;
            buffer.emplace_back(seg_id, line);
            num_bases += line.length();
            if (i != 0 and i % 1000000 == 0) {
                std::cout << "read " << i << " sequences, " << num_bases << " bases" << std::endl;
            }
        }
        is.close();
        sort_and_flush();
    }

    std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases" << std::endl;

    { /* merge sorted runs */
        assert(num_files_to_merge > 0);
        std::cout << "files to merge = " << num_files_to_merge << std::endl;

        std::vector<lines_iterator> iterators;
        std::vector<uint32_t> idx_heap;
        iterators.reserve(num_files_to_merge);
        idx_heap.reserve(num_files_to_merge);
        std::vector<mm::file_source<char>> mm_files(num_files_to_merge);

        auto heap_idx_comparator = [&](uint32_t i, uint32_t j) {
            /* must be found */
            assert(permutation.find(iterators[i].seg_id()) != permutation.cend());
            assert(permutation.find(iterators[j].seg_id()) != permutation.cend());
            return permutation[iterators[i].seg_id()] > permutation[iterators[j].seg_id()];
        };

        auto advance_heap_head = [&]() {
            auto idx = idx_heap.front();
            iterators[idx].advance_to_next();
            if (iterators[idx].has_next()) {  // percolate down the head
                uint64_t pos = 0;
                uint64_t size = idx_heap.size();
                while (2 * pos + 1 < size) {
                    uint64_t i = 2 * pos + 1;
                    if (i + 1 < size and heap_idx_comparator(idx_heap[i], idx_heap[i + 1])) ++i;
                    if (heap_idx_comparator(idx_heap[i], idx_heap[pos])) break;
                    std::swap(idx_heap[pos], idx_heap[i]);
                    pos = i;
                }
            } else {
                std::pop_heap(idx_heap.begin(), idx_heap.end(), heap_idx_comparator);
                idx_heap.pop_back();
            }
        };

        /* create the input iterators and make the heap */
        for (uint64_t i = 0; i != num_files_to_merge; ++i) {
            auto tmp_output_filename = get_tmp_output_filename(i);
            mm_files[i].open(tmp_output_filename, mm::advice::sequential);
            iterators.emplace_back(mm_files[i].data(), mm_files[i].data() + mm_files[i].size());
            idx_heap.push_back(i);
        }
        std::make_heap(idx_heap.begin(), idx_heap.end(), heap_idx_comparator);

        std::ofstream out((build_config.file_base_name + ".permuted.fa").c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open output file");

        uint64_t num_written_sequences = 0;
        while (!idx_heap.empty()) {
            auto const& it = iterators[idx_heap.front()];
            /* write permuted seqments in FASTA format for SSHash */
            out << '>' << it.seg_id() << '\n' << it.sequence() << '\n';
            num_written_sequences += 1;
            if (num_written_sequences % 1000000 == 0) {
                std::cout << "written sequences = " << num_written_sequences << "/" << num_sequences
                          << std::endl;
            }
            advance_heap_head();
        }
        std::cout << "written sequences = " << num_written_sequences << "/" << num_sequences
                  << std::endl;
        out.close();
        assert(num_written_sequences == num_sequences);

        /* remove tmp files */
        for (uint64_t i = 0; i != num_files_to_merge; ++i) {
            mm_files[i].close();
            auto tmp_output_filename = get_tmp_output_filename(i);
            std::remove(tmp_output_filename.c_str());
        }
    }
}

int permute_unitigs(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("file_base_name", "Cuttlefish input file_base_name.", "-i", true);
    parser.add("RAM",
               "RAM limit in GiB. Default value is " +
                   std::to_string(constants::default_ram_limit_in_GiB) + ".",
               "-g", false);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("verbose", "Verbose output during construction.", "--verbose", false, true);

    if (!parser.parse()) return 1;
    util::print_cmd(argc, argv);

    build_configuration build_config;
    build_config.file_base_name = parser.get<std::string>("file_base_name");
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    if (parser.parsed("RAM")) build_config.ram_limit_in_GiB = parser.get<double>("RAM");
    if (parser.parsed("verbose")) {
        assert(parser.get<bool>("verbose"));
        build_config.verbose = true;
    }

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::seconds> timer;
    timer.start();
    permute_and_write(build_config);
    timer.stop();
    std::cout << "** permuting took " << timer.elapsed() << " seconds / " << timer.elapsed() / 60
              << " minutes" << std::endl;

    return 0;
}
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "../external/CLI11.hpp"
#include "../external/sshash/include/gz/zip_stream.hpp"
#include "../external/FQFeeder/include/FastxParser.hpp"
#include "../external/FQFeeder/include/blockingconcurrentqueue.h"

using namespace fulgor;

void process_lines(
  index_type& index,
  moodycamel::BlockingConcurrentQueue<std::vector<std::vector<uint32_t>>*>& q,
  std::atomic<bool>& done,
  std::atomic<uint32_t>& in_flight,
  std::mutex& ofile_mut,
  std::ofstream& out_file) {

  std::stringstream ss;
  int32_t buff_size{0};
  constexpr int32_t buff_thresh{100};


  std::vector<uint32_t> colors;

  std::vector<std::vector<uint32_t>>* sbatch;
  while (!done or in_flight > 0) {
    if (q.try_dequeue(sbatch)) {
      in_flight -= 1;
      for (auto& cids : *sbatch) {
        index.intersect_color_ids(cids, colors);
        if (!colors.empty()) {
          //num_mapped_reads += 1;
          ss << "read_name" << "\t" << colors.size();
          for (auto c : colors) { ss << "\t" << c; }
          ss << "\n";
        } else {
          ss << "read_name" << "\t0\n";
        }
        buff_size += 1;
        if (buff_size > buff_thresh) {
          std::string outs = ss.str();
          ss.str("");
          ofile_mut.lock();
          out_file.write(outs.data(), outs.size());
          ofile_mut.unlock();
          buff_size = 0;
        }

        colors.clear();
      }
      delete sbatch;
     }
  }
  // dump anything left in the buffer
  if (buff_size > 0) {
    std::string outs = ss.str();
    ss.str("");
    ofile_mut.lock();
    out_file.write(outs.data(), outs.size());
    ofile_mut.unlock();
    buff_size = 0;
  }

}

void do_intersection(index_type& index, const std::string& query_filename, std::ofstream& output) {

  std::ifstream ifile(query_filename, std::ios::binary);

  std::atomic<bool> done{false};
  std::atomic<uint32_t> in_flight{0};

  size_t num_threads = 15;
  moodycamel::BlockingConcurrentQueue< std::vector<std::vector<uint32_t>>* > q(3*num_threads);

  std::thread producer([&q, &ifile, &done, &in_flight]() -> void {
    uint64_t ctr = 0;
    
    std::vector<std::vector<uint32_t>>* batch = new std::vector<std::vector<uint32_t>>();
    uint32_t list_len = 0;

    while (ifile.read( reinterpret_cast<char*>(&list_len), sizeof(list_len) )) {
      if (list_len > 0) {
        if (list_len > 100) { std::cerr << "list_len " << list_len << "\n"; }
        std::vector<uint32_t> working_vec;
        working_vec.resize(list_len);
        ifile.read( reinterpret_cast<char*>(working_vec.data()), list_len * sizeof(list_len) );

        batch->push_back(working_vec);
        ctr++;
        if (ctr % 10000 == 0) {
          std::cerr << "processed " << ctr << " lines\n";
        }
      } else {
        std::cerr << "should not happen\n";
      }
      if (batch->size() >= 100) {
        q.enqueue(batch);
        batch = new std::vector<std::vector<uint32_t>>();
        in_flight += 1;
      }
    }
    if (!batch->empty()) {
      q.enqueue(batch);
      in_flight += 1;
    }
    done = true;
  });

  std::mutex ofmut;
  std::vector<std::thread> workers;
  workers.reserve(num_threads);
  for (size_t i = 0; i < num_threads; ++i) {
    workers.push_back(
      std::thread([&index, &q, &done, &in_flight, &ofmut, &output]() -> void {
       process_lines(index, q, done, in_flight, ofmut, output);
    }));
  }

  producer.join();
  for(auto& t : workers) { t.join(); }

}

int do_color_map(index_type const& index, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser,
           std::atomic<uint64_t>& num_reads, std::atomic<uint64_t>& num_mapped_reads,
           std::ofstream& out_file,
           std::mutex& iomut, std::mutex& ofile_mut) {
    std::vector<uint32_t> colors;  // result of pseudo-alignment
    std::stringstream ss;
    int32_t buff_size{0};
    constexpr int32_t buff_thresh{100};

    auto rg = rparser.getReadGroup();
    while (rparser.refill(rg)) {
      for (auto const& record : rg) {

        index.pseudoalign_full_intersection_color_ids(record.seq, colors);

        buff_size += 1;
        uint32_t read_id = num_reads.fetch_add(1);
        if (!colors.empty()) {
          num_mapped_reads += 1;
          ss.write( reinterpret_cast<char*>(&read_id), sizeof(read_id) );

          uint32_t num_colors = static_cast<uint32_t>(colors.size());
          ss.write( reinterpret_cast<char*>(&num_colors), sizeof(num_colors) );
          
          for (auto c : colors) { 
            ss.write( reinterpret_cast<char*>(&c), sizeof(c) );
          }
        } else {
          ss.write( reinterpret_cast<char*>(&read_id), sizeof(read_id) );
          uint32_t num_colors = 0;
          ss.write( reinterpret_cast<char*>(&num_colors), sizeof(num_colors) );
        }
        colors.clear();
        if (num_reads > 0 and num_reads % 1000000 == 0) {
          iomut.lock();
          std::cout << "mapped " << num_reads << " reads" << std::endl;
          iomut.unlock();
        }
        if (buff_size > buff_thresh) {
          std::string outs = ss.str();
          ss.str("");
          ofile_mut.lock();
          out_file.write(outs.data(), outs.size());
          ofile_mut.unlock();
          buff_size = 0;
        }
      }
    }

    // dump anything left in the buffer
    if (buff_size > 0) {
        std::string outs = ss.str();
        ss.str("");
        ofile_mut.lock();
        out_file.write(outs.data(), outs.size());
        ofile_mut.unlock();
        buff_size = 0;
    }

    return 0;
}

void sort_file(const std::string& tmp_outname) {
  std::vector<std::vector<uint32_t>> v;
  std::ifstream ifile(tmp_outname, std::ios::binary);

  std::vector<uint32_t> vals;
  uint32_t read_num = 0;
  uint32_t c = 0;
  while (ifile.read( reinterpret_cast<char*>(&read_num), sizeof(read_num) )) {
    uint32_t num_colors = 0;
    ifile.read( reinterpret_cast<char*>(&num_colors), sizeof(num_colors) );
    vals.resize(num_colors);
    ifile.read( reinterpret_cast<char*>(vals.data()), num_colors * sizeof(num_colors) );
    if (vals.size() > 0) {
      v.push_back(vals);
      vals.clear();
    }
  }

  std::cerr << "begin sort\n";

  std::sort(v.begin(), v.end(), 
            [](const std::vector<uint32_t>& a,
               const std::vector<uint32_t>& b) -> bool {
  
      bool a_shorter = a.size() < b.size();
      auto ai = a.begin();
      auto bi = b.begin();
      while (ai != a.end() && bi != b.end()) {
          if (*ai < *bi) { return true; }
          else if (*ai > *bi) { return false; }
          ai++; bi++;
      }
      return a_shorter;
  });

  auto end = std::unique(v.begin(), v.end());

  std::cerr << "sorted!\n";
  std::cerr << "v[0][0] = " << v[0][0] << "\n";
  ifile.close();
  std::ofstream ofile(tmp_outname, std::ios::trunc | std::ios::binary);
  for (auto vec_it = v.begin(); vec_it != end; ++vec_it) {
    uint32_t s = vec_it->size();
    if (s > 100) { std::cerr << "on output, size of v is " << s << "\n"; }
    ofile.write( reinterpret_cast<char*>(&s), sizeof(s) );
    ofile.write( reinterpret_cast<char*>(vec_it->data()), sizeof(s) * s );
  }
  ofile.close();
}

int intersect_colors(int argc, char** argv) {
    std::string index_filename;
    std::string query_filename;
    std::string output_filename;
    size_t num_threads{1};

    CLI::App app{"Perform color intersections of a set of pre-computed color ID lists"};
    app.add_option("-i,--index", index_filename, "Fulgor index filename")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-q,--query", query_filename,
                   "Query filename in CLIST format")
        ->required();
    app.add_option("-o,--output", output_filename, "File where output should be written")
        ->required();
    app.add_option("-t,--threads", num_threads, "Number of threads")->default_val(1);
    CLI11_PARSE(app, argc, argv);

    index_type index;
    essentials::logger("loading index from disk...");
    essentials::load(index, index_filename.c_str());
    essentials::logger("DONE");

    // first map to get the color lists
    std::ifstream is(query_filename.c_str());
    if (!is.good()) {
        std::cerr << "error in opening the file '" + query_filename + "'" << std::endl;
        return 1;
    }

    essentials::logger("performing queries from file '" + query_filename + "'...");
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> t;
    t.start();

    std::atomic<uint64_t> num_mapped_reads{0};
    std::atomic<uint64_t> num_reads{0};

    auto query_filenames = std::vector<std::string>({query_filename});
    if (num_threads == 1) {
        num_threads += 1;
        essentials::logger(
            "1 thread was specified, but an additional thread will be allocated for parsing");
    }
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(query_filenames, num_threads,
                                                             num_threads - 1);

    rparser.start();
    std::vector<std::thread> workers;
    std::mutex iomut;
    std::mutex ofile_mut;

    char tmp_filename[] = "tmp_color_listXXXXXX";
    mkstemp(tmp_filename);

    std::string tmp_outname(tmp_filename);
    essentials::logger("writing temporary color lists to " + tmp_outname);
    std::ofstream tmp_file;
    tmp_file.open(tmp_outname, std::ios::out | std::ios::trunc);
    if (!tmp_file) {
        essentials::logger("could not open output file " + tmp_outname);
        return 1;
    }

    for (size_t i = 1; i < num_threads; ++i) {
        workers.push_back(std::thread([&index, &rparser, &num_reads, &num_mapped_reads, 
                                       &tmp_file, &iomut, &ofile_mut]() {
            do_color_map(index, rparser, num_reads, num_mapped_reads, tmp_file, iomut,
                         ofile_mut);
        }));
    }

    for (auto& w : workers) { w.join(); }
    rparser.stop();

    t.stop();
    essentials::logger("DONE");
    tmp_file.close();

    std::cout << "mapped " << num_reads << " reads" << std::endl;
    std::cout << "elapsed = " << t.elapsed() << " millisec / ";
    std::cout << t.elapsed() / 1000 << " sec / ";
    std::cout << t.elapsed() / 1000 / 60 << " min / ";
    std::cout << (t.elapsed() * 1000) / num_reads << " musec/read" << std::endl;
    std::cout << "num_mapped_reads " << num_mapped_reads << "/" << num_reads << " ("
              << (num_mapped_reads * 100.0) / num_reads << "%)" << std::endl;




    sort_file(tmp_outname);
 
    
  std::ofstream output(output_filename);
  do_intersection(index, tmp_outname, output);

  return 0;
}
 

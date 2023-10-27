#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>

#include "../external/CLI11.hpp"
#include "../external/sshash/include/gz/zip_stream.hpp"
#include "../external/FQFeeder/include/FastxParser.hpp"
#include "../external/FQFeeder/include/blockingconcurrentqueue.h"

#include "../external/unordered_dense/include/ankerl/unordered_dense.h"

using namespace fulgor;

struct color_info {
  std::vector<uint32_t>::iterator begin;
  std::vector<uint32_t>::iterator end;
  std::vector<uint32_t>::iterator curr;
  size_t size() const { return std::distance(begin, end); }
  bool is_exhausted() const { return curr >= end; }
};

struct custom_vec_hash {
    using is_avalanching = void;

    [[nodiscard]] auto operator()(std::vector<uint32_t> const& f) const noexcept -> uint64_t {
        return ankerl::unordered_dense::detail::wyhash::hash(f.data(), sizeof(uint32_t) * f.size());
    }
};

using frequent_map_t = ankerl::unordered_dense::map<std::vector<uint32_t>, std::vector<uint32_t>, custom_vec_hash>;

void next_geq(color_info& info, uint32_t tgt, uint32_t num_docs) {
  info.curr = std::lower_bound(info.curr, info.end, tgt);
}

void intersect_uncompressed(std::vector<color_info>& iterators, uint32_t num_docs, std::vector<uint32_t>& colors) {
    assert(colors.empty());

    if (iterators.empty()) return;

    std::sort(iterators.begin(), iterators.end(),
              [](auto const& x, auto const& y) { return x.size() < y.size(); });

    if (iterators[0].is_exhausted()) {
      //std::cerr << "shouldn't happen";
      return;
    }

    uint32_t candidate = *(iterators[0].curr);
    uint64_t i = 1;
    while (candidate < num_docs) {
        for (; i != iterators.size(); ++i) {
            next_geq(iterators[i], candidate, num_docs);
            uint32_t val = iterators[i].is_exhausted() ? num_docs : *(iterators[i].curr);
            if (val != candidate) {
                candidate = val;
                i = 0;
                break;
            }
        }
        if (i == iterators.size()) {
            colors.push_back(candidate);
            iterators[0].curr++;
            candidate = iterators[0].is_exhausted() ? num_docs : *(iterators[0].curr);
            i = 1;
        }
    }
}

void analyze_batch(std::vector<std::vector<uint32_t>>* batch, 
                   index_type& index,
                   frequent_map_t& frequent_map) {

  ankerl::unordered_dense::map<std::vector<uint32_t>, uint32_t, custom_vec_hash> count_map;
  std::vector<uint32_t> key;
  size_t k = 3;
  uint32_t max_elem = 0;
  for (auto& cids : *batch) {
    if (cids.size() < k) { continue; }
    for (size_t i = 0; i < cids.size() - k; i++) {
      key.assign(cids.begin() + i, cids.begin() + i + k);
      auto& val = count_map[key];
      val += 1;
      //max_elem = std::max(val, max_elem);
    }
  }

  std::vector<uint32_t> colors;
  uint32_t thresh = 2;//std::max(static_cast<uint32_t>(2), static_cast<uint32_t>(0.1 * max_elem) + 1);
  for (auto const& [key, val] : count_map) {
    auto kv = frequent_map.find(key);
    if (kv != frequent_map.end() and (val >= thresh)) {
      index.intersect_color_ids(key, 0, colors);
      frequent_map[key] = colors;
      colors.clear();
    }
  }
   
}

// extract colors that are in the freuqent_map from 
// the colors vector
std::vector<color_info> extract_frequent_colors(
  std::vector<uint32_t>& cids, 
  frequent_map_t& frequent_map) {

    size_t k = 3;
    if (cids.size() < k) { return {}; }

    std::vector<uint32_t> cids_out;
    cids_out.reserve(cids.size());

    std::vector<color_info> res;
    std::vector<uint32_t> key;
    for (size_t i = 0; i < cids.size(); i++) {
      bool not_found = true;
      if (i < cids.size() - k) {
        key.assign(cids.begin() + i, cids.begin() + i + k);
        auto kv = frequent_map.find(key);
        if (kv != frequent_map.end()) {
          not_found = false;
          res.push_back( {kv->second.begin(), kv->second.end(), kv->second.begin()} );
          i += k;
        }
      } 
      if (not_found) {
        cids_out.push_back(cids[i]);
      }
    }
    std::swap(cids, cids_out);
    return res;
}

struct pref_suf_bounds {
  uint32_t prefix_len{0};
  uint32_t suffix_len{0};
};

pref_suf_bounds find_max_prefix_suffix(std::vector<uint32_t>& prev_cid, std::vector<uint32_t>& new_cid) {
  if (prev_cid.empty() or new_cid.empty()) {
    return {0, 0};
  }
  uint32_t pctr{0};
  uint32_t sctr{0};
  {
    auto pit = prev_cid.begin();
    auto cit = new_cid.begin();
    // count the length of the LCP
    for (; (*pit == *cit) and (pit < prev_cid.end()) and (cit < new_cid.end()); ++pit, ++cit, ++pctr) { }
  }

  {
    auto pit = prev_cid.rbegin();
    auto cit = new_cid.rbegin();
    // count the length of the LCS
    for (; (*pit == *cit) and (pit < prev_cid.rend()) and (cit < new_cid.rend()); ++pit, ++cit, ++sctr) { }
  }

  return {pctr, sctr};
}

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

  frequent_map_t frequent_map;
  std::vector<uint32_t> pkey;
  std::vector<uint32_t> skey;
  std::vector<uint32_t> colors;
  std::vector<uint32_t> colors_tmp;

  std::vector<std::vector<uint32_t>>* sbatch;

  size_t total_cid_len = 0;
  size_t total_ps_len = 0;

  size_t batch_ctr = 0;
  while (!done or in_flight > 0) {
    if (q.try_dequeue(sbatch)) {
      in_flight -= 1;
      frequent_map.clear();
      /*
      if (batch_ctr >= 5) {
        frequent_map.clear();
        batch_ctr = 0;
      }
      */
      batch_ctr++;
      //analyze_batch(sbatch, index, frequent_map);

      std::vector<uint32_t> prev_cids{};
      for (auto& cids : *sbatch) {
        pkey.clear();
        skey.clear();
        
        auto ps_info = find_max_prefix_suffix(prev_cids, cids);
        if (ps_info.prefix_len > 0) {
          pkey.assign(cids.begin(), cids.begin() + ps_info.prefix_len);
        }
        if (ps_info.suffix_len > 0) {
          skey.assign(cids.rbegin(), cids.rbegin() + ps_info.suffix_len);
        }
 
        std::vector<color_info> uncompressed_res;
        auto pkey_bak = pkey;
        auto skey_bak = skey;

        auto match_key = [&](auto& key, auto& key_bak) -> auto {
          if (key.empty()) { return frequent_map.end(); }
          auto pref_it = frequent_map.find(key);
          while (!key.empty()) {
            if (pref_it == frequent_map.end()) {
              key.pop_back();
              pref_it = frequent_map.find(key);
            } else {
              uncompressed_res.push_back(
                {pref_it->second.begin(), pref_it->second.end(), pref_it->second.begin()}
              );
              break;
            }
          }

          // if we didn't match the whole prefix, then 
          // compute the result for this prefix and 
          // put it in our hash
          if ( (key.size() >= 4) or (key.size() == key_bak.size()) ) {
            std::swap(key, key_bak);
          } else {
            colors_tmp.clear();
            std::vector<uint32_t> rem(key_bak.begin() + key.size(), key_bak.end());
            index.intersect_color_ids(rem, 0, colors_tmp);
            uncompressed_res.push_back(
              {colors_tmp.begin(), colors_tmp.end(), colors_tmp.begin()}
            );
            colors.clear();
            intersect_uncompressed(uncompressed_res, index.num_docs(), colors);
            frequent_map[key_bak] = colors;
            colors.clear();
          }
          uncompressed_res.clear();
          return frequent_map.find(key_bak);
        };

        auto pit = match_key(pkey, pkey_bak);
        auto sit = match_key(skey, skey_bak);

        if (pit != frequent_map.end()) {
          uncompressed_res.push_back({pit->second.begin(), pit->second.end(), pit->second.begin()});
        }
        if (sit != frequent_map.end()) {
          uncompressed_res.push_back({sit->second.begin(), sit->second.end(), sit->second.begin()});
        }

        total_cid_len += cids.size();
        total_ps_len += pkey_bak.size() + skey_bak.size();
        prev_cids = cids;
        std::vector<uint32_t>(cids.begin() + pkey_bak.size(), cids.end() - skey_bak.size()).swap(cids);
        if (!cids.empty()) {
          colors_tmp.clear();
          index.intersect_color_ids(cids, 0, colors_tmp);
          uncompressed_res.push_back({colors_tmp.begin(), colors_tmp.end(), colors_tmp.begin()});
        }
        colors.clear();
        intersect_uncompressed(uncompressed_res, index.num_docs(), colors);
        /*
        // extract colors that are in the freuqent_map from 
        // the colors vector
        std::vector<color_info> frequent_res = extract_frequent_colors(cids, frequent_map);

        // intersect what remains
        colors_tmp.clear();
        if (!cids.empty()) { 
          index.intersect_color_ids(cids, colors_tmp);
        }
        // add the color vector we just computed
        if (!colors_tmp.empty()) {
          frequent_res.push_back( {colors_tmp.begin(), colors_tmp.end(), colors_tmp.begin()} );
        }
        // intersect_uncompressed
        intersect_uncompressed(frequent_res, index.num_docs(), colors);
        */

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

  std::cerr << "total cid len: " << total_cid_len << ", total_ps_len: " << total_ps_len << " ratio " << static_cast<double>(total_ps_len)/total_cid_len << "\n";

}


void process_lines_simple(
  index_type& index,
  moodycamel::BlockingConcurrentQueue<std::vector<std::vector<uint32_t>>*>& q,
  std::atomic<bool>& done,
  std::atomic<uint32_t>& in_flight,
  std::mutex& ofile_mut,
  std::ofstream& out_file,
  size_t thread_index) {

  std::stringstream ss;
  int32_t buff_size{0};
  constexpr int32_t buff_thresh{500};

  std::vector<uint32_t> colors;
  std::vector<std::vector<uint32_t>>* sbatch = nullptr;

  size_t skipped_aln = 0;
  size_t total_aln = 0;
  size_t batch_ctr = 0;
  while (q.try_dequeue(sbatch) or !done or (in_flight > 0)) {
    if (sbatch != nullptr) {
      in_flight -= 1;
      batch_ctr++;

      // initially for each batch, the output
      // should be empty.
      colors.clear();

      for (auto& cids : *sbatch) {
        ++total_aln;
        // if the only thing in the cids vector is 
        // 1 element (the read_id), then the color 
        // output should be exactly the same as the 
        // last call to index.instersect_color_ids
        // and so we don't recompute it here.
        if (cids.size() > 1) {
          // interesect the color ids to get the colors
          colors.clear();
          index.intersect_color_ids(cids, 1, colors);
        } else {
          ++skipped_aln;
        }

        if (!colors.empty()) {
          ss << cids.front() << "\t" << colors.size();
          for (auto c : colors) { ss << "\t" << c; }
          ss << "\n";
        } else {
          ss << cids.front() << "\t0\n";
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

      }
      delete sbatch;
      sbatch = nullptr;
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
  std::cerr << "(thread_index : " << thread_index << ") total_aln: " << total_aln << ", skipped_aln: " << skipped_aln << ", batch_ctr: " << batch_ctr << "\n";
}

void do_intersection(index_type& index, size_t num_threads_in, const std::string& query_filename, std::ofstream& output) {

  std::ifstream ifile(query_filename, std::ios::binary);

  std::atomic<bool> done{false};
  std::atomic<uint32_t> in_flight{0};

  moodycamel::BlockingConcurrentQueue< std::vector<std::vector<uint32_t>>* > q(3*num_threads_in);

  std::thread producer([&q, &ifile, &done, &in_flight]() -> void {
    uint64_t ctr = 0;
    
    std::vector<std::vector<uint32_t>>* batch = new std::vector<std::vector<uint32_t>>();
    uint32_t list_len = 0;
    size_t nbatches = 0;
    while (ifile.read( reinterpret_cast<char*>(&list_len), sizeof(list_len) )) {
      if (list_len > 0) {
        //if (list_len > 100) { std::cerr << "list_len " << list_len << "\n"; }
        std::vector<uint32_t> working_vec;
        working_vec.resize(list_len);
        ifile.read( reinterpret_cast<char*>(working_vec.data()), list_len * sizeof(list_len) );

        // if the current batch is big enough that we want to push it
        // make sure the current element isn't a duplicate (i.e. size 1)
        // and push the batch *before* we add the next element, which 
        // may be the start of a new duplicate run.
        if ((batch->size() >= 10000) and (working_vec.size() > 1)) {
          while (!q.try_enqueue(batch)) {}
          batch = new std::vector<std::vector<uint32_t>>();
          in_flight += 1;
          ++nbatches;
        }
        // always push back the current element. If it's in the middle of 
        // a duplicate run, we keep growing the previous batch, otherwise
        // this is the first element of the enw batch.
        batch->push_back(working_vec);
      } else {
        std::cerr << "should not happen\n";
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
  workers.reserve(num_threads_in);
  for (size_t i = 0; i < num_threads_in; ++i) {
    workers.push_back(
      std::thread([&index, &q, &done, &in_flight, &ofmut, &output, i]() -> void {
       process_lines_simple(index, q, done, in_flight, ofmut, output, i);
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

void sort_file(const std::string& tmp_outname, std::ofstream& output) {
  std::vector<std::vector<uint32_t>> v;
  std::ifstream ifile(tmp_outname, std::ios::binary);

  std::vector<uint32_t> vals;
  uint32_t read_num = 0;
  while (ifile.read( reinterpret_cast<char*>(&read_num), sizeof(read_num) )) {
    uint32_t num_colors = 0;
    ifile.read( reinterpret_cast<char*>(&num_colors), sizeof(num_colors) );
    
    vals.resize(num_colors + 1);
    vals[0] = read_num;

    ifile.read( reinterpret_cast<char*>(&vals[1]), num_colors * sizeof(num_colors) );
    if (vals.size() > 1) {
      v.push_back(vals);
      vals.clear();
    } else {
      // just write out the unmapped reads here
      output << read_num << "\t0\n";
    }
  }

  if (v.empty()) { return; }

  std::cerr << "begin sort\n";
  
  std::sort(v.begin(), v.end(), 
            [](const std::vector<uint32_t>& a,
               const std::vector<uint32_t>& b) -> bool {
              return std::lexicographical_compare(a.begin() + 1, a.end(), b.begin() + 1, b.end());
            });

  std::cerr << "done sort\n";

  auto curr = v.begin();
  auto next = curr++;
  size_t identical_lists = 0;
  while (next < v.end()) {
    if ((curr->size() == next->size()) and std::equal(curr->begin() + 1, curr->end(), next->begin() + 1)) {
      next->resize(1); // retain only the id.
      ++identical_lists;
      ++next;
    } else {
      curr = next;
      ++next;
    }
  }

  std::cerr << "number of identical lists = " << identical_lists << "\n";
  std::cerr << "sorted!\n";
  ifile.close();
  std::ofstream ofile(tmp_outname, std::ios::trunc | std::ios::binary);
  for (auto vec_it = v.begin(); vec_it != v.end(); ++vec_it) {
    uint32_t s = vec_it->size();
    //if (s > 100) { std::cerr << "on output, size of v is " << s << "\n"; }
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




  std::ofstream output(output_filename);
  sort_file(tmp_outname, output);
 
    
  do_intersection(index, num_threads, tmp_outname, output);
  std::remove(tmp_outname.c_str());

  return 0;
}
 

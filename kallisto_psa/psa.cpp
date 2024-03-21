#include "psa.hpp"

// The below function is meant to mimic, as closely as possible,
// the behavior of Kallisto. This function is taken directly from
// https://github.com/pachterlab/kallisto, from the KmerIndex.cpp file
// and is modified only to make use of the relevant index and data
// structure in this codebase.
//
// use:  match(s,l,idx,v)
// pre:  v is initialized
// post: v contains all equiv classes for the k-mers in s

template <typename FulgorIndex>
void match(const std::string& s, int l, FulgorIndex const* idx,
           std::vector<std::pair<projected_hits, int>>& v) {
    sshash::dictionary const& kmap = idx->get_k2u();
    int64_t k = kmap.k();
    // NOTE: have to set this first, before starting any iterator or anything that
    // uses k-mers!!
    CanonicalKmer::k(k);
    pufferfish::CanonicalKmerIterator kit(s), kit_end;

    int skip = 1;

    auto getDist = [k](projected_hits& phits) -> int {
        size_t cStartPos = phits.globalPos_ - phits.contigPos_;
        size_t cEndPos = cStartPos + phits.contigLen_;
        size_t cCurrPos = phits.globalPos_;

        int32_t ctg_skip = 1;
        // fw ori
        if (phits.contigOrientation_) {
            ctg_skip = static_cast<int64_t>(cEndPos) - (static_cast<int64_t>(cCurrPos + k));
        } else {  // rc ori
            ctg_skip = static_cast<int32_t>(phits.contigPos_);
        }
        return ctg_skip;
    };

    bool backOff = false;
    int nextPos = 0;  // nextPosition to check
    for (int i = 0; kit != kit_end; ++i, ++kit) {
        // need to check it
        auto search = kmap.lookup_advanced_uint(kit->first.fwWord());

        int pos = kit->second;

        if (search.kmer_id != sshash::constants::invalid_uint64) {
            // KmerEntry val = search->second;
            projected_hits val;
            val.globalPos_ = search.kmer_id + (search.contig_id * (k - 1));
            val.contigPos_ = search.kmer_id_in_contig;
            val.contigIdx_ = search.contig_id;
            val.contigLen_ = search.contig_size + k - 1;
            val.contigOrientation_ =
                (search.kmer_orientation == sshash::constants::forward_orientation);

            v.push_back({val, kit->second});

            // see if we can skip ahead
            // bring thisback later
            // bool forward = val.contigOrientation_;//(kit->first == search->first);
            int dist = getDist(val);

            // const int lastbp = 10;
            if (dist >= 2) {
                // where should we jump to?
                int nextPos = pos + dist;  // default jump

                if (pos + dist >= l - k) {
                    // if we can jump beyond the read, check the end
                    nextPos = l - k;
                }

                // check next position
                pufferfish::CanonicalKmerIterator kit2(kit);
                kit2.jumpTo(nextPos);
                if (kit2 != kit_end) {
                    auto search2 = kmap.lookup_advanced_uint(kit2->first.fwWord());
                    bool found2 = false;
                    int found2pos = pos + dist;
                    if (search2.kmer_id == sshash::constants::invalid_uint64) {
                        found2 = true;
                        found2pos = pos;
                    } else if (val.contigIdx_ == search2.contig_id) {
                        found2 = true;
                        found2pos = pos + dist;
                    }
                    if (found2) {
                        // great, a match (or nothing) see if we can move the k-mer forward
                        if (found2pos >= l - k) {
                            v.push_back({val, l - k});  // push back a fake position
                            break;                      //
                        } else {
                            v.push_back({val, found2pos});
                            kit = kit2;  // move iterator to this new position
                        }
                    } else {
                        // this is weird, let's try the middle k-mer
                        bool foundMiddle = false;
                        if (dist > 4) {
                            int middlePos = (pos + nextPos) / 2;
                            uint64_t middleContig = sshash::constants::invalid_uint64;
                            int found3pos = pos + dist;
                            pufferfish::CanonicalKmerIterator kit3(kit);
                            kit3.jumpTo(middlePos);
                            projected_hits val3;
                            if (kit3 != kit_end) {
                                auto search3 = kmap.lookup_advanced_uint(kit3->first.fwWord());
                                if (search3.kmer_id != sshash::constants::invalid_uint64) {
                                    val3.globalPos_ =
                                        search3.kmer_id + (search3.contig_id * (k - 1));
                                    val3.contigPos_ = search3.kmer_id_in_contig;
                                    val3.contigIdx_ = search3.contig_id;
                                    val3.contigLen_ = search3.contig_size + k - 1;
                                    val3.contigOrientation_ =
                                        (search3.kmer_orientation ==
                                         sshash::constants::forward_orientation);

                                    middleContig = search3.contig_id;
                                    if (middleContig == val.contigIdx_) {
                                        foundMiddle = true;
                                        found3pos = middlePos;
                                    } else if (middleContig == search2.contig_id) {
                                        foundMiddle = true;
                                        found3pos = pos + dist;
                                    }
                                }

                                if (foundMiddle) {
                                    v.push_back({val3, found3pos});
                                    if (nextPos >= l - k) {
                                        break;
                                    } else {
                                        kit = kit2;
                                    }
                                }
                            }
                        }

                        if (!foundMiddle) {
                            ++kit;
                            backOff = true;
                            goto donejumping;  // sue me Dijkstra!
                        }
                    }
                } else {
                    // the sequence is messed up at this point, let's just take the match
                    // v.push_back({dbGraph.ecs[val.contig], l-k});
                    break;
                }
            }
        }

    donejumping:

        if (backOff) {
            // backup plan, let's play it safe and search incrementally for the rest, until nextStop
            for (int j = 0; kit != kit_end; ++kit, ++j) {
                if (j == skip) { j = 0; }
                if (j == 0) {
                    // need to check it
                    auto search = kmap.lookup_advanced_uint(kit->first.fwWord());
                    if (search.kmer_id != sshash::constants::invalid_uint64) {
                        projected_hits tmpval;
                        tmpval.globalPos_ = search.kmer_id + (search.contig_id * (k - 1));
                        tmpval.contigPos_ = search.kmer_id_in_contig;
                        tmpval.contigIdx_ = search.contig_id;
                        tmpval.contigLen_ = search.contig_size + k - 1;
                        tmpval.contigOrientation_ =
                            (search.kmer_orientation == sshash::constants::forward_orientation);

                        // if k-mer found
                        v.push_back({tmpval, kit->second});  // add equivalence class, and position
                    }
                }

                if (kit->second >= nextPos) {
                    backOff = false;
                    break;  // break out of backoff for loop
                }
            }
        }
    }
}

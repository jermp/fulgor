#pragma once

#include "external/sshash/include/constants.hpp"
#include <iostream>

struct projected_hits {
    projected_hits()
        : contigIdx_(0)
        , contigPos_(0)
        , contigOrientation_(0)
        , contigLen_(0)
        , globalPos_(0)
        , k_(0) {}

    uint64_t contigIdx_;
    // The relative position of the k-mer inducing this hit on the
    // contig
    uint32_t contigPos_;
    // How the k-mer inducing this hit maps to the contig
    // true for fw, false for rc
    bool contigOrientation_;
    uint32_t contigLen_;
    uint64_t globalPos_;
    uint32_t k_;

    // sshash::util::contig_span refRange;
    inline bool empty() {
        return contigIdx_ == sshash::constants::invalid_uint64;
    }  // refRange.empty(); }

    inline uint32_t contig_id() const { return contigIdx_; }
    inline bool hit_fw_on_contig() const { return contigOrientation_; }

    /*
    inline ref_pos decode_hit(uint64_t v) {
        // true if the contig is fowrard on the reference
        bool contigFW = sshash::util::orientation(v);
        // we are forward with respect to the reference if :
        // (1) contigFW and contigOrientation_
        // (2) !contigFW and !contigOrientation_
        // we are reverse complement with respect to the reference if :
        // (3) configFW and !contigOrientation_
        // (4) !configFW and contigOrientation_

        // if we're in the forward orientation, then our position is
        // just the contig offset plus or relative position
        uint32_t rpos;  //{0};
        bool rfw;       //{false};
        if (contigFW and contigOrientation_) {
            // kmer   :          AGC
            // contig :      ACTTAGC
            // ref    :  GCA[ACTTAGC]CA
            rpos = sshash::util::pos(v) + contigPos_;
            rfw = true;
        } else if (contigFW and !contigOrientation_) {
            // kmer   :          GCT
            // contig :      ACTTAGC
            // ref    :  GCA[ACTTAGC]CA
            rpos = sshash::util::pos(v) + contigPos_;
            rfw = false;
        } else if (!contigFW and contigOrientation_) {
            // kmer   :          AGT
            // contig :      GCTAAGT
            // ref    :  GCA[ACTTAGC]CA
            rpos = sshash::util::pos(v) + contigLen_ - (contigPos_ + k_);
            rfw = false;
        } else {  // if (!contigFW and !contigOrientation_) {
            // kmer   :          ACT
            // contig :      GCTAAGT
            // ref    :  GCA[ACTTAGC]CA
            rpos = sshash::util::pos(v) + contigLen_ - (contigPos_ + k_);
            rfw = true;
        }

        return {rpos, rfw};
    }
    */

    // inline friend function :
    // https://stackoverflow.com/questions/381164/friend-and-inline-method-whats-the-point
    // this helps to avoid duplicate symbol error.
    inline friend std::ostream& operator<<(std::ostream& os, projected_hits& h) {
        os << "{ proj_hit : \n"
           << "\t{ contig_idx : " << h.contigIdx_ << ", "
           << "contig_pos : " << h.contigPos_ << ", "
           << "contig_ori : " << (h.contigOrientation_ ? "fw" : "rc") << ", "
           << "contig_len : " << h.contigLen_ << ", "
           << "global_pos : " << h.globalPos_ << "}\n}\n";
        //<< "ref_range_len : " << h.refRange.size() << "}\n}\n";
        return os;
    }
};

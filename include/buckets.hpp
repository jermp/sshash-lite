#pragma once

#include "util.hpp"
#include "bit_vector_iterator.hpp"
#include "ef_sequence.hpp"

namespace sshash {

struct buckets {
    uint64_t offset_to_contig_end(uint64_t offset) const {
        auto [pos, contig_end] = pieces.next_geq(offset);

        /* The following two facts hold. */
        assert(pieces.access(pos) == contig_end);
        assert(contig_end >= offset);

        if (contig_end == offset) {
            assert(pos + 1 < pieces.size());
            contig_end = pieces.access(pos + 1);
        }

        /* Now, the following fact holds. */
        assert(offset < contig_end);

        return contig_end;
    }

    std::pair<uint64_t, uint64_t> locate_bucket(uint64_t bucket_id) const {
        uint64_t begin = num_super_kmers_before_bucket.access(bucket_id) + bucket_id;
        uint64_t end = num_super_kmers_before_bucket.access(bucket_id + 1) + bucket_id + 1;
        assert(begin < end);
        return {begin, end};
    }

    bool lookup(uint64_t bucket_id, uint64_t target_kmer, uint64_t k, uint64_t m) const {
        auto [begin, end] = locate_bucket(bucket_id);
        return lookup(begin, end, target_kmer, k, m);
    }

    bool lookup(uint64_t begin, uint64_t end, uint64_t target_kmer, uint64_t k, uint64_t m) const {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            bool answer = lookup_in_super_kmer(super_kmer_id, target_kmer, k, m);
            if (answer) return true;
        }
        return false;
    }

    bool lookup_in_super_kmer(uint64_t super_kmer_id, uint64_t target_kmer, uint64_t k,
                              uint64_t m) const {
        uint64_t offset = offsets.access(super_kmer_id);
        uint64_t contig_end = offset_to_contig_end(offset);
        bit_vector_iterator bv_it(strings, 2 * offset);
        uint64_t window_size = std::min<uint64_t>(k - m + 1, contig_end - offset - k + 1);
        for (uint64_t w = 0; w != window_size; ++w) {
            uint64_t read_kmer = bv_it.read_and_advance_by_two(2 * k);
            if (read_kmer == target_kmer) return true;
        }
        return false;
    }

    bool lookup_canonical(uint64_t bucket_id, uint64_t target_kmer, uint64_t target_kmer_rc,
                          uint64_t k, uint64_t m) const {
        auto [begin, end] = locate_bucket(bucket_id);
        return lookup_canonical(begin, end, target_kmer, target_kmer_rc, k, m);
    }

    bool lookup_canonical(uint64_t begin, uint64_t end, uint64_t target_kmer,
                          uint64_t target_kmer_rc, uint64_t k, uint64_t m) const {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = offsets.access(super_kmer_id);
            uint64_t contig_end = offset_to_contig_end(offset);
            bit_vector_iterator bv_it(strings, 2 * offset);
            uint64_t window_size = std::min<uint64_t>(k - m + 1, contig_end - offset - k + 1);
            for (uint64_t w = 0; w != window_size; ++w) {
                uint64_t read_kmer = bv_it.read_and_advance_by_two(2 * k);
                if (read_kmer == target_kmer or read_kmer == target_kmer_rc) return true;
            }
        }
        return false;
    }

    uint64_t num_bits() const {
        return pieces.num_bits() + num_super_kmers_before_bucket.num_bits() +
               8 * (offsets.bytes() + strings.bytes());
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(pieces);
        visitor.visit(num_super_kmers_before_bucket);
        visitor.visit(offsets);
        visitor.visit(strings);
    }

    ef_sequence<true> pieces;
    ef_sequence<false> num_super_kmers_before_bucket;
    pthash::compact_vector offsets;
    pthash::bit_vector strings;
};

}  // namespace sshash
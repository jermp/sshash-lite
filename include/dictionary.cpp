#include "dictionary.hpp"

namespace sshash {

bool dictionary::lookup_uint64_regular_parsing(uint64_t uint64_kmer) const {
    uint64_t minimizer = util::compute_minimizer(uint64_kmer, m_k, m_m, m_seed);
    uint64_t bucket_id = m_minimizers.lookup(minimizer);

    if (m_skew_index.empty()) return m_buckets.lookup(bucket_id, uint64_kmer, m_k, m_m);

    auto [begin, end] = m_buckets.locate_bucket(bucket_id);
    uint64_t num_super_kmers_in_bucket = end - begin;
    uint64_t log2_bucket_size = util::ceil_log2_uint32(num_super_kmers_in_bucket);
    if (log2_bucket_size > m_skew_index.min_log2) {
        uint64_t pos = m_skew_index.lookup(uint64_kmer, log2_bucket_size);
        /* It must hold pos < num_super_kmers_in_bucket for the kmer to exist. */
        if (pos < num_super_kmers_in_bucket) {
            return m_buckets.lookup_in_super_kmer(begin + pos, uint64_kmer, m_k, m_m);
        }
        return false;
    }

    return m_buckets.lookup(begin, end, uint64_kmer, m_k, m_m);
}

bool dictionary::lookup_uint64_canonical_parsing(uint64_t uint64_kmer) const {
    uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, m_k);
    uint64_t minimizer = util::compute_minimizer(uint64_kmer, m_k, m_m, m_seed);
    uint64_t minimizer_rc = util::compute_minimizer(uint64_kmer_rc, m_k, m_m, m_seed);
    uint64_t bucket_id = m_minimizers.lookup(std::min<uint64_t>(minimizer, minimizer_rc));

    if (m_skew_index.empty()) {
        return m_buckets.lookup_canonical(bucket_id, uint64_kmer, uint64_kmer_rc, m_k, m_m);
    }

    auto [begin, end] = m_buckets.locate_bucket(bucket_id);
    uint64_t num_super_kmers_in_bucket = end - begin;
    uint64_t log2_bucket_size = util::ceil_log2_uint32(num_super_kmers_in_bucket);
    if (log2_bucket_size > m_skew_index.min_log2) {
        uint64_t pos = m_skew_index.lookup(uint64_kmer, log2_bucket_size);
        if (pos < num_super_kmers_in_bucket) {
            bool answer = m_buckets.lookup_in_super_kmer(begin + pos, uint64_kmer, m_k, m_m);
            if (answer) return true;
        }
        uint64_t pos_rc = m_skew_index.lookup(uint64_kmer_rc, log2_bucket_size);
        if (pos_rc < num_super_kmers_in_bucket) {
            bool answer = m_buckets.lookup_in_super_kmer(begin + pos_rc, uint64_kmer_rc, m_k, m_m);
            return answer;
        }
        return false;
    }

    return m_buckets.lookup_canonical(begin, end, uint64_kmer, uint64_kmer_rc, m_k, m_m);
}

bool dictionary::lookup_advanced_uint64(uint64_t uint64_kmer,
                                        bool check_reverse_complement_too) const {
    if (m_canonical_parsing) return lookup_uint64_canonical_parsing(uint64_kmer);
    bool answer = lookup_uint64_regular_parsing(uint64_kmer);
    if (check_reverse_complement_too and answer == false) {
        uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, m_k);
        answer = lookup_uint64_regular_parsing(uint64_kmer_rc);
    }
    return answer;
}

bool dictionary::is_member(char const* string_kmer, bool check_reverse_complement_too) const {
    uint64_t uint64_kmer = util::string_to_uint64_no_reverse(string_kmer, m_k);
    return lookup_advanced_uint64(uint64_kmer, check_reverse_complement_too);
}
bool dictionary::is_member_uint64(uint64_t uint64_kmer, bool check_reverse_complement_too) const {
    return lookup_advanced_uint64(uint64_kmer, check_reverse_complement_too);
}

uint64_t dictionary::num_bits() const {
    return 8 * (sizeof(m_size) + sizeof(m_seed) + sizeof(m_k) + sizeof(m_m) +
                sizeof(m_canonical_parsing)) +
           m_minimizers.num_bits() + m_buckets.num_bits() + m_skew_index.num_bits();
}

}  // namespace sshash
#pragma once

#include "util.hpp"
#include "minimizers.hpp"
#include "buckets.hpp"
#include "skew_index.hpp"

namespace sshash {

struct dictionary {
    dictionary() : m_size(0), m_seed(0), m_k(0), m_m(0), m_canonical_parsing(0) {}

    void build(std::string const& filename, build_configuration const& build_config);

    uint64_t size() const { return m_size; }
    uint64_t seed() const { return m_seed; }
    uint64_t k() const { return m_k; }
    uint64_t m() const { return m_m; }
    uint64_t num_contigs() const { return m_buckets.pieces.size() - 1; }
    bool canonicalized() const { return m_canonical_parsing; }

    /* Membership queries. */
    bool is_member(char const* string_kmer, bool check_reverse_complement_too = true) const;
    bool is_member_uint64(uint64_t uint64_kmer, bool check_reverse_complement_too = true) const;

    /* Streaming queries. */
    friend struct streaming_query_canonical_parsing;
    friend struct streaming_query_regular_parsing;
    streaming_query_report streaming_query_from_file(std::string const& filename,
                                                     bool multiline) const;

    uint64_t num_bits() const;
    void print_info() const;
    void print_space_breakdown() const;

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_size);
        visitor.visit(m_seed);
        visitor.visit(m_k);
        visitor.visit(m_m);
        visitor.visit(m_canonical_parsing);
        visitor.visit(m_minimizers);
        visitor.visit(m_buckets);
        visitor.visit(m_skew_index);
    }

private:
    uint64_t m_size;
    uint64_t m_seed;
    uint16_t m_k;
    uint16_t m_m;
    uint16_t m_canonical_parsing;
    minimizers m_minimizers;
    buckets m_buckets;
    skew_index m_skew_index;

    bool lookup_advanced_uint64(uint64_t uint64_kmer,
                                bool check_reverse_complement_too = true) const;
    bool lookup_uint64_regular_parsing(uint64_t uint64_kmer) const;
    bool lookup_uint64_canonical_parsing(uint64_t uint64_kmer) const;
};

}  // namespace sshash
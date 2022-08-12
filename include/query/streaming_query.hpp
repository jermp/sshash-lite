#pragma once

#include "../dictionary.hpp"
#include "../util.hpp"

#include "../gz/zip_stream.hpp"
#include "streaming_query_canonical_parsing.hpp"
#include "streaming_query_regular_parsing.hpp"

namespace sshash {

template <typename Query, bool FASTQ>
streaming_query_report streaming_query_from_fastx_file(dictionary const* dict, std::istream& is,
                                                       double threshold) {
    /*

    We assume the query file is well-formed, that is:
    - there are exactly 2 lines per query in case of .fasta file
      (one header line, one read line);
    - there are exactly 4 lines per query in case of .fastq file
      (one header line, one read line, '+' line, quality score line').

    If the flag "SSHASH_QUERY_VERBOSE_OUTPUT" is defined
    (compile the code with -DSSHASH_QUERY_VERBOSE_OUTPUT=On),
    then we print the pair (read-id, answer) for each read.

    */

    streaming_query_report report;
    std::string line;
    uint64_t k = dict->k();
    Query query(dict);

    while (!is.eof()) {
        query.start();
        std::getline(is, line);  // skip first header line
        std::getline(is, line);
        if (line.size() < k) continue;

        uint64_t num_positive_kmers_in_read = 0;
        uint64_t num_negative_kmers_in_read = 0;
        uint64_t num_kmers_in_read = line.size() - k + 1;
        uint64_t min_num_positive_kmers_in_read =
            static_cast<double>(num_kmers_in_read) * threshold;
        uint64_t max_num_negative_kmers_in_read =
            num_kmers_in_read - min_num_positive_kmers_in_read;

#ifdef SSHASH_QUERY_VERBOSE_OUTPUT
        std::cerr << report.num_reads << '\t';
#endif

        for (uint64_t i = 0; i != num_kmers_in_read; ++i) {
            char const* kmer = line.data() + i;
            bool answer = query.is_member(kmer);
            if (answer) {
                num_positive_kmers_in_read += 1;
                report.num_positive_kmers += 1;
                if (num_positive_kmers_in_read >= min_num_positive_kmers_in_read) {
                    report.num_positive_reads += 1;
#ifdef SSHASH_QUERY_VERBOSE_OUTPUT
                    std::cerr << '1';
                    std::cerr << '\n';
#endif
                    break;
                }
            } else {
                num_negative_kmers_in_read += 1;
                if (num_negative_kmers_in_read >= max_num_negative_kmers_in_read) {
#ifdef SSHASH_QUERY_VERBOSE_OUTPUT
                    std::cerr << '0';
                    std::cerr << '\n';
#endif
                    break;
                }
            }
        }

        report.num_reads += 1;
        report.num_kmers += num_kmers_in_read;

        if constexpr (FASTQ) {
            std::getline(is, line);  // skip '+'
            std::getline(is, line);  // skip score
        }
    }

    report.num_searches = query.num_searches();
    report.num_extensions = query.num_extensions();

    return report;
}

streaming_query_report dictionary::streaming_query_from_file(std::string const& filename,
                                                             double threshold) const {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    streaming_query_report report;

    if (util::ends_with(filename, ".fa.gz") or util::ends_with(filename, ".fasta.gz")) {
        zip_istream zis(is);
        if (canonicalized()) {
            report = streaming_query_from_fastx_file<streaming_query_canonical_parsing, false>(
                this, zis, threshold);
        } else {
            report = streaming_query_from_fastx_file<streaming_query_regular_parsing, false>(
                this, zis, threshold);
        }
    } else if (util::ends_with(filename, ".fq.gz") or util::ends_with(filename, ".fastq.gz")) {
        zip_istream zis(is);
        if (canonicalized()) {
            report = streaming_query_from_fastx_file<streaming_query_canonical_parsing, true>(
                this, zis, threshold);
        } else {
            report = streaming_query_from_fastx_file<streaming_query_regular_parsing, true>(
                this, zis, threshold);
        }
    } else if (util::ends_with(filename, ".fa") or util::ends_with(filename, ".fasta")) {
        if (canonicalized()) {
            report = streaming_query_from_fastx_file<streaming_query_canonical_parsing, false>(
                this, is, threshold);
        } else {
            report = streaming_query_from_fastx_file<streaming_query_regular_parsing, false>(
                this, is, threshold);
        }
    } else if (util::ends_with(filename, ".fq") or util::ends_with(filename, ".fastq")) {
        if (canonicalized()) {
            report = streaming_query_from_fastx_file<streaming_query_canonical_parsing, true>(
                this, is, threshold);
        } else {
            report = streaming_query_from_fastx_file<streaming_query_regular_parsing, true>(
                this, is, threshold);
        }
    }

    is.close();
    return report;
}

}  // namespace sshash
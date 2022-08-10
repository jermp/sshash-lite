#pragma once

#include "../dictionary.hpp"
#include "../util.hpp"

#include "../gz/zip_stream.hpp"
#include "streaming_query_canonical_parsing.hpp"
#include "streaming_query_regular_parsing.hpp"

namespace sshash {

inline bool is_positive_read(uint64_t num_positive_kmers_in_read, uint64_t num_kmers_in_read,
                             float threshold) {
    return num_positive_kmers_in_read > num_kmers_in_read * threshold;
}

template <typename Query>
streaming_query_report streaming_query_from_fasta_file_multiline(dictionary const* dict,
                                                                 std::istream& is,
                                                                 float threshold) {
    streaming_query_report report;
    buffered_lines_iterator it(is);
    std::string buffer;
    uint64_t k = dict->k();
    Query query(dict);
    query.start();
    uint64_t num_positive_kmers_in_read = 0;
    uint64_t num_kmers_in_read = 0;
    while (!it.eof()) {
        bool empty_line_was_read = it.fill_buffer(buffer);
        num_kmers_in_read += buffer.size() - k + 1;
        for (uint64_t i = 0; i != buffer.size() - k + 1; ++i) {
            char const* kmer = buffer.data() + i;
            bool answer = query.is_member(kmer);
            report.num_positive_kmers += answer;
            num_positive_kmers_in_read += answer;
        }
        report.num_kmers += buffer.size() - k + 1;
        if (empty_line_was_read) { /* re-start the kmers' buffer */
            report.num_reads += 1;
            report.num_positive_reads +=
                is_positive_read(num_positive_kmers_in_read, num_kmers_in_read, threshold);
            num_positive_kmers_in_read = 0;
            num_kmers_in_read = 0;
            buffer.clear();
            query.start();
        } else {
            if (buffer.size() > k - 1) {
                std::copy(buffer.data() + buffer.size() - k + 1, buffer.data() + buffer.size(),
                          buffer.data());
                buffer.resize(k - 1);
            }
        }
    }
    report.num_searches = query.num_searches();
    report.num_extensions = query.num_extensions();
    return report;
}

template <typename Query>
streaming_query_report streaming_query_from_fasta_file(dictionary const* dict, std::istream& is,
                                                       float threshold) {
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
        uint64_t num_kmers_in_read = line.size() - k + 1;
        for (uint64_t i = 0; i != num_kmers_in_read; ++i) {
            char const* kmer = line.data() + i;
            bool answer = query.is_member(kmer);
            report.num_positive_kmers += answer;
            num_positive_kmers_in_read += answer;
        }
        report.num_reads += 1;
        report.num_positive_reads +=
            is_positive_read(num_positive_kmers_in_read, num_kmers_in_read, threshold);
        report.num_kmers += num_kmers_in_read;
    }
    report.num_searches = query.num_searches();
    report.num_extensions = query.num_extensions();
    return report;
}

template <typename Query>
streaming_query_report streaming_query_from_fastq_file(dictionary const* dict, std::istream& is,
                                                       float threshold) {
    streaming_query_report report;
    std::string line;
    uint64_t k = dict->k();
    Query query(dict);
    while (!is.eof()) {
        query.start();
        /* We assume the file is well-formed, i.e., there are exactly 4 lines per read. */
        std::getline(is, line);  // skip first header line
        std::getline(is, line);
        if (line.size() < k) continue;
        uint64_t num_positive_kmers_in_read = 0;
        uint64_t num_kmers_in_read = line.size() - k + 1;
        for (uint64_t i = 0; i != num_kmers_in_read; ++i) {
            char const* kmer = line.data() + i;
            bool answer = query.is_member(kmer);
            report.num_positive_kmers += answer;
            num_positive_kmers_in_read += answer;
        }
        report.num_reads += 1;
        report.num_positive_reads +=
            is_positive_read(num_positive_kmers_in_read, num_kmers_in_read, threshold);
        report.num_kmers += num_kmers_in_read;
        std::getline(is, line);  // skip '+'
        std::getline(is, line);  // skip score
    }
    report.num_searches = query.num_searches();
    report.num_extensions = query.num_extensions();
    return report;
}

template <typename Query>
streaming_query_report streaming_query_from_fasta_file(dictionary const* dict, std::istream& is,
                                                       bool multiline, float threshold) {
    if (multiline) return streaming_query_from_fasta_file_multiline<Query>(dict, is, threshold);
    return streaming_query_from_fasta_file<Query>(dict, is, threshold);
}

streaming_query_report dictionary::streaming_query_from_file(std::string const& filename,
                                                             bool multiline,
                                                             float threshold) const {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    streaming_query_report report;

    if (util::ends_with(filename, ".fa.gz") or util::ends_with(filename, ".fasta.gz")) {
        zip_istream zis(is);
        if (canonicalized()) {
            report = streaming_query_from_fasta_file<streaming_query_canonical_parsing>(
                this, zis, multiline, threshold);
        } else {
            report = streaming_query_from_fasta_file<streaming_query_regular_parsing>(
                this, zis, multiline, threshold);
        }
    } else if (util::ends_with(filename, ".fq.gz") or util::ends_with(filename, ".fastq.gz")) {
        if (multiline) {
            std::cout << "==> Warning: option 'multiline' is only valid for FASTA files, not FASTQ."
                      << std::endl;
        }
        zip_istream zis(is);
        if (canonicalized()) {
            report = streaming_query_from_fastq_file<streaming_query_canonical_parsing>(this, zis,
                                                                                        threshold);
        } else {
            report = streaming_query_from_fastq_file<streaming_query_regular_parsing>(this, zis,
                                                                                      threshold);
        }
    } else if (util::ends_with(filename, ".fa") or util::ends_with(filename, ".fasta")) {
        if (canonicalized()) {
            report = streaming_query_from_fasta_file<streaming_query_canonical_parsing>(
                this, is, multiline, threshold);
        } else {
            report = streaming_query_from_fasta_file<streaming_query_regular_parsing>(
                this, is, multiline, threshold);
        }
    } else if (util::ends_with(filename, ".fq") or util::ends_with(filename, ".fastq")) {
        if (multiline) {
            std::cout << "==> Warning: option 'multiline' is only valid for FASTA files, not FASTQ."
                      << std::endl;
        }
        if (canonicalized()) {
            report = streaming_query_from_fastq_file<streaming_query_canonical_parsing>(this, is,
                                                                                        threshold);
        } else {
            report = streaming_query_from_fastq_file<streaming_query_regular_parsing>(this, is,
                                                                                      threshold);
        }
    }

    is.close();
    return report;
}

}  // namespace sshash
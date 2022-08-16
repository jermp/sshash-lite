#pragma once

#include "../include/gz/zip_stream.hpp"

namespace sshash {

void random_kmer(char* kmer, uint64_t k) {
    for (uint64_t i = 0; i != k; ++i) kmer[i] = "ACGT"[rand() % 4];
}

bool check_correctness_membership(std::istream& is, dictionary const& dict) {
    uint64_t k = dict.k();
    uint64_t n = dict.size();

    std::string kmer_str(k, 0);
    std::cout << "checking correctness of membership..." << std::endl;

    std::string line;
    uint64_t num_kmers = 0;

    while (!is.eof()) {
        std::getline(is, line);  // header sequence
        std::getline(is, line);  // DNA sequence
        if (line.size() < k) continue;

        for (uint64_t i = 0; i + k <= line.size(); ++i) {
            assert(util::is_valid(line.data() + i, k));
            uint64_t uint64_kmer = util::string_to_uint64_no_reverse(line.data() + i, k);

            if (num_kmers != 0 and num_kmers % 5000000 == 0) {
                std::cout << "checked " << num_kmers << " kmers" << std::endl;
            }
            if ((num_kmers & 1) == 0) {
                /* transform 50% of the kmers into their reverse complements */
                uint64_kmer = util::compute_reverse_complement(uint64_kmer, k);
            }

            util::uint64_to_string_no_reverse(uint64_kmer, kmer_str.data(), k);
            bool answer = dict.is_member(kmer_str.c_str());

            if (!answer) {
                std::cout << "ERROR: got false but expected true for kmer '" << kmer_str << "'"
                          << std::endl;
            }
            assert(answer);

            ++num_kmers;
        }
    }
    std::cout << "checked " << num_kmers << " kmers" << std::endl;

    std::cout << "EVERYTHING OK!" << std::endl;

    std::cout << "checking correctness of negative membership with random kmers..." << std::endl;
    uint64_t num_queries = std::min<uint64_t>(1000000, n);
    for (uint64_t i = 0; i != num_queries; ++i) {
        random_kmer(kmer_str.data(), k);
        /*
            We could use a std::unordered_set to check if kmer is really absent,
            but that would take much more memory...
        */
        bool answer = dict.is_member(kmer_str.c_str());
        if (answer) std::cerr << "kmer '" << kmer_str << "' found!" << std::endl;
    }

    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

/*
   The input file must be the one the index was built from.
*/
bool check_correctness_membership(dictionary const& dict, std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        good = check_correctness_membership(zis, dict);
    } else {
        good = check_correctness_membership(is, dict);
    }
    is.close();
    return good;
}

}  // namespace sshash
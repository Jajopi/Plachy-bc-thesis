#pragma once

#include <vector>
#include <iostream>
#include <limits>

#include "../kmers.h"
#include "loac_compute.h"

#define DEFAULT_PRECISION 100

template <typename kmer_t, typename size_t_max>
void compute_with_loac(std::vector<kmer_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty, size_t precision){

    std::sort(kMers.begin(), kMers.end());

    auto loac = LeafOnlyAC<kmer_t, size_t_max>(
        kMers, size_t_max(k), complements);

    if (run_penalty == 0) run_penalty = log2(kMers.size());
    if (precision == 0) precision = DEFAULT_PRECISION;
    if (precision > sizeof(size_t_max) * 8) precision = sizeof(size_t_max) * 8;

    loac.set_search_parameters(run_penalty, precision);
    loac.compute_result();

    os << ">superstring k=" << k << std::endl;
    loac.print_result(os);
}

template <typename kmer_t>
void set_limit_and_compute_with_loac(std::vector<kmer_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty, size_t precision){
    try {
        if (kMers.empty()) {
            throw std::invalid_argument("Input cannot be empty.");
        }

        // Add complementary k-mers.
        size_t n = kMers.size();
        if (complements){
            kMers.resize(n * 2);
            for (size_t i = 0; i < n; ++i) kMers[i + n] = ReverseComplement(kMers[i], k);
        }

        size_t limit = kMers.size();

        if      (limit < (size_t(1) << 15))
            compute_with_loac<kmer_t, uint16_t>(kMers, os, k, complements, run_penalty, precision);
        else if (limit < (size_t(1) << 31))
            compute_with_loac<kmer_t, uint32_t>(kMers, os, k, complements, run_penalty, precision);
        else
            compute_with_loac<kmer_t, uint64_t>(kMers, os, k, complements, run_penalty, precision);
    }
    catch (const std::exception& e){
        std::cerr << std::endl << "Exception was thrown: " << e.what() << std::endl;
    }
}

template <typename kmer_t>
void LOAC(std::vector<kmer_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty = 0, size_t precision = 0);

void LOAC(std::vector<kmer64_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty = 0, size_t precision = 0) {
    set_limit_and_compute_with_loac<kmer64_t>(kMers, os, k, complements, run_penalty, precision);
}
void LOAC(std::vector<kmer128_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty = 0, size_t precision = 0) {
    set_limit_and_compute_with_loac<kmer128_t>(kMers, os, k, complements, run_penalty, precision);
}
void LOAC(std::vector<kmer256_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty = 0, size_t precision = 0) {
    set_limit_and_compute_with_loac<kmer256_t>(kMers, os, k, complements, run_penalty, precision);
}

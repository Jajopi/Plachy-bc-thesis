#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <list>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "../kmers.h"

#include "optimize.h"

template <typename kmer_t>
size_t decode_and_print_indexes(const std::vector<kmer_t>& kMers, const std::vector<size_t>& indexes, std::ostream& os, size_t k){
    size_t total_length = 0;
    size_t run_count = 1;
    
    kmer_t actual_kmer = kMers[indexes[0]], new_kmer = 0;
    for (size_t i = 1; i < indexes.size(); ++i){
        new_kmer = kMers[indexes[i]];
        
        size_t ov = max_overlap_length(actual_kmer, new_kmer, k);
        print_kmer_masked(actual_kmer, k, os, size_t(k - ov));
        
        total_length += k - ov;
        if (ov < k - 1) ++run_count;

        actual_kmer = new_kmer;
    }
    print_kmer_masked(new_kmer, k, os, k);
    
    total_length += k;

    // std::cerr << std::endl;
    // std::cerr << "Total length: " << total_length << std::endl;
    // std::cerr << "Run count: " << run_count << std::endl;

    return total_length;
}

template <typename kmer_t>
void GlobalILP(std::vector<kmer_t>& kMers, std::ostream& os, size_t k, bool complements) {
    if (kMers.empty()) {
		throw std::invalid_argument("Input cannot be empty.");
	}
    
    // Add complementary k-mers.
    size_t N = kMers.size();
    if (complements){
        kMers.resize(N * 2);
        for (size_t i = 0; i < N; ++i) {
            kMers[i + N] = ReverseComplement(kMers[i], k);
        }
    }

    std::vector<size_t> indexes = compute_indexes(kMers, k, complements);

    os << ">superstring k=" << k << std::endl;
    size_t total_length = decode_and_print_indexes(kMers, indexes, os, k);
    
    std::cerr << total_length << " / " << kMers.size() * k / (complements ? 2 : 1) << std::endl;
}


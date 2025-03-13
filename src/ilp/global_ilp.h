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
#include "distance_functions.h"

template <typename kmer_t>
void GlobalILP(std::vector<kmer_t>& kMers, std::ostream& os, size_t k, bool complements) {
    if (kMers.empty()) {
		throw std::invalid_argument("input cannot be empty");
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

    std::cerr << indexes.size() << std::endl;

    size_t total_length = decode_and_print_indexes(kMers, indexes, os, k);
    
    std::cerr << total_length << " / " << kMers.size() * k / (complements ? 2 : 1) << std::endl;
}


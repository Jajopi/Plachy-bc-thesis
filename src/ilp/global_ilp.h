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
void GlobalILP(std::vector<kmer_t>& kMers, std::ostream& os, size_t k, bool complements, size_t lower_bound = 0) {
    if (kMers.empty()) {
		throw std::invalid_argument("input cannot be empty");
	}
    
    // Add complementary k-mers.
    size_t n = kMers.size();
    kMers.resize(n * (1 + complements));
    if (complements){
        for (size_t i = 0; i < n; ++i) {
            kMers[i + n] = ReverseComplement(kMers[i], k);
        }
    }

    std::vector<size_t> indexes = optimize_indexes(kMers, trivial_distance, k, complements, lower_bound);

    size_t total_length = decode_and_print_indexes(kMers, indexes, os, k);
    os << std::endl;
    
    std::cout << total_length << " / " << kMers.size() * k / (complements ? 2 : 1) << std::endl;
}


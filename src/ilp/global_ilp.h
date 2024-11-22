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
void GlobalILP(std::vector<kmer_t> kMers, std::ostream& of, size_t k) {
	of << k << " " << kMers.size() << std::endl;

    if (kMers.empty()) {
		throw std::invalid_argument("input cannot be empty");
	}
    /*size_t n = kMers.size();
    // Add complementary k-mers.
    kMers.resize(n * (1 + complements));
    if (complements) for (size_t i = 0; i < n; ++i) {
        kMers[i + n] = ReverseComplement(kMers[i]);
    }*/

    std::vector<size_t> indexes = optimize_indexes(kMers, trivial_distance, k);

    for (auto index : indexes){
        of << index << " ";
    }
    of << std::endl;
}


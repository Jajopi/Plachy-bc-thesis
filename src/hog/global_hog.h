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

#include "hog.h"

template <typename kmer_t, typename size_t_max>
void compute_with_hog(size_t_max type_num,
        std::vector<kmer_t>& kMers, std::ostream& os, size_t k, bool complements){
    auto hog = HOGConstructer<kmer_t, size_t_max>(kMers, k); // TODO handle complements
    hog.create();

    os << "ok";
}

template <typename kmer_t>
void GlobalHOG(std::vector<kmer_t>& kMers, std::ostream& os, size_t k, bool complements) {
    if (kMers.empty()) {
		throw std::invalid_argument("Input cannot be empty");
	}

    // Add complementary k-mers.
    size_t n = kMers.size();
    kMers.resize(n * (1 + complements));
    if (complements){
        for (size_t i = 0; i < n; ++i) {
            kMers[i + n] = ReverseComplement(kMers[i], k);
        }
    }
    size_t limit = kMers.size() * k;
    if (limit < (1 << 16)) compute_with_hog((uint16_t)0, kMers, os, k, complements);
    else if (limit < (1 << 32)) compute_with_hog((uint32_t)0, kMers, os, k, complements);
    else compute_with_hog((uint64_t)0, kMers, os, k, complements);
}


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

template <typename kmer_t>
int trivial_distance(const std::vector<kmer_t>& kmers, size_t k, size_t index1, size_t index2){
    if (index1 >= kmers.size() || index2 >= kmers.size())
        throw std::invalid_argument("index " + std::to_string(std::max(index1, index2)) + " outside of bounds [0, " + std::to_string(kmers.size() - 1) + "]");
    
    kmer_t kmer1 = kmers[index1], kmer2 = kmers[index2];
    return k - max_overlap_length(kmer1, kmer2, k);
}

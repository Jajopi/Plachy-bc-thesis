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
        std::cout << index << " ";
    }
    std::cout << std::endl;

    for (size_t i = 0; i < indexes.size(); ++i){
        print_kmer(kMers[indexes[i]], k);
        std::cout << " ";
    }
    std::cout << std::endl;

    size_t total_length = 0;
    kmer_t actual_kmer = kMers[indexes[0]];
    for (size_t c = 0; c < k; c++){
        of << NucleotideAtIndex(actual_kmer, k, c);
        total_length++;
    }
    kmer_t last_kmer;
    for (size_t i = 1; i < indexes.size(); ++i){
        last_kmer = actual_kmer;
        actual_kmer = kMers[indexes[i]];

        size_t ov = compute_overlap(last_kmer, actual_kmer, k);

        for (size_t c = ov; c < k; ++c){
            of << NucleotideAtIndex(actual_kmer, k, c);
            total_length++;
        }
    }
    of << std::endl;

    std::cout << total_length << " / " << kMers.size() * k << std::endl;
}


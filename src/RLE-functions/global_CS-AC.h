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

#include "CS-AC.h"

template <typename kmer_t, typename size_t_max>
void compute_with_cs_ac(size_t_max type_num,
        std::vector<kmer_t>& kMers, std::ostream& os, size_t k, bool complements){
    auto csac = CuttedSortedAC<kmer_t, size_t_max>(kMers, k, k / 2); // TODO handle complements
    csac.construct_graph();
    //csac.print_topological();
    csac.construct_leaf_ranges();
    //csac.print_sorted();
    csac.print_stats();

    //auto indexes = csac.compute_ordering();
    //decode_and_print_indexes(kMers, indexes, os, k);
}

template <typename kmer_t>
void GlobalCS_AC(std::vector<kmer_t>& kMers, std::ostream& os, size_t k, bool complements) {
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
    if (limit < (size_t(1) << 15))
        compute_with_cs_ac(uint16_t(0), kMers, os, k, complements);
    else if (limit < (size_t(1) << 31))
        compute_with_cs_ac(uint32_t(0), kMers, os, k, complements);
    else
        compute_with_cs_ac(uint64_t(0), kMers, os, k, complements);
}


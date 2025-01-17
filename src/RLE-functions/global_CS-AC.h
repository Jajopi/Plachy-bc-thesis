#pragma once

#include <vector>
#include <iostream>
#include <limits>

#include "../kmers.h"

#include "CS-AC.h"

#define RESERVED_MEMORY_GB 6

// getting available memory according to https://stackoverflow.com/questions/2513505/how-to-get-available-memory-c-g
#ifdef USING_WINDOWS // NOT TESTED
    #include <windows.h>
    size_t getTotalSystemMemory()
    {
        MEMORYSTATUSEX status;
        status.dwLength = sizeof(status);
        GlobalMemoryStatusEx(&status);
        return status.ullTotalPhys;
    }
#else // NOT TESTED EITHER except for my own computer
    #include <unistd.h>
    size_t getTotalSystemMemory()
    {
        size_t pages = sysconf(_SC_PHYS_PAGES);
        size_t page_size = sysconf(_SC_PAGE_SIZE);
        return pages * page_size;
    }
#endif

template <typename kmer_t, typename size_t_max>
size_t_max compute_max_depth(size_t_max number_size_type, std::vector<kmer_t>& kMers, bool complements){
    size_t available_memory = getTotalSystemMemory() - RESERVED_MEMORY_GB * size_t(1 << 30);
    std::cout << "Available memory: " << available_memory << std::endl;

    size_t memory_used_per_kmer = sizeof(kmer_t) // kMers
                                + sizeof(kmer_t) + sizeof(size_t_max) // failures
                                + sizeof(CS_AC_Node<size_t_max>) * 2 // top two levels of nodes
                                + sizeof(CS_AC_Node<size_t_max>); // current_nodes
    size_t memory_used_by_kmers = memory_used_per_kmer * kMers.size() * (1 + complements);
    
    size_t expected_memory_used_by_solving = 0; // TODO add
    
    if (available_memory < memory_used_by_kmers + expected_memory_used_by_solving){
        throw std::invalid_argument("Not enough memory for computation with complements!");
    }

    size_t available_nodes = (available_memory - memory_used_by_kmers - expected_memory_used_by_solving) / sizeof(CS_AC_Node<size_t_max>);
    size_t available_depth = available_nodes / kMers.size() + 2; // plus top two levels

    return std::min(available_depth,
                    size_t(std::numeric_limits<size_t_max>::max()));
}

template <typename kmer_t, typename size_t_max>
void compute_with_cs_ac(size_t_max number_size_type,
        std::vector<kmer_t>& kMers, std::ostream& os, size_t k, bool complements){
    
    size_t_max max_depth = compute_max_depth(number_size_type, kMers, complements);
    std::cout << "Maximal available depth: " << max_depth << std::endl;
    size_t_max depth_cutoff = 0;
    if (max_depth < k) depth_cutoff = k - max_depth;

    auto csac = CuttedSortedAC<kmer_t, size_t_max>(kMers, k, depth_cutoff, complements);
    csac.construct_graph();
    //csac.print_topological();
    //csac.construct_leaf_ranges();
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


#pragma once

#include <vector>
#include <iostream>
#include <limits>

#include "../kmers.h"

#include "CS-AC.h"

#define RESERVED_MEMORY_GB 4

// getting available memory according to https://stackoverflow.com/questions/2513505/how-to-get-available-memory-c-g
#ifdef USING_WINDOWS // NOT TESTED
    #include <windows.h>
    size_t getTotalSystemMemory(){
        MEMORYSTATUSEX status;
        status.dwLength = sizeof(status);
        GlobalMemoryStatusEx(&status);
        return status.ullTotalPhys;
    }
#else // NOT TESTED EITHER except for my own computer
    #include <unistd.h>
    size_t getTotalSystemMemory(){
        size_t pages = sysconf(_SC_PHYS_PAGES);
        size_t page_size = sysconf(_SC_PAGE_SIZE);
        return pages * page_size;
    }
#endif

template <typename kmer_t, typename size_t_max>
size_t_max compute_max_depth(size_t_max number_size_type, std::vector<kmer_t>& kMers, bool complements){
    size_t available_memory = getTotalSystemMemory() - RESERVED_MEMORY_GB * size_t(1 << 30);
    std::cout << "Available memory: " << available_memory << std::endl;

    size_t storing_memory_reserved_per_kmer = sizeof(kmer_t);
    
    size_t constructing_and_searching_memory_reserved_per_kmer = std::max((
            2 * sizeof(size_t_max)       // failures
        ),(
            3 * sizeof(size_t_max)       // hq
            + sizeof(size_t_max)         // previous
            + sizeof(size_t_max)         // unionfind
            + sizeof(size_t_max)         // indexes - may reuse space after hq?
        ));
    
    size_t memory_reserved_per_kmer = storing_memory_reserved_per_kmer + constructing_and_searching_memory_reserved_per_kmer;
    size_t memory_reserved_for_kmers = memory_reserved_per_kmer * kMers.size() * (1 + complements);
    
    if (available_memory < memory_reserved_for_kmers){
        throw std::invalid_argument("Not enough memory for computation with complements!");
    }

    size_t available_nodes = (available_memory - memory_reserved_for_kmers) / sizeof(CS_AC_Node<size_t_max, K_SIZE_CONSTANT>);
    size_t available_depth = available_nodes / kMers.size();

    if (available_depth < 3){ // at least two top levels + current_nodes
        throw std::invalid_argument("Not enough memory for computation!");
    }

    return std::min(available_depth,
                    size_t(std::numeric_limits<size_t_max>::max()));
}

template <typename kmer_t, typename size_t_max>
void compute_with_cs_ac(
        std::vector<kmer_t>& kMers, std::ostream& os, size_t_max k, bool complements){
    
    size_t_max max_depth = compute_max_depth(size_t_max(0), kMers, complements);
    std::cout << "Maximal available depth: " << max_depth << std::endl;
    size_t_max depth_cutoff = 0;
    if (max_depth < k) depth_cutoff = k - max_depth;

    std::sort(kMers.begin(), kMers.end()); // TODO move into class

    auto csac = CuttedSortedAC<kmer_t, size_t_max>(kMers, k, depth_cutoff, complements);
    csac.construct_graph();
    csac.print_stats();
    // csac.print_topological();
    csac.convert_to_searchable_representation();
    //csac.print_sorted();

    auto indexes = csac.compute_indexes(k - depth_cutoff);
    size_t_max length = decode_and_print_indexes(kMers, indexes, os, k);//, true, true);
    os << std::endl;
    std::cout << length << std::endl;
}

template <typename kmer_t>
void GlobalCS_AC(std::vector<kmer_t>& kMers, std::ostream& os, size_t k, bool complements) {
    try {
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
            compute_with_cs_ac(kMers, os, uint16_t(k), complements);
        else if (limit < (size_t(1) << 31))
            compute_with_cs_ac(kMers, os, uint32_t(k), complements);
        else
            compute_with_cs_ac(kMers, os, uint64_t(k), complements);
    }
    catch (const std::exception& e){
        std::cerr << "Exception was thrown: " << e.what() << std::endl;
    }
}


#pragma once

#include <vector>
#include <iostream>
#include <limits>

#include "../kmers.h"

#include "CS-AC_construction.h"
#include "CS-AC_searching.h"

#define RESERVED_MEMORY_GB 4
#define DEFAULT_PRECISION 15
// #define DEBUG_FAST_COMPILATION

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

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
size_t compute_max_depth(size_t kmer_count){
    size_t available_memory = getTotalSystemMemory() - (RESERVED_MEMORY_GB * size_t(1 << 30));
    std::cerr << "Available memory: " << available_memory << std::endl;

    size_t storing_memory_per_kmer = sizeof(kmer_t);
    
    size_t constructing_memory_per_kmer = 2 * sizeof(size_t_max) * 2;   // failures, copied when shortening
    size_t searching_memory_per_kmer = sizeof(size_t_max)               // unionfind
                                     + sizeof(size_t_max) * 3           // stack, 3 numbers per element, up to N elements
                                     + sizeof(size_t_max)               // backtracks
                                     + sizeof(size_t_max)               // backtrack_indexes
                                     + sizeof(size_t_max)               // previous
                                     + sizeof(size_t_max);              // remaining_priorities
    
    size_t memory_reserved_per_kmer = storing_memory_per_kmer + std::max(constructing_memory_per_kmer, searching_memory_per_kmer);
    size_t memory_reserved_for_kmers = memory_reserved_per_kmer * kmer_count;

    // std::cerr << "Memory reserved for kmers: " << memory_reserved_for_kmers << " ( " << kmer_count << " )" << std::endl;
    
    if (available_memory < memory_reserved_for_kmers){
        throw std::invalid_argument("Not enough memory for storing the data.");
    }

    size_t available_memory_for_nodes = (available_memory - memory_reserved_for_kmers);
    size_t available_nodes = (available_memory_for_nodes) / sizeof(CS_AC_Node<size_t_max, K_BIT_SIZE>);
    size_t available_depth = available_nodes / kmer_count;
    // std::cerr << "Available nodes: " << available_nodes << ' ' << sizeof(CS_AC_Node<size_t_max, K_BIT_SIZE>) << std::endl;

    if (available_depth < 2){ // at least top two levels
        throw std::invalid_argument("Not enough memory for computation.");
    }
    std::cerr << "Maximal available depth: " << available_depth << std::endl;
    return available_depth;
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
void compute_with_cs_ac(std::vector<kmer_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty, size_t precision){

    size_t max_depth = compute_max_depth<kmer_t, size_t_max, K_BIT_SIZE>(kMers.size());

    size_t depth_cutoff = 0;
    if (max_depth < k) depth_cutoff = k - max_depth;
    size_t practical_depth = max_depth;

    size_t_max exponent = 1;
    while (size_t(1 << (2 * exponent)) < kMers.size()) ++exponent;
    if (exponent > depth_cutoff){
        depth_cutoff = 0;
        ++practical_depth;
    }

    auto csac = CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>(
        kMers, size_t_max(k), size_t_max(depth_cutoff), size_t_max(practical_depth), complements);
    csac.construct_graph();
    csac.convert_to_searchable_representation();

    if (run_penalty == 0) run_penalty = log2(kMers.size());
    if (precision == 0) precision = DEFAULT_PRECISION;
    if (precision >= sizeof(size_t_max) * 8) precision = std::numeric_limits<size_t_max>::max();
    csac.set_search_parameters(run_penalty, 1, precision);

    csac.compute_result();
    csac.print_result(os);
}

template <typename kmer_t, size_t K_BIT_SIZE>
void set_limit_and_compute_with_cs_ac(std::vector<kmer_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty, size_t precision){
    try {
        if (kMers.empty()) {
            throw std::invalid_argument("Input cannot be empty.");
        }

        // Add complementary k-mers.
        size_t n = kMers.size();
        kMers.resize(n * (1 + complements));
        if (complements){
            for (size_t i = 0; i < n; ++i) {
                kMers[i + n] = ReverseComplement(kMers[i], k);
            }
        }

        size_t limit = kMers.size() * (size_t(1) << K_BIT_SIZE);

        if      (limit < (size_t(1) << 15)) compute_with_cs_ac<kmer_t, uint16_t, K_BIT_SIZE>(kMers, os, k, complements, run_penalty, precision);
        else if (limit < (size_t(1) << 31)) compute_with_cs_ac<kmer_t, uint32_t, K_BIT_SIZE>(kMers, os, k, complements, run_penalty, precision);
        else                                compute_with_cs_ac<kmer_t, uint64_t, K_BIT_SIZE>(kMers, os, k, complements, run_penalty, precision);
    }
    catch (const std::exception& e){
        std::cerr << std::endl << "Exception was thrown: " << e.what() << std::endl;
    }
}

template <typename kmer_t>
void GlobalCS_AC(std::vector<kmer_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty = 0, size_t precision = 0);

void GlobalCS_AC(std::vector<kmer32_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty = 0, size_t precision = 0) {
    set_limit_and_compute_with_cs_ac<kmer32_t, 4>(kMers, os, k, complements, run_penalty, precision);
}
void GlobalCS_AC(std::vector<kmer64_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty = 0, size_t precision = 0) {
    set_limit_and_compute_with_cs_ac<kmer64_t, 5>(kMers, os, k, complements, run_penalty, precision);
}
void GlobalCS_AC(std::vector<kmer128_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty = 0, size_t precision = 0) {
    set_limit_and_compute_with_cs_ac<kmer128_t, 6>(kMers, os, k, complements, run_penalty, precision);
}
void GlobalCS_AC(std::vector<kmer256_t>& kMers, std::ostream& os, size_t k,
        bool complements, size_t run_penalty = 0, size_t precision = 0) {
    set_limit_and_compute_with_cs_ac<kmer256_t, 7>(kMers, os, k, complements, run_penalty, precision);
}

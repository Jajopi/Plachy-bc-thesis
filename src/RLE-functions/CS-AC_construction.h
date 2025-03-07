#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>
#include <map>

#include "../kmers.h"

#define MAX_DEPTH_WIDTH 3

template<typename size_t_max>
class UnionFind {
    std::vector<size_t_max> roots;
    size_t_max component_count;
public:
    UnionFind(size_t_max size) : roots(size), component_count(size) {
        for (size_t_max i = 0; i < size; ++i) roots[i] = i;
    };
    UnionFind() = default;
    
    inline size_t_max find(size_t_max x) {
        size_t_max root = roots[x];
        if (roots[root] == root) return root;
        while (roots[root] != root) root = roots[root];
        
        while (x != root){
            size_t_max new_x = roots[x];
            roots[x] = root;
            x = new_x;
        }
        return root;
    }

    inline bool are_connected(size_t_max x, size_t_max y){
        return find(x) == find(y);
    }

    inline void connect(size_t_max to, size_t_max from){ // Second one points to the first one - points to the begining of a chain
        if (are_connected(from, to)) return;
        roots[from] = to;
        --component_count;
    }

    inline size_t_max count() const { return component_count; };
};


template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
class CuttedSortedAC {
    static inline size_t_max INVALID_LEAF() { return std::numeric_limits<size_t_max>::max(); };
    static inline size_t_max INVALID_NODE() { return std::numeric_limits<size_t_max>::max(); };

    std::vector<kmer_t> kMers;  // sorted
    std::vector<size_t_max> chains;
    // std::vector<size_t_max> chain_beginnings;
    
    size_t_max K;               // kmer-length
    size_t_max N;               // number of kmers, number of leaves
    size_t_max DEPTH_CUTOFF;    // first not-to-be-reached depth
    size_t_max PRACTICAL_DEPTH;
    size_t_max SEARCH_CUTOFF;   // fraction of uncompleted leaves to switch to fast mode
    bool COMPLEMENTS;           // whether or not complements are used
    
    bool COMPUTED_RESULT = false;
    size_t_max RUN_PENALTY = INVALID_NODE();
    
    std::vector<std::tuple<size_t_max, size_t_max, size_t_max, size_t_max>> stack;
    UnionFind<size_t_max> components;
    std::vector<size_t_max> complements;
    std::vector<size_t_max> backtracks;
    std::vector<size_t_max> backtrack_indexes;
    std::vector<size_t_max> previous;
    std::vector<size_t_max> next;
    std::vector<size_t_max> remaining_priorities;
    std::vector<bool> used;

    void sort_kmers();

    bool try_complete_leaf(size_t_max leaf_to_connect, size_t_max priority_drop_limit);
    void push_failure_of_node_into_stack(size_t_max priority, size_t_max leaf_index, size_t_max chain_depth, size_t_max priority_drop_limit, size_t_max last_leaf);
    void squeeze_uncompleted_leaves(std::vector<size_t_max>& unclompleted_leaves);
    void set_backtrack_path_for_leaf(size_t_max origin_leaf, size_t_max next_leaf);
    
    size_t_max find_complement_kmer_index(size_t_max kmer_index);
    size_t_max find_failure_kmer_index(kmer_t kmer, size_t_max suffix = 0);
public:
    CuttedSortedAC(const std::vector<kmer_t>& kmers, size_t_max K, size_t_max DEPTH_CUTOFF = 0, bool complements = false) :
        kMers(kmers), K(K), N(kMers.size()), DEPTH_CUTOFF(DEPTH_CUTOFF), PRACTICAL_DEPTH(K - DEPTH_CUTOFF),
        COMPLEMENTS(complements) { sort_kmers(); };
    CuttedSortedAC(std::vector<kmer_t>&& kmers, size_t_max K, size_t_max DEPTH_CUTOFF = 0, bool complements = false) :
        kMers(std::move(kmers)), K(K), N(kMers.size()), DEPTH_CUTOFF(DEPTH_CUTOFF), PRACTICAL_DEPTH(K - DEPTH_CUTOFF),
        COMPLEMENTS(complements) { sort_kmers(); };

    void construct_graph();

    void convert_to_searchable_representation();
    void set_search_parameters(size_t_max run_penalty,
                               size_t_max precision);
    void compute_result();
    size_t_max print_result(std::ostream& os);

    std::ostream& LOG_STREAM = std::cerr;
};

// Constructing

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::construct_graph() {
    if (!chains.empty()){
        throw std::invalid_argument("Graph has already been constructed.");
    }

    chains.reserve(N * (K - DEPTH_CUTOFF));
    // chain_beginnings.resize(N + 1);

    for (size_t_max i = 0; i < N; ++i){
        kmer_t current = kMers[i];
        // chain_beginnings[i] = chains.size();

        for (size_t_max depth = K - 1; depth > DEPTH_CUTOFF; --depth){
            size_t_max failure_index = find_failure_kmer_index(current, depth);
            chains.push_back(failure_index);
        }

        if (COMPLEMENTS) complements[i] = find_complement_kmer_index(i);
    }

    // chain_beginnings[N] = chains.size();
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::set_search_parameters(
            size_t_max run_penalty, size_t_max precision) {
    RUN_PENALTY = run_penalty;
    if (precision >= sizeof(size_t_max) * 8) SEARCH_CUTOFF = 0; // Infinite precision (no early ending)
    else SEARCH_CUTOFF = N / (1 << precision);

    LOG_STREAM << "Run penalty: " << run_penalty << std::endl;
    // LOG_STREAM << "Precision: " << precision << std::endl;
}

// Sorting

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::sort_kmers() {
    std::sort(kMers.begin(), kMers.end());
}

// Internal functions

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline size_t_max CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::find_complement_kmer_index(size_t_max kmer_index) {
    kmer_t complement = ReverseComplement(kMers[kmer_index], K);
    
    size_t_max begin = 0, end = N - 1;
    while (begin < end){
        size_t_max middle = (begin + end) / 2;
        if (kMers[middle] == complement) return middle;
        else if (kMers[middle] < complement) begin = middle + 1;
        else end = middle - 1;
    }
    return begin; // Complement should always be present exactly once
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline size_t_max CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::find_failure_kmer_index(kmer_t kmer, size_t_max suffix) {
    kmer_t searched = BitSuffix(kmer, suffix);
    
    size_t_max begin = 0, end = N - 1;
    while (begin < end){
        size_t_max middle = (begin + end) / 2;
        kmer_t current = BitSuffix(kMers[middle], suffix);

        if (current == searched) return middle;
        else if (current < searched) begin = middle + 1;
        else end = middle - 1;
    }
    if (BitSuffix(kMers[begin], suffix) != searched) return INVALID_LEAF();
    return begin;
}

#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>

#include "../kmers.h"
#include "CS-AC_construction.h"

// Searching

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::compute_result() {
    if (!CONVERTED_TO_SEARCHABLE){
        throw std::invalid_argument("Graph has not been converted to searchable yet.");
    }
    if (NEW_RUN_SCORE == INVALID_NODE()){
        throw std::invalid_argument("Mandatory parameter NEW_RUN_SCORE was not set.");
    }
    if (BASE_EXTENSION_SCORE == INVALID_NODE()){
        throw std::invalid_argument("Mandatory parameter BASE_EXTENSION_SCORE was not set.");
    }

    LOG_STREAM << "Computing indexes..." << std::endl;

    components = UnionFind(N);
    backtracks.reserve(N); // Chains for leaves where backtracking is needed
    backtrack_indexes.resize(N, INVALID_LEAF()); // Indexes into backtracks for each leaf
    previous.resize(N);
    remaining_priorities.resize(N);

    std::vector<size_t_max> uncompleted_leaves(N);
    for (size_t_max i = 0; i < N; ++i) uncompleted_leaves[i] = i;

    size_t_max max_priority_drop = (K - DEPTH_CUTOFF) * BASE_EXTENSION_SCORE + NEW_RUN_SCORE;

    size_t_max remaining_iterations = max_priority_drop / BASE_EXTENSION_SCORE;
    LOG_STREAM << "Leaves: " << N << ", iterations: " << remaining_iterations << std::endl;
    LOG_STREAM << std::setw(10) << N << ' ' << std::setw(4) << remaining_iterations; LOG_STREAM.flush();

    for (size_t_max priority_drop_limit = BASE_EXTENSION_SCORE;
                    priority_drop_limit <= max_priority_drop;
                    priority_drop_limit += BASE_EXTENSION_SCORE){
        size_t_max uncompleted_leaf_count = uncompleted_leaves.size();
        LOG_STREAM << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::setw(10) << uncompleted_leaf_count
            << ' ' << std::setw(4) << remaining_iterations--; LOG_STREAM.flush();
        
        if (uncompleted_leaf_count < SEARCH_CUTOFF) break;
        
        for (size_t_max i = 0; i < uncompleted_leaf_count; ++i){
            size_t_max leaf_index = uncompleted_leaves[i];

            if (COMPLEMENTS && nodes[leaf_index].complement_completed()){
                uncompleted_leaves[i] = INVALID_LEAF();
                continue;
            }

            bool result = try_complete_leaf(leaf_index, priority_drop_limit);
            
            if (result){
                uncompleted_leaves[i] = INVALID_LEAF();
                if (COMPLEMENTS) nodes[nodes[leaf_index].complement_index].set_complement_completed();
            }
        }

        squeeze_uncompleted_leaves(uncompleted_leaves);
    }

    LOG_STREAM << std::endl;

    COMPUTED_RESULT = true;
}

// Internal functions

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline bool CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::try_complete_leaf(
        size_t_max leaf_index, size_t_max priority_drop_limit) {
    
    if (priority_drop_limit == BASE_EXTENSION_SCORE){
        Node& leaf_node = nodes[leaf_index];
        Node& failure_node = nodes[leaf_node.failure];
        if (failure_node.depth() < K - 1) return false;
        
        size_t_max leaf_range_end = nodes[leaf_node.failure + 1].original_leaf_range_begin();
        if (leaf_range_end < failure_node.leaf_range_begin) leaf_range_end = N;
        
        for (size_t_max i = failure_node.leaf_range_begin; i < leaf_range_end; ++i){
            if (!nodes[i].used() && !components.are_connected(leaf_index, i)){
                nodes[i].set_used();
                leaf_node.set_next(i);
                components.connect(leaf_index, i); // Second one pointing to the first one

                if (COMPLEMENTS){
                    nodes[nodes[i].complement_index].set_used();
                    nodes[leaf_node.complement_index].set_complement_completed();
                }
                return true;
            }
        }
        return false;
    }
    
    size_t_max root_node_index = nodes.size() - 1;
    Node& leaf_node = nodes[leaf_index];

    size_t_max minimal_priority_limit = INVALID_NODE() - priority_drop_limit;
    hq.clear(); // Stores priority, current node_index, last_leaf index
    push_failure_of_node_into_hq(INVALID_NODE(), leaf_index, minimal_priority_limit, leaf_index); // priority, node, limit, last_leaf

    while (!hq.empty()){
        auto t = hq.back(); hq.pop_back();
        size_t_max priority = std::get<0>(t);
        size_t_max node_index = std::get<1>(t);
        size_t_max last_leaf = std::get<2>(t);
        Node& node = nodes[node_index];
        
        if (node_index == root_node_index) return false;

        size_t_max leaf_range_end = nodes[node_index + 1].original_leaf_range_begin();
        if (leaf_range_end < node.leaf_range_begin) leaf_range_end = N;

        while (node.leaf_range_begin != leaf_range_end &&
                (nodes[node.leaf_range_begin].used() || // Was already used as next for other leaf / complement was used
                (COMPLEMENTS && nodes[node.leaf_range_begin].complement_completed()) || // Complement was completed
                components.are_connected(leaf_index, node.leaf_range_begin))){ // Is from the same chain as current leaf trying to be completed
            ++node.leaf_range_begin;
        }
        if (node.leaf_range_begin != leaf_range_end){ // We found at least one suitable leaf to complete the current one
            nodes[node.leaf_range_begin].set_used();
            leaf_node.set_next(node.leaf_range_begin);
            components.connect(leaf_index, node.leaf_range_begin);

            if (COMPLEMENTS){
                nodes[nodes[node.leaf_range_begin].complement_index].set_used();
                nodes[leaf_node.complement_index].set_complement_completed();
            }

            if (last_leaf != leaf_index){
                previous[node.leaf_range_begin] = last_leaf;
                set_backtrack_path_for_leaf(leaf_index, node.leaf_range_begin);
            }

            return true;
        }
        else {
            node.leaf_range_begin = node.original_leaf_range_begin(); // Reset the counter -- improves precision AND speed 
        }

        for (size_t_max i = node.original_leaf_range_begin(); i < leaf_range_end; ++i){ // Search also through already used leaves
            if (i == leaf_index) continue;
            if (nodes[i].resigned()) continue; // But not the ones pointing nowhere

            if (remaining_priorities[i] >= priority - minimal_priority_limit){
                continue; // There is nothing to be found with remaining steps
            }
            remaining_priorities[i] = priority - minimal_priority_limit;
            
            previous[i] = last_leaf;
            push_failure_of_node_into_hq(priority, i, minimal_priority_limit, i); // Add failure of that leaf
        }

        push_failure_of_node_into_hq(priority, node_index, minimal_priority_limit, last_leaf); // Add failure of current node
    }

    return false;
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::push_failure_of_node_into_hq(
        size_t_max priority, size_t_max node_index, size_t_max minimal_priority_limit, size_t_max last_leaf) {
    Node& node = nodes[node_index];
    Node& failure_node = nodes[node.failure];

    size_t_max node_depth = (node_index < N) ? K : node.depth(); // Leaves have always depth of K and complement_index is stored in depth field
    
    size_t_max failure_priority = priority - (node_depth - failure_node.depth()) * BASE_EXTENSION_SCORE;
    if (node_depth == K - 1 || (node_depth == K && failure_node.depth() < K - 1)){ // Run will be interrupted
        failure_priority -= NEW_RUN_SCORE;
        if (minimal_priority_limit == 0) return; // Special value indicating no interruptions are expected
    }

    if (failure_priority < minimal_priority_limit) return;

    hq.push_back(std::make_tuple(failure_priority, node.failure, last_leaf));
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::squeeze_uncompleted_leaves(std::vector<size_t_max> &unclompleted_leaves) {
    size_t_max count = unclompleted_leaves.size(), shift = 0;
    
    for (size_t_max i = 0; i < count; ++i){
        if (unclompleted_leaves[i] == INVALID_LEAF()) ++shift;
        else unclompleted_leaves[i - shift] = unclompleted_leaves[i];
    }

    unclompleted_leaves.resize(count - shift);
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::set_backtrack_path_for_leaf(size_t_max origin_leaf, size_t_max next_leaf) {
    size_t_max actual = next_leaf;
    while (actual != origin_leaf){
        backtracks.push_back(actual);
        actual = previous[actual];
    }
    backtrack_indexes[origin_leaf] = backtracks.size() - 1;
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline size_t_max CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::print_result(std::ostream& os) {
    if (!COMPUTED_RESULT){
        throw std::invalid_argument("Result has not been computed yet.");
    }

    LOG_STREAM << "Printing..." << std::endl;
    size_t_max total_length = 0;
    size_t_max run_count = 1;

    bool first = true;
    size_t_max actual = INVALID_LEAF();
    kmer_t last_kmer = 0;
    for (size_t_max i = 0; i < N; ++i){
        if (components.find(i) != i) continue;
        if (COMPLEMENTS && nodes[i].complement_completed()) continue;

        actual = i;
        while (actual != INVALID_LEAF()){
            if (first){
                first = false;
                last_kmer = kMers[i];
                actual = nodes[i].next();
                continue;
            }

            kmer_t actual_kmer = kMers[actual];
            size_t_max ov = compute_max_overlap(last_kmer, actual_kmer, K);
            print_kmer_masked(last_kmer, K, os, size_t_max(K - ov));

            last_kmer = actual_kmer;
            total_length += K - ov;
            if (ov < K - 1) ++run_count;

            if (backtrack_indexes[actual] != INVALID_LEAF()){
                size_t_max backtrack_index = backtrack_indexes[actual];
                size_t_max actual_backtrack = backtracks[backtrack_index];
                size_t_max next = nodes[actual].next();
                while (actual_backtrack != next){
                    kmer_t actual_kmer = kMers[actual_backtrack];
                    size_t_max ov = compute_max_overlap(last_kmer, actual_kmer, K);

                    print_kmer_masked(last_kmer, K, os, size_t_max(K - ov));

                    last_kmer = actual_kmer;
                    total_length += K - ov;
                    if (ov < K - 1) ++run_count;

                    --backtrack_index;
                    actual_backtrack = backtracks[backtrack_index];
                }
            }

            actual = nodes[actual].next();
        }        
    }
    print_kmer_masked(last_kmer, K, os);
    total_length += K;

    LOG_STREAM << std::endl;
    LOG_STREAM << "Total length: " << total_length << std::endl;
    LOG_STREAM << "Run count: " << run_count << std::endl;

    return total_length * BASE_EXTENSION_SCORE + run_count * NEW_RUN_SCORE;
}

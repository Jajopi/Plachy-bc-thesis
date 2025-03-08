#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>
#include <random>

#include "../kmers.h"
#include "CS-AC_construction.h"

#define RANDOM_SEED 0
#define MAX_COUNT_WIDTH 12
#define MAX_ITERS_WIDTH 3

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::compute_result() {
    if (RUN_PENALTY == INVALID_NODE()){
        throw std::invalid_argument("Mandatory parameter RUN_PENALTY was not set.");
    }

    components = UnionFind(N);
    backtracks.reserve(N); // Chains for leaves where backtracking is needed
    backtrack_indexes.resize(N, INVALID_NODE()); // Indexes into backtracks for each leaf
    previous.resize(N, INVALID_NODE());
    next.resize(N, INVALID_NODE());
    remaining_priorities.resize(N, 0);
    used.resize(N, false);

    std::vector<size_t_max> uncompleted_leaves(N);
    for (size_t_max i = 0; i < N; ++i) uncompleted_leaves[i] = i;
    std::shuffle(uncompleted_leaves.begin(), uncompleted_leaves.end(), std::default_random_engine(RANDOM_SEED));
    size_t_max next_preffered_leaf = INVALID_NODE();
    
    size_t_max max_priority_drop = (K - 1) + RUN_PENALTY;
    
    size_t_max remaining_iterations = max_priority_drop;
    LOG_STREAM << std::setw(MAX_COUNT_WIDTH) << N << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations << std::endl;
    LOG_STREAM << std::setw(MAX_COUNT_WIDTH) << N << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations; LOG_STREAM.flush();
    
    for (size_t_max priority_drop_limit = 1;
        priority_drop_limit <= max_priority_drop;
        ++priority_drop_limit){
            size_t_max uncompleted_leaf_count = uncompleted_leaves.size();
            LOG_STREAM << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::setw(MAX_COUNT_WIDTH) << uncompleted_leaf_count
            << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations--; LOG_STREAM.flush();
            
        if (uncompleted_leaf_count < SEARCH_CUTOFF) break;
        
        for (size_t_max i = 0; i < uncompleted_leaf_count; ++i){
            size_t_max leaf_index = uncompleted_leaves[i];
            
            if (next_preffered_leaf != INVALID_NODE()){
                leaf_index = next_preffered_leaf;
                --i;
                next_preffered_leaf = INVALID_NODE();

                if (next[leaf_index] != INVALID_NODE()) continue;
    
                bool result = try_complete_leaf(leaf_index, priority_drop_limit);
                if (result){
                    next_preffered_leaf = next[leaf_index];
                }
            }
            else{
                if (leaf_index == INVALID_NODE()) continue;
                if (next[leaf_index] != INVALID_NODE()){
                    uncompleted_leaves[i] = INVALID_NODE();
                    continue;
                }
    
                bool result = try_complete_leaf(leaf_index, priority_drop_limit);
                if (result){
                    next_preffered_leaf = next[leaf_index];
                    uncompleted_leaves[i] = INVALID_NODE();
                }
            }
            
        }

        squeeze_uncompleted_leaves(uncompleted_leaves);
    }

    LOG_STREAM << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::setw(MAX_COUNT_WIDTH) << uncompleted_leaves.size()
        << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations; LOG_STREAM.flush();

    COMPUTED_RESULT = true;
}

// Internal functions

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline bool CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::try_complete_leaf(
        size_t_max leaf_to_complete, size_t_max priority_drop_limit) {

    print_kmer(kMers[leaf_to_complete], K, LOG_STREAM, K);
    LOG_STREAM << ' ' << leaf_to_complete << ' ' << priority_drop_limit;

    if (priority_drop_limit == 1){
        size_t_max chain_begin = leaf_to_complete * DEPTH;
        size_t_max first_failure_leaf = chains[chain_begin];

        if (first_failure_leaf == INVALID_NODE()){
            LOG_STREAM << std::endl;
            return false;
        }
        size_t_max leaf_complement = COMPLEMENTS ? complements[leaf_to_complete] : INVALID_NODE();

        for (size_t_max i = first_failure_leaf; i < N &&
                BitPrefix(kMers[first_failure_leaf], K, K - 1) == BitPrefix(kMers[i], K, K - 1);
                ++i){
            if ((COMPLEMENTS && leaf_complement == i) ||
                used[i] ||
                components.are_connected(leaf_to_complete, i)) continue;

            used[i] = true;
            next[leaf_to_complete] = i;
            components.connect(leaf_to_complete, i); // Second one pointing at the first one

            if (COMPLEMENTS){ // Connect complements inversely
                used[leaf_complement] = true;
                next[complements[i]] = leaf_complement;
                components.connect(complements[i], leaf_complement);
            }

            LOG_STREAM << " -> " << i << ' '; print_kmer(kMers[i], K, LOG_STREAM, K); LOG_STREAM << std::endl;
            return true;
        }

        LOG_STREAM << std::endl;
        return false;
    }

    stack.clear(); // Stores priority, current leaf_index, current chain depth, last_leaf index
    push_failure_of_node_into_stack(priority_drop_limit, leaf_to_complete, K, leaf_to_complete);

    while (!stack.empty()){
        auto t = stack.back(); stack.pop_back();
        size_t_max priority = std::get<0>(t);
        size_t_max leaf_index = std::get<1>(t);
        size_t_max chain_depth = std::get<2>(t);
        size_t_max last_leaf = std::get<3>(t);
        LOG_STREAM << std::endl << priority << ' ' << leaf_index << ' ' << chain_depth << ' ' << last_leaf << ' ';
        print_kmer(kMers[leaf_index], K, LOG_STREAM, K);

        if (chain_depth == 0){
            continue;
        }

        size_t_max leaf_complement = COMPLEMENTS ? complements[leaf_to_complete] : INVALID_NODE();

        for (size_t_max i = leaf_index; i < N &&
                BitPrefix(kMers[leaf_index], K, chain_depth) == BitPrefix(kMers[i], K, chain_depth);
                ++i){
            if ((COMPLEMENTS && leaf_complement == i) ||
                used[i] ||
                components.are_connected(leaf_to_complete, i)) continue;

            used[i] = true;
            next[leaf_to_complete] = i;
            components.connect(leaf_to_complete, i); // Second one pointing at the first one

            if (COMPLEMENTS){ // Connect complements inversely
                used[leaf_complement] = true;
                next[complements[i]] = leaf_complement;
                components.connect(complements[i], leaf_complement);
            }

            if (last_leaf != leaf_to_complete){
                previous[i] = last_leaf;
                set_backtrack_path_for_leaf(leaf_to_complete, i);
            }

            LOG_STREAM << " -> " << i << ' '; print_kmer(kMers[i], K, LOG_STREAM, K); LOG_STREAM << std::endl;
            return true;
        }

        for (size_t_max i = leaf_index; i < N &&
            BitPrefix(kMers[leaf_index], K, chain_depth) == BitPrefix(kMers[i], K, chain_depth);
            ++i){
            if (i == leaf_index || i == leaf_to_complete ||
                (COMPLEMENTS && i == leaf_complement)) continue;

            if (remaining_priorities[i] >= priority) continue;
            remaining_priorities[i] = priority;

            previous[i] = last_leaf;
            push_failure_of_node_into_stack(priority, i, K, i); // Add failure of that leaf
            LOG_STREAM << ' ' << 'p' << ' ' << i;
        }

        push_failure_of_node_into_stack(priority, leaf_index, chain_depth, last_leaf); // Add failure of current node
    }
    LOG_STREAM << std::endl;
    return false;
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::push_failure_of_node_into_stack(
        size_t_max priority, size_t_max node_index, size_t_max node_depth, size_t_max last_leaf) {
    size_t_max failure_depth = node_depth - 1;
    size_t_max failure_index = INVALID_NODE();

    while (failure_depth > 0){
        if (failure_depth > K - 1 - DEPTH){
            failure_index = chains[node_index * DEPTH + K - 1 - failure_depth];
        }
        else {
            failure_index = find_first_failure_leaf(kMers[node_index], failure_depth);
        }
        
        if (failure_index != INVALID_NODE()) break;
        --failure_depth;
    }
    if (failure_depth == 0) return;
    
    if (priority <= node_depth - failure_depth) return;
    priority -= (node_depth - failure_depth);
    if (node_depth == K - 1 || (node_depth == K && failure_depth < K - 1)){ // Run will be interrupted
        if (priority <= RUN_PENALTY) return;
        priority -= RUN_PENALTY;
    }

    stack.push_back(std::make_tuple(priority, failure_index, failure_depth, last_leaf));
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::squeeze_uncompleted_leaves(std::vector<size_t_max> &uncompleted_leaves) {
    size_t_max count = uncompleted_leaves.size(), shift = 0;
    
    for (size_t_max i = 0; i < count; ++i){
        if (uncompleted_leaves[i] == INVALID_NODE()) ++shift;
        else uncompleted_leaves[i - shift] = uncompleted_leaves[i];
    }

    uncompleted_leaves.resize(count - shift);
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::set_backtrack_path_for_leaf(size_t_max origin_leaf, size_t_max next_leaf) {
    // LOG_STREAM << origin_leaf << ' ' << next_leaf << std::endl;
    size_t_max actual = next_leaf;
    while (actual != origin_leaf){
        backtracks.push_back(actual);
        actual = previous[actual];
        // LOG_STREAM << actual << std::endl;
        if (actual == INVALID_NODE()) break;
    }
    backtrack_indexes[origin_leaf] = backtracks.size() - 1;
    // LOG_STREAM << "x" << std::endl;
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline size_t_max CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::print_result(std::ostream& os) {
    if (!COMPUTED_RESULT){
        throw std::invalid_argument("Result has not been computed yet.");
    }
    // LOG_STREAM << std::endl << "Printing..."; LOG_STREAM.flush();

    size_t_max total_length = 0;
    size_t_max run_count = 1;

    bool first = true;
    size_t_max actual = INVALID_NODE();
    kmer_t last_kmer = 0;
    for (size_t_max i = 0; i < N; ++i){
        if (components.find(i) != i){
            // other_count++;
            continue;
        }
        if (COMPLEMENTS){
            components.connect(i, components.find(complements[i]));
            // Prevent complementary chain from being printed later
        }

        actual = i;
        while (actual != INVALID_NODE()){
            if (first){
                first = false;
                last_kmer = kMers[i];
                actual = next[i];
                continue;
            }

            kmer_t actual_kmer = kMers[actual];
            size_t_max ov = max_overlap_length(last_kmer, actual_kmer, K);
            print_kmer_masked(last_kmer, K, os, size_t_max(K - ov));

            last_kmer = actual_kmer;
            total_length += K - ov;
            if (ov < K - 1) ++run_count;

            if (backtrack_indexes[actual] != INVALID_NODE()){
                size_t_max backtrack_index = backtrack_indexes[actual];
                size_t_max actual_backtrack = backtracks[backtrack_index];
                size_t_max next_node = next[actual];
                while (actual_backtrack != next_node){
                    // other_count++;
                    kmer_t actual_kmer = kMers[actual_backtrack];
                    size_t_max ov = max_overlap_length(last_kmer, actual_kmer, K);
                    print_kmer_masked(last_kmer, K, os, size_t_max(K - ov));

                    last_kmer = actual_kmer;
                    total_length += K - ov;
                    if (ov < K - 1) ++run_count;

                    --backtrack_index;
                    actual_backtrack = backtracks[backtrack_index];
                }
            }

            if (next[actual] == actual) break;
            actual = next[actual];
        }        
    }
    print_kmer_masked(last_kmer, K, os, K);
    total_length += K;

    LOG_STREAM << std::endl;

    return total_length + run_count * RUN_PENALTY;
}

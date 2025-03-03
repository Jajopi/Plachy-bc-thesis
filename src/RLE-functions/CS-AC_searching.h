#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>
#include <random>

#include "../kmers.h"
#include "CS-AC_construction.h"

// Searching

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::compute_result() {
    if (!CONVERTED_TO_SEARCHABLE){
        throw std::invalid_argument("Graph has not been converted to searchable yet.");
    }
    if (RUN_PENALTY == INVALID_NODE()){
        throw std::invalid_argument("Mandatory parameter RUN_PENALTY was not set.");
    }
    if (EXTENSION_PENALTY == INVALID_NODE()){
        throw std::invalid_argument("Mandatory parameter EXTENSION_PENALTY was not set.");
    }

    // LOG_STREAM << "Computing indexes..." << std::endl;

    components = UnionFind(N);
    backtracks.reserve(N); // Chains for leaves where backtracking is needed
    backtrack_indexes.resize(N, INVALID_LEAF()); // Indexes into backtracks for each leaf
    previous.resize(N);
    remaining_priorities.resize(N);

    std::vector<size_t_max> uncompleted_leaves(N);
    for (size_t_max i = 0; i < N; ++i) uncompleted_leaves[i] = i;
    std::shuffle(uncompleted_leaves.begin(), uncompleted_leaves.end(), std::default_random_engine(0));
    size_t_max next_preffered_leaf = INVALID_LEAF();

    size_t_max max_priority_drop = (K - 1) * EXTENSION_PENALTY + RUN_PENALTY;

    size_t_max remaining_iterations = max_priority_drop / EXTENSION_PENALTY;
    LOG_STREAM << "Leaves, iterations:" << std::endl;
    LOG_STREAM << std::setw(12) << N << ' ' << std::setw(4) << remaining_iterations << std::endl;
    LOG_STREAM << std::setw(12) << N << ' ' << std::setw(4) << remaining_iterations; LOG_STREAM.flush();

    for (size_t_max priority_drop_limit = EXTENSION_PENALTY;
                    priority_drop_limit <= max_priority_drop;
                    priority_drop_limit += EXTENSION_PENALTY){
        size_t_max uncompleted_leaf_count = uncompleted_leaves.size();
        LOG_STREAM << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::setw(12) << uncompleted_leaf_count
            << ' ' << std::setw(4) << remaining_iterations--; LOG_STREAM.flush();
        
        if (uncompleted_leaf_count < SEARCH_CUTOFF) break;
        
        for (size_t_max i = 0; i < uncompleted_leaf_count; ++i){
            
            size_t_max leaf_index = uncompleted_leaves[i];
            
            if (next_preffered_leaf != INVALID_LEAF()){
                leaf_index = next_preffered_leaf;
                --i;
                next_preffered_leaf = INVALID_LEAF();

                if (nodes[leaf_index].next() != INVALID_LEAF()) continue;
                // if (COMPLEMENTS && nodes[leaf_index].complement_chosen()) continue;
    
                bool result = try_complete_leaf(leaf_index, priority_drop_limit);
                if (result){
                    next_preffered_leaf = nodes[leaf_index].next();
                }
            }
            else{
                if (leaf_index == INVALID_LEAF()) continue;
                if (nodes[leaf_index].next() != INVALID_LEAF()){
                    uncompleted_leaves[i] = INVALID_LEAF();
                    continue;
                }
                // if (COMPLEMENTS && nodes[leaf_index].complement_chosen()) continue;
    
                bool result = try_complete_leaf(leaf_index, priority_drop_limit);
                if (result){
                    next_preffered_leaf = nodes[leaf_index].next();
                    uncompleted_leaves[i] = INVALID_LEAF();
                }
            }
            
        }

        squeeze_uncompleted_leaves(uncompleted_leaves);
    }

    LOG_STREAM << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::setw(12) << uncompleted_leaves.size()
        << ' ' << std::setw(4) << remaining_iterations << std::endl;

    COMPUTED_RESULT = true;
}

// Internal functions

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline bool CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::try_complete_leaf(
        size_t_max leaf_index, size_t_max priority_drop_limit) {
    
    if (priority_drop_limit == EXTENSION_PENALTY){
        Node& leaf_node = nodes[leaf_index];
        Node& failure_node = nodes[leaf_node.failure];
        if (failure_node.depth() < K - 1) return false;
        
        size_t_max leaf_range_end = nodes[leaf_node.failure + 1].original_leaf_range_begin();
        if (leaf_range_end < failure_node.leaf_range_begin) leaf_range_end = N;
        
        for (size_t_max i = failure_node.leaf_range_begin; i < leaf_range_end; ++i){
            if (nodes[i].used() ||
                components.are_connected(leaf_index, i) ||
                (COMPLEMENTS && leaf_node.complement_index == i))
                // (COMPLEMENTS && nodes[i].complement_chosen()))
                    continue;
            nodes[i].set_used();
            leaf_node.set_next(i);
            components.connect(leaf_index, i); // Second one pointing at the first one
            // print_kmer(kMers[leaf_index], K, LOG_STREAM, K); LOG_STREAM << '-';
            // print_kmer(kMers[i], K, LOG_STREAM, K); LOG_STREAM << std::endl;

            if (COMPLEMENTS){
                // nodes[nodes[i].complement_index].set_complement_chosen();
                // nodes[leaf_node.complement_index].set_complement_chosen();
                nodes[leaf_node.complement_index].set_used();
                nodes[nodes[i].complement_index].set_next(leaf_node.complement_index);
                components.connect(nodes[i].complement_index, leaf_node.complement_index);
            }
            // LOG_STREAM << i << std::endl;
            return true;
        }
        // LOG_STREAM << "x" << std::endl;
        return false;
    }
    
    size_t_max root_node_index = nodes.size() - 1;
    Node& leaf_node = nodes[leaf_index];

    size_t_max minimal_priority_limit = INVALID_NODE() - priority_drop_limit;
    stack.clear(); // Stores priority, current node_index, last_leaf index
    push_failure_of_node_into_stack(INVALID_NODE(), leaf_index, minimal_priority_limit, leaf_index); // priority, node, limit, last_leaf

    while (!stack.empty()){
        auto t = stack.back(); stack.pop_back();
        size_t_max priority = std::get<0>(t);
        size_t_max node_index = std::get<1>(t);
        size_t_max last_leaf = std::get<2>(t);
        Node& node = nodes[node_index];
        
        if (node_index == root_node_index) return false;

        size_t_max leaf_range_end = nodes[node_index + 1].original_leaf_range_begin();
        if (leaf_range_end < node.leaf_range_begin) leaf_range_end = N;
        size_t_max old_begin = node.leaf_range_begin;

        while (node.leaf_range_begin != leaf_range_end &&
                (nodes[node.leaf_range_begin].used() || // Was already used as next for other leaf
                components.are_connected(leaf_index, node.leaf_range_begin) || // Is from the same chain as current leaf trying to be completed
                (COMPLEMENTS && leaf_node.complement_index == node.leaf_range_begin))){ // Complement was chosen to be used before
            ++node.leaf_range_begin;
        }
        if (node.leaf_range_begin != leaf_range_end){ // We found at least one suitable leaf to complete the current one
            nodes[node.leaf_range_begin].set_used();
            leaf_node.set_next(node.leaf_range_begin);
            components.connect(leaf_index, node.leaf_range_begin);
            // print_kmer(kMers[leaf_index], K, LOG_STREAM, K); LOG_STREAM << ' ';
            // print_kmer(kMers[node.leaf_range_begin], K, LOG_STREAM, K); LOG_STREAM << std::endl;

            if (COMPLEMENTS){
                nodes[leaf_node.complement_index].set_used();
                nodes[nodes[node.leaf_range_begin].complement_index].set_next(leaf_node.complement_index);
                components.connect(nodes[node.leaf_range_begin].complement_index, leaf_node.complement_index);
            }

            if (last_leaf != leaf_index){
                previous[node.leaf_range_begin] = last_leaf;
                set_backtrack_path_for_leaf(leaf_index, node.leaf_range_begin);
            }

            // LOG_STREAM << node.leaf_range_begin << std::endl;
            return true;
        }
        else {
            if (old_begin != node.original_leaf_range_begin()){ // If we didn't skip all the leaves in single search
                node.leaf_range_begin = node.original_leaf_range_begin(); // Reset the counter -- improves precision AND speed
            }
        }

        for (size_t_max i = node.original_leaf_range_begin(); i < leaf_range_end; ++i){ // Search also through already used leaves
            if (i == leaf_index) continue;

            if (remaining_priorities[i] >= priority - minimal_priority_limit){
                continue; // There is nothing to be found with remaining steps
            }
            remaining_priorities[i] = priority - minimal_priority_limit;
            
            previous[i] = last_leaf;
            push_failure_of_node_into_stack(priority, i, minimal_priority_limit, i); // Add failure of that leaf
        }

        push_failure_of_node_into_stack(priority, node_index, minimal_priority_limit, last_leaf); // Add failure of current node
    }

    // LOG_STREAM << "x" << std::endl;
    return false;
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::push_failure_of_node_into_stack(
        size_t_max priority, size_t_max node_index, size_t_max minimal_priority_limit, size_t_max last_leaf) {
    Node& node = nodes[node_index];
    Node& failure_node = nodes[node.failure];

    size_t_max node_depth = (node_index < N) ? K : node.depth(); // Leaves have always depth of K and complement_index is stored in depth field
    
    size_t_max failure_priority = priority - (node_depth - failure_node.depth()) * EXTENSION_PENALTY;
    if (node_depth == K - 1 || (node_depth == K && failure_node.depth() < K - 1)){ // Run will be interrupted
        failure_priority -= RUN_PENALTY;
        if (minimal_priority_limit == 0) return; // Special value indicating no interruptions are expected
    }

    if (failure_priority < minimal_priority_limit) return;

    stack.push_back(std::make_tuple(failure_priority, node.failure, last_leaf));
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::squeeze_uncompleted_leaves(std::vector<size_t_max> &uncompleted_leaves) {
    size_t_max count = uncompleted_leaves.size(), shift = 0;
    
    for (size_t_max i = 0; i < count; ++i){
        if (uncompleted_leaves[i] == INVALID_LEAF()) ++shift;
        // || (
                // COMPLEMENTS && nodes[uncompleted_leaves[i]].complement_chosen()))
        else uncompleted_leaves[i - shift] = uncompleted_leaves[i];
    }

    uncompleted_leaves.resize(count - shift);
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

    // LOG_STREAM << "Printing..." << std::endl;
    size_t_max total_length = 0;
    size_t_max run_count = 1;
    // size_t_max count = 0, other_count = 0;

    bool first = true;
    size_t_max actual = INVALID_LEAF();
    kmer_t last_kmer = 0;
    for (size_t_max i = 0; i < N; ++i){
        if (components.find(i) != i){
            // other_count++;
            continue;
        }
        if (COMPLEMENTS){
            components.connect(i, components.find(nodes[i].complement_index));
            // Prevent complementary chain from being printed later
        }
        // if (COMPLEMENTS && nodes[i].complement_chosen()){
        //     // count++;
        //     continue;
        // }
        // other_count++;

        actual = i;
        while (actual != INVALID_LEAF()){
            // other_count++;
            // if (COMPLEMENTS && nodes[actual].complement_chosen()){
            //     // count++; break;
            // }

            if (first){
                first = false;
                last_kmer = kMers[i];
                actual = nodes[i].next();
                continue;
            }

            kmer_t actual_kmer = kMers[actual];
            size_t_max ov = max_overlap_length(last_kmer, actual_kmer, K);
            print_kmer_masked(last_kmer, K, os, size_t_max(K - ov));

            last_kmer = actual_kmer;
            total_length += K - ov;
            if (ov < K - 1) ++run_count;

            if (backtrack_indexes[actual] != INVALID_LEAF()){
                size_t_max backtrack_index = backtrack_indexes[actual];
                size_t_max actual_backtrack = backtracks[backtrack_index];
                size_t_max next = nodes[actual].next();
                while (actual_backtrack != next){
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

            actual = nodes[actual].next();
        }        
    }
    print_kmer_masked(last_kmer, K, os, K);
    total_length += K;

    // LOG_STREAM << std::endl;
    // LOG_STREAM << "Total length: " << total_length << std::endl;
    // LOG_STREAM << "Run count: " << run_count << std::endl;
    // /*if (COMPLEMENTS)*/ LOG_STREAM << count << ' ' << other_count << std::endl;

    return total_length * EXTENSION_PENALTY + run_count * RUN_PENALTY;
}

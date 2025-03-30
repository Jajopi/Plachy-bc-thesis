#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <random>
#include <queue>

#include "../kmers.h"
#include "../unionfind.h"

#define RANDOM_SEED 0
#define MAX_COUNT_WIDTH 12
#define MAX_ITERS_WIDTH 3

template<typename kmer_t, typename size_t_max>
class FailureIndex {
    const std::vector<kmer_t> &kMers;
    size_t_max N;
    size_t_max K;
    
    size_t_max FIRST_ROW_COUNT;
    std::vector<size_t_max> first_rows;
    std::vector<std::unordered_map<size_t_max, size_t_max>> cached;
    size_t_max CACHE_LIMIT;

    std::vector<size_t_max> search_speedup;
    size_t_max SPEEDUP_DEPTH;

    static inline size_t_max INVALID_NODE() { return std::numeric_limits<size_t_max>::max(); };

    void construct_index(){
        size_t_max N = kMers.size();
        
        SPEEDUP_DEPTH = log2(N) / 2;
        size_t_max speedup_size = (size_t_max(1) << 2 * SPEEDUP_DEPTH);
        search_speedup.resize(speedup_size + 1);
        size_t_max index = 0;
        for (kmer_t k = 0; k < kmer_t(speedup_size); ++k){
            search_speedup[k] = index;
            while (index < N && BitPrefix(kMers[index], K, SPEEDUP_DEPTH) == k) ++index;
        }
        search_speedup[speedup_size] = N;

        FIRST_ROW_COUNT = log2(K);
        first_rows.resize(N * FIRST_ROW_COUNT);
        for (size_t_max row = 0; row < FIRST_ROW_COUNT; ++row){
            for (size_t_max i = 0; i < N; ++i){
                if (i > 0 && BitSuffix(kMers[i], K - row - 1) == BitSuffix(kMers[i - 1], K - row - 1)){
                    first_rows[row * N + i] = first_rows[row * N + i - 1];
                    continue;
                }
                first_rows[row * N + i] = binary_search(i, K - row - 1);
            }
        }

        CACHE_LIMIT = N / (K - FIRST_ROW_COUNT);
        cached.resize(K - FIRST_ROW_COUNT);
    }

    size_t_max binary_search(size_t_max index, size_t_max depth){
        kmer_t searched = BitSuffix(kMers[index], depth);
        size_t_max begin = 0, end = N - 1;

        if (depth >= SPEEDUP_DEPTH){
            size_t_max speedup_index = BitPrefix(searched, depth, SPEEDUP_DEPTH);
            begin = search_speedup[speedup_index];
            end = search_speedup[speedup_index + 1];
        }

        while (begin < end){
            size_t_max middle = (begin + end) / 2;
            kmer_t current = BitPrefix(kMers[middle], K, depth);

            if (current == searched) end = middle;
            else if (current < searched) begin = middle + 1;
            else end = middle - 1;
        }
        return (BitPrefix(kMers[begin], K, depth) != searched) ? INVALID_NODE() : begin;
    }
public:
    FailureIndex(const std::vector<kmer_t> &kmers, size_t_max k) :
            kMers(kmers), N(kmers.size()), K(k) {
        std::cerr << "Constructing index..."; std::cerr.flush();
        construct_index();
        std::cerr << std::endl;
    }

    size_t_max find_first_failure_leaf(size_t_max index, size_t_max depth){
        if (depth >= (K - FIRST_ROW_COUNT)){
            return first_rows[(K - depth - 1) * N + index];
        }

        auto c = cached[K - FIRST_ROW_COUNT - depth - 1];
        
        auto it = c.find(index);
        if (it != c.end()){
            return it->second;
        }
        
        size_t_max result = binary_search(index, depth);
        
        c.emplace(std::make_pair(index, result));

        if (c.size() >= CACHE_LIMIT){
            std::cerr << ' ' << c.size() << std::endl;
            c.clear();
        }

        return result;
    }
};


template <typename kmer_t, typename size_t_max>
class LeafOnlyAC {
    static inline size_t_max INVALID_NODE() { return std::numeric_limits<size_t_max>::max(); };

    std::vector<kmer_t> kMers;  // sorted
    
    size_t_max K;               // kmer-length
    size_t_max N;               // number of kmers, number of leaves
    bool COMPLEMENTS;           // whether or not complements are used
    bool COMPUTED_RESULT = false;
    
    size_t_max RUN_PENALTY = INVALID_NODE();
    size_t_max SEARCH_CUTOFF;   // fraction of uncompleted leaves to switch to fast mode
    
    FailureIndex<kmer_t, size_t_max> fi;
    
    std::vector<size_t_max> complements;
    UnionFind<size_t_max> components;
    std::vector<std::tuple<size_t_max, size_t_max, size_t_max, size_t_max>> stack;
    // std::queue<std::tuple<size_t_max, size_t_max, size_t_max, size_t_max>> que;
    std::vector<size_t_max> backtracks;
    std::vector<size_t_max> backtrack_indexes;
    std::vector<size_t_max> previous;
    std::vector<size_t_max> next;
    std::vector<size_t_max> remaining_priorities;
    std::vector<size_t_max> skip_to;
    std::vector<bool> used;

    bool try_complete_leaf(size_t_max leaf_to_connect, size_t_max priority_drop_limit);
    void push_failure_of_node_into_stack(size_t_max priority, size_t_max node_index, size_t_max node_depth, size_t_max last_leaf);
    void squeeze_uncompleted_leaves(std::vector<size_t_max>& unclompleted_leaves);
    void set_backtrack_path_for_leaf(size_t_max origin_leaf, size_t_max next_leaf);
public:
    std::ostream& LOG_STREAM = std::cerr;

    LeafOnlyAC(const std::vector<kmer_t>& kmers, size_t_max k, bool complements = false) :
        kMers(kmers), K(k), N(kMers.size()), COMPLEMENTS(complements), fi(kMers, K) {
            find_complements();
        };
    LeafOnlyAC(std::vector<kmer_t>&& kmers, size_t_max k, bool complements = false) :
        kMers(std::move(kmers)), K(k), N(kMers.size()), COMPLEMENTS(complements), fi(kMers, K) {
            find_complements();
        };

    void find_complements();
    void set_search_parameters(size_t_max run_penalty,
                               size_t_max precision);
    void compute_result();
    size_t_max print_result(std::ostream& os);
};

// Constructing

template <typename kmer_t, typename size_t_max>
inline void LeafOnlyAC<kmer_t, size_t_max>::find_complements() {
    if (!COMPLEMENTS) return;

    complements.resize(N);

    std::vector<std::pair<kmer_t, size_t_max>> complement_kmers(N);
    for (size_t_max i = 0; i < N; ++i){
        complement_kmers[i] = std::make_pair(ReverseComplement(kMers[i], K), i);
    }

    std::sort(complement_kmers.begin(), complement_kmers.end());
    
    for (size_t_max i = 0; i < N; ++i){
        complements[i] = complement_kmers[i].second;
    }
    
}

template <typename kmer_t, typename size_t_max>
inline void LeafOnlyAC<kmer_t, size_t_max>::set_search_parameters(
            size_t_max run_penalty, size_t_max precision) {
    RUN_PENALTY = run_penalty;
    if (precision >= sizeof(size_t_max) * 8) SEARCH_CUTOFF = 0; // Infinite precision (no early ending)
    else SEARCH_CUTOFF = N / (1 << precision);

    LOG_STREAM << "Run penalty: " << run_penalty << std::endl;
    // LOG_STREAM << "Precision: " << precision << std::endl;
}

template <typename kmer_t, typename size_t_max>
inline void LeafOnlyAC<kmer_t, size_t_max>::compute_result() {
    if (RUN_PENALTY == INVALID_NODE()){
        throw std::invalid_argument("Mandatory parameter RUN_PENALTY was not set.");
    }

    components = UnionFind(N);
    backtracks.reserve(N); // Chains for leaves where backtracking is needed
    backtrack_indexes.resize(N, INVALID_NODE()); // Indexes into backtracks for each leaf
    previous.resize(N, INVALID_NODE());
    next.resize(N, INVALID_NODE());
    remaining_priorities.resize(N, 0);
    skip_to.resize(N);
    for (size_t_max i = 0; i < N; ++i) skip_to[i] = i + 1;
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
        size_t_max unc = uncompleted_leaf_count;
        
        for (size_t_max i = 0; i < uncompleted_leaf_count; ++i){
            size_t_max leaf_index = uncompleted_leaves[i];
            
            if (next_preffered_leaf != INVALID_NODE()){
                leaf_index = next_preffered_leaf;
                --i;
                next_preffered_leaf = INVALID_NODE();

                if (next[leaf_index] != INVALID_NODE()){
                    continue;
                }
    
                bool result = try_complete_leaf(leaf_index, priority_drop_limit);
                if (result){
                    next_preffered_leaf = next[leaf_index];

                    if (remaining_iterations <= 5){
                        LOG_STREAM << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::setw(MAX_COUNT_WIDTH) << unc--
                        << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations; LOG_STREAM.flush();
                    }
                }
            }
            else{
                if (leaf_index == INVALID_NODE()) continue;
                if (next[leaf_index] != INVALID_NODE()){
                    uncompleted_leaves[i] = INVALID_NODE();
                    if (remaining_iterations <= 5){
                        LOG_STREAM << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::setw(MAX_COUNT_WIDTH) << unc--
                        << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations; LOG_STREAM.flush();
                    }
                    continue;
                }
    
                bool result = try_complete_leaf(leaf_index, priority_drop_limit);
                if (result){
                    next_preffered_leaf = next[leaf_index];
                    uncompleted_leaves[i] = INVALID_NODE();

                    if (remaining_iterations <= 5){
                        LOG_STREAM << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::setw(MAX_COUNT_WIDTH) << unc--
                        << ' ' << std::setw(MAX_ITERS_WIDTH) << remaining_iterations; LOG_STREAM.flush();
                    }
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

template <typename kmer_t, typename size_t_max>
inline bool LeafOnlyAC<kmer_t, size_t_max>::try_complete_leaf(
        size_t_max leaf_to_complete, size_t_max priority_drop_limit) {

    // LOG_STREAM << '\t';
    // print_kmer(kMers[leaf_to_complete], K, LOG_STREAM, K);
    // LOG_STREAM << ' ' << leaf_to_complete << ' ' << priority_drop_limit;
    // LOG_STREAM << std::endl;

    if (priority_drop_limit == 1){
        size_t_max first_failure_leaf = fi.find_first_failure_leaf(leaf_to_complete, K - 1);

        if (first_failure_leaf == INVALID_NODE()){
            // LOG_STREAM << std::endl;
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

            // print_kmer(kMers[i], K, LOG_STREAM, K);
            // LOG_STREAM << '\t' << i << std::endl;

            // LOG_STREAM << " -> " << i << ' '; print_kmer(kMers[i], K, LOG_STREAM, K); LOG_STREAM << std::endl;
            return true;
        }

        // LOG_STREAM << std::endl;
        return false;
    }

    stack.clear();
    // que = {}; // Stores priority, current leaf_index, current chain depth, last_leaf index
    push_failure_of_node_into_stack(priority_drop_limit, leaf_to_complete, K, leaf_to_complete);

    while (!stack.empty()){
        auto t = stack.back(); stack.pop_back();
        size_t_max priority = std::get<0>(t);
        size_t_max leaf_index = std::get<1>(t);
        size_t_max chain_depth = std::get<2>(t);
        size_t_max last_leaf = std::get<3>(t);
        // LOG_STREAM << std::endl << priority << ' ' << leaf_index << ' ' << chain_depth << ' ' << last_leaf << ' ';
        // print_kmer(kMers[leaf_index], K, LOG_STREAM, K);

        // if (chain_depth == 0) continue;

        size_t_max leaf_complement = COMPLEMENTS ? complements[leaf_to_complete] : INVALID_NODE();

        for (size_t_max i = leaf_index; i < N &&
                BitPrefix(kMers[leaf_index], K, chain_depth) == BitPrefix(kMers[i], K, chain_depth);
                ++i){
            // LOG_STREAM << ' ' << i << std::endl;
            if ((COMPLEMENTS && leaf_complement == i) ||
                used[i] ||
                components.are_connected(leaf_to_complete, i)) continue;

            // print_kmer(kMers[i], K, LOG_STREAM, K);
            // LOG_STREAM << '\t' << i << std::endl;

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

            // LOG_STREAM << " -> " << i << ' '; print_kmer(kMers[i], K, LOG_STREAM, K); LOG_STREAM << std::endl;
            return true;
        }

        for (size_t_max i = leaf_index; i < N &&
                BitPrefix(kMers[leaf_index], K, chain_depth) == BitPrefix(kMers[i], K, chain_depth);
                ++i){
            if (i == leaf_to_complete) continue;

            if (priority_drop_limit >= RUN_PENALTY){
                if (remaining_priorities[i] >= priority){
                    std::vector<size_t_max> skipped;
                    size_t_max j = i;
                    while (j < N && remaining_priorities[j] >= priority){
                        skipped.push_back(j);
                        j = skip_to[j];
                    }
                    for (size_t_max s : skipped){
                        skip_to[s] = j;
                    }
                    i = j;
                    continue;
                }
                skip_to[i] = i + 1;
                remaining_priorities[i] = priority;
            }
            else {
                if (remaining_priorities[i] >= priority) continue;
                remaining_priorities[i] = priority;
            }

            previous[i] = last_leaf;
            push_failure_of_node_into_stack(priority, i, K, i); // Add failure of that leaf
        }

        push_failure_of_node_into_stack(priority, last_leaf, chain_depth, last_leaf); // Add failure of current node
    }
    // LOG_STREAM << std::endl;
    return false;
}

template <typename kmer_t, typename size_t_max>
inline void LeafOnlyAC<kmer_t, size_t_max>::push_failure_of_node_into_stack(
        size_t_max priority, size_t_max node_index, size_t_max node_depth, size_t_max last_leaf) {
    size_t_max failure_depth = node_depth - 1;
    size_t_max failure_index = INVALID_NODE();

    while (failure_depth > 0){
        failure_index = fi.find_first_failure_leaf(node_index, failure_depth);
        
        if (failure_index != INVALID_NODE()) break;
        --failure_depth;
    }
    if (failure_depth == 0) return;
    
    if (priority < node_depth - failure_depth) return;
    priority -= (node_depth - failure_depth);

    if (node_depth == K - 1 || (node_depth == K && failure_depth < K - 1)){ // Run will be interrupted
        if (priority < RUN_PENALTY) return;
        priority -= RUN_PENALTY;
    }

    // LOG_STREAM << " p " << failure_index;
    stack.push_back(std::make_tuple(priority, failure_index, failure_depth, last_leaf));
}

template <typename kmer_t, typename size_t_max>
inline void LeafOnlyAC<kmer_t, size_t_max>::squeeze_uncompleted_leaves(std::vector<size_t_max> &uncompleted_leaves) {
    size_t_max count = uncompleted_leaves.size(), shift = 0;
    
    for (size_t_max i = 0; i < count; ++i){
        if (uncompleted_leaves[i] == INVALID_NODE()) ++shift;
        else uncompleted_leaves[i - shift] = uncompleted_leaves[i];
    }

    uncompleted_leaves.resize(count - shift);
}

template <typename kmer_t, typename size_t_max>
inline void LeafOnlyAC<kmer_t, size_t_max>::set_backtrack_path_for_leaf(size_t_max origin_leaf, size_t_max next_leaf) {
    size_t_max last_size = backtracks.size();

    size_t_max actual = next_leaf;
    while (actual != origin_leaf){
        backtracks.push_back(actual);
        // if (actual == previous[actual]) break; // Bad things happen
        actual = previous[actual];
    }
    backtrack_indexes[origin_leaf] = backtracks.size() - 1;

    if (COMPLEMENTS){
        size_t_max new_size = backtracks.size();
        size_t_max count = new_size - last_size;
        for (size_t_max i = 0; i < count; ++i) backtracks.push_back(INVALID_NODE());
        
        size_t_max index = new_size + count - 1;
        size_t_max actual = previous[next_leaf];
        while (index >= last_size + count){
            backtracks[index--] = complements[actual];
            // if (actual == previous[actual]) break; // Bad things happen
            actual = previous[actual];
        }
        backtrack_indexes[complements[next_leaf]] = backtracks.size() - 1;
    }
}

template <typename kmer_t, typename size_t_max>
inline size_t_max LeafOnlyAC<kmer_t, size_t_max>::print_result(std::ostream& os) {
    if (!COMPUTED_RESULT){
        throw std::invalid_argument("Result has not been computed yet.");
    }
    LOG_STREAM << std::endl << "Printing..."; LOG_STREAM.flush();

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

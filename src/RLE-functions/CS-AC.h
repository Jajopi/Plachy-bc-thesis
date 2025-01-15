#pragma once

#include <vector>
#include <queue>
#include <algorithm>
#include <limits>

#include "../kmers.h"


template<typename size_t_max>
struct Node {
    size_t_max kmer_index;
    size_t_max depth;
    size_t_max parent;
    size_t_max failure;
    size_t_max child_range_begin;
    size_t_max child_range_end;

    static inline size_t_max INVALID_NODE() { return std::numeric_limits<size_t_max>::max(); };

    Node(size_t_max kmer_index, size_t_max depth, size_t_max parent) :
        kmer_index(kmer_index), depth(depth), parent(parent),
        failure(INVALID_NODE()),
        child_range_begin(INVALID_NODE()), child_range_end(INVALID_NODE()) {};

    Node(size_t_max kmer_index, size_t_max depth, size_t_max first_child, size_t_max last_child) :
        kmer_index(kmer_index), depth(depth), parent(INVALID_NODE()),
        failure(INVALID_NODE()),
        child_range_begin(first_child), child_range_end(last_child) {};

    void print(std::ostream& os) const;
};

template <typename size_t_max>
inline void Node<size_t_max>::print(std::ostream &os) const {
    os << depth << '/';
    if (kmer_index == INVALID_NODE()) os << "INV"; else os << kmer_index;
    os << ",\tP: ";
    if (parent == INVALID_NODE()) os << "INV"; else os << parent;
    os << ",\tF: ";
    if (failure == INVALID_NODE()) os << "INV"; else os << failure;
    os << ",\tCH: [ ";
    if (child_range_begin == INVALID_NODE()) os << "INV"; else os << child_range_begin;
    os << " - ";
    if (child_range_end == INVALID_NODE()) os << "INV"; else os << child_range_end;
    os << " )";
}


template <typename kmer_t, typename size_t_max>
class CuttedSortedAC {
    std::vector<kmer_t> kMers;
    std::vector<Node<size_t_max>> nodes;
    //std::vector<std::pair<size_t_max, size_t_max>> leaf_intervals;
    
    size_t_max K; // kmer-length
    size_t_max N; // number of kmers, number of leaves
    size_t_max DEPTH_CUTOFF; // last (minimal) depth to be reached
    
    static inline size_t_max INVALID_NODE() { return std::numeric_limits<size_t_max>::max(); };
    static inline kmer_t BIGGEST_KMER() { return std::numeric_limits<kmer_t>::max(); };

    template<typename array_element>
    void sort(std::vector<array_element>& array, size_t_max starting_base = 0);
    template<typename array_element>
    void radix_sort(std::vector<array_element>& array, size_t_max begin, size_t_max end, size_t_max starting_base = 0);

    void resort_and_shorten_failures(std::vector<std::pair<kmer_t, size_t_max>> &array, size_t_max depth);
public:
    CuttedSortedAC(const std::vector<kmer_t>& kMers, size_t_max K, size_t_max DEPTH_CUTOFF = 0) :
        kMers(kMers), K(K), N(kMers.size()), DEPTH_CUTOFF(DEPTH_CUTOFF) {};

    CuttedSortedAC(std::vector<kmer_t>&& kMers, size_t_max K, size_t_max depth_cutoff = 0) :
        kMers(std::move(kMers)), K(K), N(kMers.size()), DEPTH_CUTOFF(DEPTH_CUTOFF) {};

    void construct_graph();
    void construct_leaf_ranges();
    std::vector<size_t_max> compute_ordering(size_t_max new_run_score, size_t_max base_score = 1);

    void print_stats(std::ostream& os = std::cout);
    void print_sorted(std::ostream& os = std::cout);
    void print_topological(std::ostream& os = std::cout);
    void print_topological(std::ostream& os, size_t_max root, size_t_max depth);
};

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::construct_graph() {
    if (!nodes.empty()) throw std::invalid_argument("Graph has already been constructed.");

    sort(kMers);

    nodes.reserve((K - DEPTH_CUTOFF) * N);

    std::vector<Node<size_t_max>> current_nodes;
    current_nodes.reserve(N);
    
    std::vector<std::pair<kmer_t, size_t_max>> failures; // suffix, last index
    failures.resize(N);
    
    // Add leaves
    for (size_t_max i = 0; i < N; ++i){
        current_nodes.emplace_back(size_t_max(N - i - 1), K, i, i + 1);
        failures[i] = std::make_pair(kMers[N - i - 1], i);
    }
    sort(failures);

    // Add other nodes
    for (size_t_max depth = K; depth > DEPTH_CUTOFF; --depth){
        std::cout << depth << ' '; std::cout.flush();
        
        size_t_max last_node_count = nodes.size();
        for (Node node: current_nodes) nodes.push_back(std::move(node));
        current_nodes.clear();
        size_t_max node_count = nodes.size();
        size_t_max new_node_index = node_count;

        resort_and_shorten_failures(failures, depth);
        size_t_max failure_index = failures.size() - 1;

        for (size_t_max i = last_node_count; i < node_count; ++i){
            kmer_t current_prefix = BitPrefix(kMers[nodes[i].kmer_index], K, depth - 1);
            while (failure_index > 0 && current_prefix < failures[failure_index].first) --failure_index;

            bool new_node_on_failure_path = (current_prefix == failures[failure_index].first);

            size_t_max j = i;
            while (j < node_count && current_prefix == BitPrefix(kMers[nodes[j].kmer_index], K, depth - 1)){
                nodes[j].parent = new_node_index;
                ++j;
            }

            current_nodes.emplace_back(nodes[i].kmer_index,
                                       size_t_max(depth - 1),
                                       size_t_max(i),
                                       size_t_max(j));
            
            if (new_node_on_failure_path){
                bool first_similar_failure = true;
                while (current_prefix == failures[failure_index].first){
                    nodes[failures[failure_index].second].failure = new_node_index;
                    
                    if (first_similar_failure){
                        first_similar_failure = false;
                        failures[failure_index].second = new_node_index;
                    }
                    else {
                        failures[failure_index] = std::make_pair(BIGGEST_KMER(), INVALID_NODE()); // Will be last, lost
                    }

                    if (failure_index == 0) break;
                    --failure_index;
                }
            }
            
            ++new_node_index;
            i = j - 1;
        }
    }
    std::cout << std::endl;

    nodes.emplace_back(0, 0, 0, nodes.size());
    size_t_max root_node = nodes.size() - 1;
    for (size_t_max i = root_node - 1; i >= N; --i){
        if (nodes[i].parent == INVALID_NODE()) nodes[i].parent = root_node;
        else break;
    }

    nodes.shrink_to_fit();
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::construct_leaf_ranges() {
    size_t_max node_count = nodes.size();
    for (size_t_max i = N; i < node_count; ++i){
        nodes[i].child_range_begin = nodes[nodes[i].child_range_begin].child_range_begin;
        nodes[i].child_range_end = nodes[nodes[i].child_range_end - 1].child_range_end;
    }
}

// Sorting

template <typename kmer_t, typename size_t_max>
template <typename array_element>
inline void CuttedSortedAC<kmer_t, size_t_max>::sort(std::vector<array_element> &array, size_t_max starting_base) {
    if (array.size() > 2 * size_t_max(1 << (K - starting_base))){
        radix_sort(array, 0, array.size(), starting_base);
    }
    else {
        std::sort(array.begin(), array.end());
    }
}

//in-place radix sort as in https://stackoverflow.com/questions/463105/in-place-radix-sort
template <typename kmer_t, typename size_t_max>
template <typename array_element>
inline void CuttedSortedAC<kmer_t, size_t_max>::radix_sort(std::vector<array_element> &array, size_t_max begin, size_t_max end, size_t_max starting_base) {
    if (begin >= end) return;

    size_t_max AEnd = begin, TBegin = end;
    size_t_max i = begin;
    while (i < TBegin){
        char nucleotide = NucleotideAtIndex(array[i], K, starting_base);
        if (nucleotide == 'A') std::swap(array[i++], array[AEnd++]);
        else if (nucleotide == 'T') std::swap(array[i], array[--TBegin]);
        else ++i;
    }
    size_t_max CEnd = i = AEnd;
    while (i < TBegin){
        char nucleotide = NucleotideAtIndex(array[i], K, starting_base);
        if (nucleotide == 'C') std::swap(array[i], array[CEnd++]);
        ++i;
    }

    if (starting_base < K){
        radix_sort(array, begin, AEnd, (size_t_max)(starting_base + 1));  // A
        radix_sort(array, AEnd, CEnd, (size_t_max)(starting_base + 1));   // C
        radix_sort(array, CEnd, TBegin, (size_t_max)(starting_base + 1)); // G
        radix_sort(array, TBegin, end, (size_t_max)(starting_base + 1));  // T
    }
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::resort_and_shorten_failures(std::vector<std::pair<kmer_t, size_t_max>> &array, size_t_max depth) {
    size_t_max count = array.size();

    size_t_max skipped = 0;
    for (size_t_max i = 0; i < count; ++i){
        if (array[i].second == INVALID_NODE()) ++skipped;
        else array[i - skipped] = array[i];
    }

    count -= skipped;
    array.resize(count);
    std::vector<std::pair<kmer_t, size_t_max>> sorted_array;
    sorted_array.resize(count);

    size_t_max end_indexes[4] = {0, 0, 0, 0}; // Starts of: A, C, G, T
    for (size_t_max i = 0; i < count; ++i){
        uint8_t base_index = NucleotideIndexAtIndex(array[i], K, K - depth);
        ++end_indexes[base_index];

        array[i].first = BitSuffix(array[i].first, depth - 1);
    }

    size_t_max start_indexes[4] = {0, 0, 0, 0};
    for (uint8_t i = 1; i < 4; ++i){
        end_indexes[i] += end_indexes[i - 1];
        start_indexes[i] = end_indexes[i - 1];
    }

    for (size_t_max s = 0; s < count; ++s){
        uint8_t best_i = 4;
        for (uint8_t i = 0; i < 4; ++i){
            if (start_indexes[i] == end_indexes[i]) continue;
            if (best_i == 4 || array[start_indexes[best_i]].first > array[start_indexes[i]].first){
                best_i = i;
            }
        }
        sorted_array[s] = array[start_indexes[best_i]++];
    }

    array = std::move(sorted_array);
}

// Printing

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::print_stats(std::ostream &os)
{
    os << "Stats:" << std::endl;
    std::vector<size_t_max> depth_count(K + 1, 0);
    std::vector<size_t_max> real_depth_count(K + 1, 0);
    std::vector<size_t_max> real_depths(nodes.size(), 0);
    for (size_t_max i = nodes.size() - 1; i >= 0; --i){
        depth_count[nodes[i].depth]++;
        if (i != nodes.size() - 1) real_depths[i] = real_depths[nodes[i].parent] + 1;
        real_depth_count[real_depths[i]]++;
        if (i == 0) break;
    }
    os << "Depths:" << std::endl;
    for (size_t_max i = 0; i <= K; i++){
        os << i << ":\t" << depth_count[i] << std::endl;
    }
    os << "Real depths:" << std::endl;
    for (size_t_max i = 0; i <= K; i++){
        os << i << ":\t" << real_depth_count[i] << std::endl;
    }
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::print_sorted(std::ostream& os) {
    size_t_max i = 0;
    for (Node node : nodes){
        os << i++ << ":\t";
        node.print(os);
        os << ":\t";
        for (size_t_max c = 0; c < node.depth; ++c) os << NucleotideAtIndex(kMers[node.kmer_index], K, c);
        os << std::endl;
    }
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::print_topological(std::ostream& os) {
    print_topological(os, nodes.size() - 1, 0);
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::print_topological(std::ostream& os, size_t_max root, size_t_max depth) {
    Node node = nodes[root];
    
    for (size_t_max i = 0; i < depth; ++i) os << "|  ";
    os << "|->";
    
    os << root << ":\t";
    node.print(os);
    os << ":\t";
    for (size_t_max c = 0; c < node.depth; ++c) os << NucleotideAtIndex(kMers[node.kmer_index], K, c);

    if (node.parent != INVALID_NODE()){
        Node parent = nodes[node.parent];
        if (parent.child_range_begin > root || parent.child_range_end <= root) os << "xxx";
    }
    
    os << std::endl;

    if (root < N) return;
    
    for (size_t_max i = node.child_range_begin; i < node.child_range_end; ++i){
        if (i == INVALID_NODE()) break;
        print_topological(os, i, depth + 1);
    }
}

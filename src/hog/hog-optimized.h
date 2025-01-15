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
    
    Node() : kmer_index(INVALID_NODE()) {}; // Special invalid type of node

    inline bool is_invalid () const {
        return kmer_index == INVALID_NODE();
    };

    void print(std::ostream& os) const;
};

template <typename size_t_max>
inline void Node<size_t_max>::print(std::ostream &os) const {
    if (is_invalid()){
        os << "INV";
        return;
    }
    os << depth << '/';
    if (kmer_index == INVALID_NODE()) os << "INV"; else os << kmer_index;
    os << ",\tP: ";
    if (parent == INVALID_NODE()) os << "INV"; else os << parent;
    os << ",\tF: ";
    if (failure == INVALID_NODE()) os << "INV"; else os << failure;
    os << ",\tCH:[ ";
    if (child_range_begin == INVALID_NODE()) os << "INV"; else os << child_range_begin;
    os << " - ";
    if (child_range_end == INVALID_NODE()) os << "INV"; else os << child_range_end;
    os << " )";
}

template <typename kmer_t, typename size_t_max>
class HOGConstructer {
    std::vector<kmer_t>& kMers;
    std::vector<Node<size_t_max>> nodes;
    std::vector<std::pair<size_t_max, size_t_max>> leaf_intervals;
    
    size_t_max k; // kmer-length
    size_t_max n; // number of kmers, number of leaves
    
    static inline size_t_max INVALID_NODE() { return std::numeric_limits<size_t_max>::max(); };

    inline uint8_t get_edge_nucleotide(size_t_max child_index){
        return NucleotideIndexAtIndex(kMers[nodes[child_index].kmer_index], k, nodes[child_index].depth - 1);
    }
    void remove_redundant_edges(std::vector<bool>& keep);
    void recreate_BFS_indexes();
public:
    size_t VERBOSE = 0;
    HOGConstructer(std::vector<kmer_t>& kMers, size_t_max k) :
        kMers(kMers), k(k), n(kMers.size()) {};

    void construct_EHOG();
    //void convert_EHOG_to_HOG();
    void convert_to_leaf_ranges();
    // TODO convert child ranges to leaf ranges - or skip HOG and do it in construction? Is it possible?

    template<typename array_element>
    void sort(std::vector<array_element>& array, size_t_max starting_base = 0);
    //void radix_sort(std::vector<std::pair<kmer_t, size_t_max>>& array, size_t_max begin, size_t_max end, size_t_max starting_base = 0);
    void sort_by_base(std::vector<std::pair<kmer_t, size_t_max>> &array, size_t_max base);

    std::vector<size_t_max> compute_ordering(size_t_max new_run_score, size_t_max base_score = 1);

    void print_stats(std::ostream& os = std::cout);
    void print_sorted(std::ostream& os = std::cout);
    void print_topological(std::ostream& os = std::cout);
    void print_topological(std::ostream& os, size_t_max root, size_t_max depth);
};

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::construct_EHOG() {
    if (!nodes.empty()) throw std::invalid_argument("EHOG has already been constructed");

    sort(kMers);

    nodes.reserve(k * kMers.size()); // TODO decrease

    std::vector<Node<size_t_max>> current_nodes, last_nodes;
    last_nodes.reserve(n); // Might not be enought? Or?
    
    std::vector<std::pair<kmer_t, size_t_max>> failures; // suffix, last index
    
    // Add leaves
    current_nodes.resize(n);
    failures.resize(n);
    for (size_t_max i = 0; i < n; ++i){
        current_nodes[i] = Node(size_t_max(n - i - 1), k, INVALID_NODE(), INVALID_NODE());
        failures[i] = std::make_pair(kMers[n - i - 1], i);
    }
    sort(failures);
    
    // Add other nodes
    for (size_t_max depth = k; depth > 0; --depth){
        last_nodes = std::move(current_nodes);
        //std::cout << depth << ' ';
        //std::cout.flush();

        size_t_max saved_node_count = nodes.size();
        size_t_max last_node_count = last_nodes.size();
        size_t_max new_node_index = saved_node_count + last_node_count;
        size_t_max shift_from_removed = 0; // from last_nodes
        
        size_t_max failure_count = failures.size();
        for (size_t_max i = 0; i < failure_count; ++i) failures[i].first = BitSuffix(failures[i].first, depth - 1);
        //sort_by_base(failures, k - depth + 1);
        sort(failures);
        while (failures.back().second == INVALID_NODE()) failures.pop_back();
        size_t_max failure_index = failures.size() - 1;
        /*for (size_t_max i = failure_count; i > 0; --i){
            auto f = failures[i - 1];
            print_kmer(f.first, depth - 1, std::cout); std::cout << ' ' << f.second << std::endl;
        }*/

        std::vector<size_t_max> delayed_failure_indexes(last_node_count, INVALID_NODE());
        std::vector<size_t_max> updated_failures;

        for (size_t_max i = 0; i < last_node_count; ++i){
            kmer_t current_prefix = BitPrefix(kMers[last_nodes[i].kmer_index], k, depth - 1);
            while (failure_index > 0 && current_prefix < failures[failure_index].first) --failure_index;

            bool creating_new_node = (current_prefix == failures[failure_index].first);

            size_t_max j = i;
            while (j < last_node_count &&
                   current_prefix == BitPrefix(kMers[last_nodes[j].kmer_index], k, depth - 1)){
                //std::cout << i << ' ' << j << '+' << saved_node_count << std::endl;
                last_nodes[j].parent = new_node_index;

                
                if (delayed_failure_indexes[j] != INVALID_NODE()){
                    //std::cout << j << "->" << delayed_failure_indexes[j] << std::endl;
                    last_nodes[j].failure = delayed_failure_indexes[j]; // "Setting later"
                    if (creating_new_node){
                        updated_failures.push_back(saved_node_count + j);
                    }
                    else {
                        updated_failures.push_back(new_node_index + j - i);
                    }
                }

                ++j;
                //std::cout << j << std::endl;
                //std::cout << kMers[last_nodes[i].kmer_index] << ' ' << kMers[last_nodes[j].kmer_index] << ' ' << max_diff << std::endl;
            }
            
            //std::cout << i << ' ' << j << std::endl;
            if (creating_new_node){
                //std::cout << "+" << std::endl;
                current_nodes.push_back(Node(last_nodes[i].kmer_index,
                                             size_t_max(depth - 1),
                                             size_t_max(i - shift_from_removed + saved_node_count),
                                             size_t_max(j - shift_from_removed + saved_node_count)));

                bool first_similar_failure = true;
                while (current_prefix == failures[failure_index].first){
                    size_t_max last_failure_index = failures[failure_index].second;

                    if (last_failure_index < saved_node_count){
                        //std::cout << "old" << std::endl;
                        nodes[last_failure_index].failure = new_node_index;
                        updated_failures.push_back(last_failure_index);
                    }
                    else if (last_failure_index - saved_node_count < j){
                        //std::cout << "last" << std::endl;
                        last_nodes[last_failure_index - saved_node_count].failure = new_node_index;
                        updated_failures.push_back(last_failure_index);
                    }
                    else {
                        //std::cout << "delay" << std::endl;
                        delayed_failure_indexes[last_failure_index - saved_node_count] = new_node_index; // Set later
                    }
                    
                    if (first_similar_failure){
                        first_similar_failure = false;
                        failures[failure_index].second = new_node_index;
                    }
                    else {
                        failures[failure_index] = std::make_pair(std::numeric_limits<kmer_t>::max(), INVALID_NODE()); // Will be last, lost
                    }

                    if (failure_index == 0) break;
                    --failure_index;
                }
                
                ++new_node_index;
            }
            else {
                //std::cout << i + saved_node_count << ' ' << j << new_node_index << std::endl;
                if (depth == k){
                    for (size_t_max copy_index = i; copy_index < j; ++copy_index){
                        Node copy_node = last_nodes[copy_index];
                        current_nodes.push_back(Node(copy_node.kmer_index, copy_node.depth, copy_index, size_t_max(copy_index + 1)));
                        //current_nodes.back().failure = copy_node.failure;
                        failures[copy_index].second = new_node_index;
                        ++new_node_index;
                    }
                }
                else {
                    for (size_t_max swap_index = i; swap_index < j; ++swap_index){
                        current_nodes.push_back(Node<size_t_max>());
                        std::swap(current_nodes.back(), last_nodes[swap_index]);
                    }
                    shift_from_removed += j - i;
                    new_node_index += j - i;
                }
            }

            i = j - 1;
        }

        for (size_t_max update_index : updated_failures){
            if (update_index < saved_node_count){
                nodes[update_index].failure -= shift_from_removed;
            }
            else if (update_index < saved_node_count + last_node_count){
                last_nodes[update_index - saved_node_count].failure -= shift_from_removed;
            }
            else {
                current_nodes[update_index - saved_node_count - last_node_count].failure -= shift_from_removed;
            }
        }
        
        for (Node node: last_nodes){
            if (!node.is_invalid()){
                node.parent -= shift_from_removed;
                
                nodes.push_back(node);
            }
        }

        for (size_t_max i = 0; i < n; ++i){
            if (failures[i].second != INVALID_NODE() && failures[i].second >= saved_node_count + last_node_count){
                failures[i].second -= shift_from_removed;
            }
        }

        /*size_t_max i = 0;
        std::cout << "nodes:" << std::endl;
        for (Node node : nodes){
            std::cout << i++ << ":\t";
            node.print(std::cout);
            std::cout << ":\t";
            for (size_t_max c = 0; c < node.depth; ++c) std::cout << to_upper(NucleotideAtIndex(kMers[node.kmer_index], k, c));
            std::cout << std::endl;
        }*/
        /*std::cout << "last:" << std::endl;
        for (Node node : last_nodes){
            std::cout << i++ << ":\t";
            node.print(std::cout);
            std::cout << ":\t";
            for (size_t_max c = 0; c < node.depth; ++c) std::cout << to_upper(NucleotideAtIndex(kMers[node.kmer_index], k, c));
            std::cout << std::endl;
        }*/
        /*std::cout << "current:" << std::endl;
        for (Node node : current_nodes){
            std::cout << i++ << ":\t";
            node.print(std::cout);
            std::cout << ":\t";
            for (size_t_max c = 0; c < node.depth; ++c) std::cout << to_upper(NucleotideAtIndex(kMers[node.kmer_index], k, c));
            std::cout << std::endl;
        }
        std::cout << std::endl;*/
    }

    for (Node node: current_nodes) nodes.push_back(node);

    nodes.shrink_to_fit();
    std::cout << std::endl;

    size_t_max node_count = nodes.size();
    for (size_t_max i = n; i < node_count; ++i){
        //if (nodes[i].depth == k) continue;
        for (size_t_max c = nodes[i].child_range_begin; c < nodes[i].child_range_end; ++c){
            nodes[c].parent = i;
        }
    }
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::convert_to_leaf_ranges() {
    size_t_max node_count = nodes.size();
    for (size_t_max i = 0; i < node_count; ++i){
        if (nodes[i].depth == k){
            nodes[i].child_range_begin = i;
            nodes[i].child_range_end = i + 1;
        }
        else {
            nodes[i].child_range_begin = nodes[nodes[i].child_range_begin].child_range_begin;
            nodes[i].child_range_end = nodes[nodes[i].child_range_end - 1].child_range_end;
        }
    }
}

template <typename kmer_t, typename size_t_max>
template <typename array_element>
inline void HOGConstructer<kmer_t, size_t_max>::sort(std::vector<array_element> &array, size_t_max starting_base) {
    /*if (array.size() > size_t(1 << (k - starting_base))){
        radix_sort(array, 0, array.size(), starting_base);
    }
    else {
        std::sort(array.begin(), kMers.end());
    }*/
    std::sort(array.begin(), array.end());
}

/*
//in-place radix sort as in https://stackoverflow.com/questions/463105/in-place-radix-sort
template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::radix_sort(std::vector<std::pair<kmer_t, size_t_max>> &array, size_t_max begin, size_t_max end, size_t_max starting_base) {
    if (begin >= end) return;

    size_t_max AEnd = begin, TBegin = end;
    size_t_max i = begin;
    while (i < TBegin){
        char nucleotide = NucleotideAtIndex(array[i].first, k, starting_base);
        if (nucleotide == 'A') std::swap(array[i++], array[AEnd++]);
        else if (nucleotide == 'T') std::swap(array[i], array[--TBegin]);
        else ++i;
    }
    size_t_max CEnd = i = AEnd;
    while (i < TBegin){
        char nucleotide = NucleotideAtIndex(array[i].first, k, starting_base);
        if (nucleotide == 'C') std::swap(array[i], array[CEnd++]);
        ++i;
    }

    if (starting_base < k){
        radix_sort(array, begin, AEnd, (size_t_max)(starting_base + 1));  // A
        radix_sort(array, AEnd, CEnd, (size_t_max)(starting_base + 1));   // C
        radix_sort(array, CEnd, TBegin, (size_t_max)(starting_base + 1)); // G
        radix_sort(array, TBegin, end, (size_t_max)(starting_base + 1));  // T
    }
}*/

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::sort_by_base(std::vector<std::pair<kmer_t, size_t_max>> &array, size_t_max starting_base) {
    size_t_max start_indexes[5] = {0, 0, 0, 0, 0}; // Starts of: A, C, G, T, INV
    for (auto p : array){
        if (p.second == INVALID_NODE()) continue; // Not counting invalid nodes
        uint8_t base_index = NucleotideIndexAtIndex(p.first, k, starting_base);
        ++start_indexes[base_index + 1];
    }

    size_t_max end_indexes[5] = {0, 0, 0, 0, 0};
    for (uint8_t i = 0; i < 4; ++i){
        start_indexes[i + 1] += start_indexes[i];
        end_indexes[i] = start_indexes[i + 1];
    }

    for (uint8_t i = 0; i < 4; ++i){
        for (size_t_max s = start_indexes[i]; s < end_indexes[i]; ++s){
            uint8_t base_index;
            if (array[s].second == INVALID_NODE()) base_index = 4;
            else base_index = NucleotideIndexAtIndex(array[s].first, k, starting_base);

            if (base_index != i) std::swap(array[s--], array[start_indexes[base_index]++]);
        }
    }
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::print_stats(std::ostream &os)
{
    os << "Stats:" << std::endl;
    std::vector<size_t_max> depth_count(k + 1, 0);
    std::vector<size_t_max> real_depth_count(k + 1, 0);
    std::vector<size_t_max> real_depths(nodes.size(), 0);
    for (size_t_max i = nodes.size() - 1; i >= 0; --i){
        depth_count[nodes[i].depth]++;
        if (i != nodes.size() - 1) real_depths[i] = real_depths[nodes[i].parent] + 1;
        real_depth_count[real_depths[i]]++;
        if (i == 0) break;
    }
    os << "Depths:" << std::endl;
    for (size_t_max i=0; i <= k; i++){
        os << i << ":\t" << depth_count[i] << std::endl;
    }
    os << "Real depths:" << std::endl;
    for (size_t_max i=0; i <= k; i++){
        os << i << ":\t" << real_depth_count[i] << std::endl;
    }
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::print_sorted(std::ostream& os) {
    size_t_max i = 0;
    for (Node node : nodes){
        os << i++ << ":\t";
        node.print(os);
        os << ":\t";
        for (size_t_max c = 0; c < node.depth; ++c) os << to_upper(NucleotideAtIndex(kMers[node.kmer_index], k, c));
        os << std::endl;
    }
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::print_topological(std::ostream& os) {
    print_topological(os, nodes.size() - 1, 0);
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::print_topological(std::ostream& os, size_t_max root, size_t_max depth) {
    Node node = nodes[root];
    
    for (size_t_max i = 0; i < depth; ++i) os << "|  ";
    os << "|->";
    
    os << root << ":\t";
    node.print(os);
    os << ":\t";
    for (size_t_max c = 0; c < node.depth; ++c) os << to_upper(NucleotideAtIndex(kMers[node.kmer_index], k, c));

    if (node.parent != INVALID_NODE()){
        Node parent = nodes[node.parent];
        if (parent.child_range_begin > root || parent.child_range_end <= root) os << "xxx";
    }
    
    os << std::endl;

    if (root < n) return;
    
    for (size_t_max i = node.child_range_begin; i < node.child_range_end; ++i){
        if (i == INVALID_NODE()) break;
        print_topological(os, i, depth + 1);
    }
}

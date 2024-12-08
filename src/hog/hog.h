#pragma once
/*
#define kmer_index unsigned int
#define node_index unsigned int
#define kmar_max_len unsigned int
*/
#include <vector>
#include <queue>

#include "../kmers.h"

//in-place radix sort as in https://stackoverflow.com/questions/463105/in-place-radix-sort
template <typename kmer_t>
void radix_sort(std::vector<kmer_t>& kMers, size_t k){
    if (kMers.empty()) return;

    size_t n = kMers.size();
    for (size_t pos = 0; pos < k; ++pos){
        size_t ACEnd = 0, TBegin = n;
        size_t i = 0;
        while (i < TBegin){
            char nucleotide = NucleotideAtIndex(kMers[i], k, pos);
            if (nucleotide == 'A') std::swap(kMers[i++], kMers[ACEnd++]);
            else if (nucleotide == 'T') std::swap(kMers[i], kMers[--TBegin]);
            else ++i;
        }
        i = ACEnd;
        while (i < TBegin){
            char nucleotide = NucleotideAtIndex(kMers[i], k, pos);
            if (nucleotide == 'C') std::swap(kMers[i], kMers[ACEnd++]);
            ++i;
        }
    }

    /*size_t n = kMers.size();
    std::vector<kmer_t>[4] buckets;
    for (size_t i = 0; i < k; ++i){
        buckets[0].clear(); buckets[1].clear(); buckets[2].clear(); buckets[3].clear();
        for (size_t j = n - 1; j >= 0; --j){
            size_t bucket = NucleotideIndexAtIndex(kMers[j], k, i);
            if (i % 2) bucket = 3 - bucket;
            buckets[bucket].push_back(kMers[j]);
        }
        size_t j = 0;
        for (auto kmer : buckets[0]) kMers[j++] = kmer;
        for (auto kmer : buckets[1]) kMers[j++] = kmer;
        for (auto kmer : buckets[2]) kMers[j++] = kmer;
        for (auto kmer : buckets[3]) kMers[j++] = kmer;
    }*/
}

template<typename size_t_max>
struct Node {
    size_t_max kmer_index;
    size_t_max depth;
    size_t_max index; // Might not need it in the end
    size_t_max parent;
    size_t_max failure;
    bool has_basic_child[4]; // Only used when constructing AC
    std::vector<size_t_max> children;

    Node(size_t_max kmer_index, size_t_max depth, size_t_max index, size_t_max parent) :
        kmer_index(kmer_index), depth(depth), index(index), parent(parent),
        failure(0), has_child({false, false, false, false}, children(4)) {};
};

template <typename kmer_t, typename size_t_max>
struct HOGConstructer {
    std::vector<kmer_t>& kMers;
    std::vector<Node<size_t_max>> nodes;
    std::vector<size_t_max> leaves;
    std::vector<size_t_max> BFS_indexes;
    size_t_max k;
    size_t_max n;

    //constexpr size_t_max INVALID_NODE = std::numeric_limits<size_t_max>::max();
    //static bool is_node_valid(size_t_max node_index) { return node_index != INVALID_NODE};

    short get_edge_nucleotide(size_t_max child_index){
        return NucleotideIndexAtIndex(kMers[nodes[child_index].kmer_index], k, nodes[child_index].depth);
    }

    HOGConstructer(size_t_max k) :
        k(k), n(0) {};
    HOGConstructer(size_t_max k, std::vector<kmer_t>& kMers) :
        k(k), kMers(kMers), n(kMers.size()) {};
    void construct_AC();
    void convert_AC_to_EHOG();
    void convert_EHOG_to_HOG();
};

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::construct_AC() {
    if (!nodes.empty()) throw std::invalid_argument("AC has already been constructed");

    // Add nodes and forward edges
    nodes.resize(n); // TODO better guess
    nodes.emplace_back<Node<size_t_max>>(0, 0, 0, 0);
    size_t_max node_count = 1;

    leaves.resize(n);

    for (size_t_max kmer_index = 0; kmer_index < n; ++kmer_index){
        size_t_max node_index = 0;
        for (size_t_max pos = 0; pos < k; ++pos){
            short nucleotide_index = NucleotideIndexAtIndex(kMers[kmer_index], k, pos);
            if (nodes[node_index].has_basic_child[nucleotide_index]){
                // Position of new node = node count
                nodes.emplace_back<Node<size_t_max>>(kmer_index, pos + 1, node_count, node_index);
                // Add new child to parent
                nodes[node_index].children.push_back(node_count);
                nodes[node_index].has_basic_child[get_edge_nucleotide(node_count)] = true;
                // Increase node count
                node_count++;

                // Next step
                node_index = node_count - 1;
            }
            else node_index = nodes[node_index].children[nucleotide_index];
        }
        leaves.push_back(node_index);
    }

    // Add fake goto links to root
    for (short i = 0; i < 4; ++i){
        if (!nodes[0].has_basic_child[i]){
            nodes[0].has_basic_child[i];
        }
    }

    BFS_indexes.resize(node_count);

    // Add failure links
    std::queue<size_t_max> queue;
    for (size_t_max child_index : nodes[0].children) queue.push(child_index);

    while (!queue.empty()){
        size_t_max index = queue.front();
        queue.pop();

        // Add FL to all children of nodes[index]
        for (size_t_max child_index : nodes[index].children){
            short goto_char = get_edge_nucleotide(child_index);

            // Search for first node with valid goto edge
            size_t_max failure_index = nodes[index].failure;
            while(!nodes[failure_index].has_basic_child(goto_char))
                failure_index = nodes[failure_index].failure;

            // Find its child with right edge
            bool found = false;
            for (size_t_max failure_child : nodes[failure_index].children){
                if (goto_char == get_edge_nucleotide(failure_child)){
                    nodes[child_index].failure = failure_child;
                    found = true;
                    break;
                }
            }
            // Only possible in the root
            if (!found) nodes[child_index].failure = 0;

            BFS_indexes.push_back(index);

            queue.push(child_index);
        }
    }

    nodes.shrink_to_fit();
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::convert_AC_to_EHOG() {
    if (nodes.empty) throw std::invalid_argument("Cannot convert empty AC to EHOG");

    // Mark nodes to be keeped
    size_t_max node_count = nodes.size();
    std::vector<bool> keep(node_count, false);
    for (size_t_max leaf_index : leaves){
        keep[leaf_index] = true;
        size_t_max failure_index = nodes[leaf_index].failure;
        while (failure_index != 0){
            keep[failure_index] = true;
            failure_index = nodes[failure_index].failure;
        }
    }

    // Contract edges from and to non-keeped nodes
    // Ignore failure links (those we care about are already good)
    ;
    // Replace up-facing edges in BFS order
    for (size_t_max i : BFS_indexes){
        if (keep[i]){
            // Add this node as child of its new parent
            nodes[nodes[i].parent].children.push_back(i);

            // Remove its old children
            nodes[i].children.clear();
            // (Also removes loop from root)
        }
        else {
            // Tell all its children it is no more their parent
            for (size_t_max child_index : nodes[i].children) nodes[child_index].parent = nodes[i].parent;
        }   
    }

    // Remove non-keeped nodes
    // Move keeped nodes to the beginning
    std::vector<size_t_max> new_indexes(node_count, 0);
    size_t_max keeped = 0;
    for (size_t_max i = 0; i < node_count; ++i){
        if (keep[i]){
            std::swap(nodes[i], nodes[keeped++]);
            new_indexes[i] = keeped - 1;
        }
    }
    // Move all the edges to point to correct locations
    for (size_t_max i = 0; i < keeped; ++i){
        nodes[i].index = new_indexes[i];
        nodes[i].parent = new_indexes[nodes[i].parent];
        nodes[i].failure = new_indexes[nodes[i].failure];
        
        size_t_max child_count = nodes[i].children.size();
        for (size_t_max j = 0; j < child_count; ++j)
            nodes[i].children[j] = new_indexes[nodes[i].children[j]];
    }
    
    // Get rid of non-keeped nodes
    nodes.resize(keeped);
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::convert_EHOG_to_HOG() {
    if (nodes.empty) throw std::invalid_argument("Cannot convert empty EHOG to HOG");

    size_t_max node_count = nodes.size();
    std::vector<bool> keep(node_count, false);
    keep[0] = true;
}

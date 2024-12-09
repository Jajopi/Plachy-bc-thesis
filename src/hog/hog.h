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
    //size_t_max index; // Might not need it in the end
    size_t_max parent;
    size_t_max failure;
    std::vector<size_t_max> children;

    Node(size_t_max kmer_index, size_t_max depth, size_t_max index, size_t_max parent) :
        kmer_index(kmer_index), depth(depth), index(index), parent(parent),
        failure(0), has_child({false, false, false, false}, children(4)) {};
};

template <typename kmer_t, typename size_t_max>
class HOGConstructer {
    std::vector<kmer_t>& kMers;
    std::vector<Node<size_t_max>> nodes;
    std::vector<size_t_max> leaves;
    std::vector<size_t_max> BFS_indexes;
    size_t_max k;
    size_t_max n;

    inline uint8_t get_edge_nucleotide(size_t_max child_index){
        return NucleotideIndexAtIndex(kMers[nodes[child_index].kmer_index], k, nodes[child_index].depth);
    }
    void remove_redundant_edges(std::vector<bool>& keep);
public:
    /*HOGConstructer(size_t_max k) :
        k(k), n(0) {};*/
    HOGConstructer(std::vector<kmer_t>& kMers, size_t_max k) :
        k(k), kMers(kMers), n(kMers.size()) {};

    void construct_AC();
    void convert_AC_to_EHOG();
    void convert_EHOG_to_HOG();
    
    inline void create() {
        construct_AC();
        convert_AC_to_EHOG();
        convert_EHOG_to_HOG();
    }
};


// Remove non-keeped nodes
template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::remove_redundant_edges(std::vector<bool>& keep)
{
    size_t_max node_count = nodes.size();
    if (node_count != keep.size()) throw std::invalid_argument("Invalid removing mask passed as argument");

    // Move keeped nodes to the beginning
    std::vector<size_t_max> new_indexes(node_count, 0);
    size_t_max keeped = 0;
    for (size_t_max i = 0; i < node_count; ++i){
        if (keep[i]){
            std::swap(nodes[i], nodes[keeped++]);
            new_indexes[i] = keeped - 1;
        }
    }

    // Set all the edges to point to correct locations
    for (size_t_max i = 0; i < keeped; ++i){
        //nodes[i].index = new_indexes[i];
        nodes[i].parent = new_indexes[nodes[i].parent];
        nodes[i].failure = new_indexes[nodes[i].failure];
        
        size_t_max child_count = nodes[i].children.size();
        for (size_t_max j = 0; j < child_count; ++j)
            nodes[i].children[j] = new_indexes[nodes[i].children[j]];
    }

    // Set correct pointers from leaves and BFS_indexes
    for (size_t_max i = 0; i < n; ++i) leaves[i] = new_indexes[leaves[i]];
    for (size_t_max i = 0; i < keeped; ++i) BFS_indexes[i] = new_indexes[BFS_indexes[i]];
    
    // Get rid of non-keeped nodes
    nodes.resize(keeped);
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::construct_AC()
{
    if (!nodes.empty()) throw std::invalid_argument("AC has already been constructed");

    // Add nodes and forward edges
    nodes.resize(n); // TODO better guess
    nodes.emplace_back<Node<size_t_max>>(0, 0, 0, 0);
    size_t_max node_count = 1;

    leaves.resize(n);

    std::vector<bool> has_basic_child(4 * node_count, false);

    for (size_t_max kmer_index = 0; kmer_index < n; ++kmer_index){
        size_t_max node_index = 0;
        for (size_t_max pos = 0; pos < k; ++pos){
            uint8_t nucleotide_index = NucleotideIndexAtIndex(kMers[kmer_index], k, pos);
            if (has_basic_child[nodes[node_index] * 4 + nucleotide_index]){
                // Position of new node = node count
                nodes.emplace_back<Node<size_t_max>>(kmer_index, pos + 1, node_count, node_index);
                // Add new child to parent
                nodes[node_index].children.push_back(node_count);
                has_basic_child[nodes[node_index] * 4 + get_edge_nucleotide(node_count)] = true;
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
    for (uint8_t i = 0; i < 4; ++i) has_basic_child[i] = true;

    BFS_indexes.resize(node_count);

    // Add failure links
    std::queue<size_t_max> queue;
    for (size_t_max child_index : nodes[0].children) queue.push(child_index);

    while (!queue.empty()){
        size_t_max index = queue.front();
        queue.pop();

        // Add FL to all children of nodes[index]
        for (size_t_max child_index : nodes[index].children){
            uint8_t goto_char = get_edge_nucleotide(child_index);

            // Search for first node with valid goto edge
            size_t_max failure_index = nodes[index].failure;
            while(!has_basic_child[nodes[failure_index] * 4 + goto_char])
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

    remove_redundant_edges(keep);
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::convert_EHOG_to_HOG() {
    if (nodes.empty) throw std::invalid_argument("Cannot convert empty EHOG to HOG");

    // Initialize resulting array
    size_t_max node_count = nodes.size();
    std::vector<bool> keep(node_count, false);
    keep[0] = true;

    // Initialize helper arrays
    std::vector<bool> black(node_count, 0);
    std::vector<uint8_t> child_count(node_count, 0);
    std::vector<size_t_max> favorite_proper_ancestor(node_count, 0);
    std::vector<size_t_max> favorite_descendant(node_count, 0);
    std::vector<size_t_max> modified_vertices;
    // Child counts
    for (size_t_max i = 0; i < node_count; ++i) child_count[i] = nodes[i].children.size();
    // FavPAnc
    for (size_t_max i : BFS_indexes){
        size_t_max proper_ancestor = nodes[i].parent;
        if (child_count[proper_ancestor] == 1)
            proper_ancestor = favorite_proper_ancestor[proper_ancestor];
        
        favorite_proper_ancestor[i] = proper_ancestor;
    }
    // FavDesc
    for (size_t_max i = node_count - 1; i >= 0; --i){
        size_t_max descendant = BFS_indexes[i];
        if (child_count[descendant] == 1) descendant = favorite_descendant[descendant];
        
        favorite_descendant[BFS_indexes[i]] = descendant;
    }
    // Mark all leaves white
    for (size_t_max i = 0; i < n; ++i) black[leaves[i]] = true;

    for (size_t_max x : leaves){
        keep[x] = true;
        size_t_max v = nodes[x].failure;
        modified_vertices.clear();
        while (v != 0){
            size_t_max u = favorite_descendant[v];
            if (child_count[u] != 0){
                child_count[u] = 0;
                keep[x] = true;

                while (u != 0){
                    modified_vertices.push_back(u);
                    if (u != favorite_descendant[v]){
                        --child_count[u];
                        if (child_count[u] != 0) break;
                    }

                    u = favorite_proper_ancestor[u];
                }
            }

            v = nodes[v].failure;
        }
        if (!modified_vertices.empty()){
            for (size_t_max v : modified_vertices) ++child_count[v];// = nodes[v].children.size();
            child_count[modified_vertices[0]] = nodes[modified_vertices[0]].children.size();
        }
    }

    remove_redundant_edges(keep);
}

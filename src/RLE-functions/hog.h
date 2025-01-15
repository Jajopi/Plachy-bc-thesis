#pragma once
/*
#define kmer_index unsigned int
#define node_index unsigned int
#define kmar_max_len unsigned int
*/
#include <vector>
#include <queue>
#include <algorithm>

#include "../kmers.h"

template<typename size_t_max>
struct Node {
    size_t_max kmer_index;
    size_t_max depth;
    //size_t_max index; // Might not need it in the end
    size_t_max parent;
    size_t_max failure;
    std::vector<size_t_max> children;

    Node(size_t_max kmer_index, size_t_max depth, /*size_t_max index,*/ size_t_max parent) :
        kmer_index(kmer_index), /*index(index),*/ depth(depth), parent(parent),
        failure(0) {};
    Node() { throw std::invalid_argument("Creating empty node"); };

    void print(std::ostream& os) const;
};

template <typename size_t_max>
inline void Node<size_t_max>::print(std::ostream &os) const {
    os << depth << '/' << kmer_index << ", P: " << parent << ", F: " << failure << ", CH: [ ";
    for (auto child : children) os << child << ' ';
    os << ']';
}

template <typename kmer_t, typename size_t_max>
class HOGConstructer {
    std::vector<kmer_t>& kMers;
    std::vector<Node<size_t_max>> nodes;
    std::vector<size_t_max> leaves; // TODO remove (not needed: BFS_indexes[-n:]) ... Or maybe not
    std::vector<size_t_max> BFS_indexes; // TODO sort nodes in this order effectively and get rid of this
    size_t_max k;
    size_t_max n;

    inline uint8_t get_edge_nucleotide(size_t_max child_index){
        return NucleotideIndexAtIndex(kMers[nodes[child_index].kmer_index], k, nodes[child_index].depth - 1);
    }
    void remove_redundant_edges(std::vector<bool>& keep);
    void recreate_BFS_indexes();
public:
    size_t VERBOSE = 0;
    /*HOGConstructer(size_t_max k) :
        k(k), n(0) {};*/
    HOGConstructer(std::vector<kmer_t>& kMers, size_t_max k) :
        kMers(kMers), k(k), n(kMers.size()) {};

    void construct_AC();
    void convert_AC_to_EHOG();
    void convert_EHOG_to_HOG();
    
    inline void create() {
        if (VERBOSE > 0) std::cout << "Constructing AC..." << std::endl;
        construct_AC();
        if (VERBOSE > 1) print();
        if (VERBOSE > 0) std::cout << "Converting to EHOG..." << std::endl;
        convert_AC_to_EHOG();
        if (VERBOSE > 1 )print();
        if (VERBOSE > 0) std::cout << "Converting to HOG..." << std::endl;
        convert_EHOG_to_HOG();
        if (VERBOSE > 1) print();
        if (VERBOSE > 0) std::cout << "Done: " << nodes.size() << " nodes" << std::endl;
    }

    std::vector<size_t_max> compute_ordering(size_t_max new_run_score, size_t_max base_score = 1);

    void print_stats(std::ostream& os = std::cout);
    void print(std::ostream& os = std::cout);
};


// Remove non-keeped nodes
template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::remove_redundant_edges(std::vector<bool>& keep)
{
    size_t_max node_count = nodes.size();
    if (node_count != keep.size()) throw std::invalid_argument("Invalid removing mask passed as argument");

    if (VERBOSE > 0) std::cout << "Removing redundant...";

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

    // Move keeped nodes to the beginning
    std::vector<size_t_max> new_indexes(node_count, 0);
    size_t_max keeped = 0;
    for (size_t_max i = 0; i < node_count; ++i){
        if (keep[i]){
            std::swap(nodes[i], nodes[keeped++]);
            new_indexes[i] = keeped - 1;
        }
    }

    // Update all the edges to point to correct locations
    for (size_t_max i = 0; i < keeped; ++i){
        //nodes[i].index = new_indexes[i];
        nodes[i].parent = new_indexes[nodes[i].parent];
        nodes[i].failure = new_indexes[nodes[i].failure];
        for (size_t_max& child : nodes[i].children) child = new_indexes[child];
    }

    // Update pointers from leaves
    for (size_t_max& leaf : leaves) leaf = new_indexes[leaf];
    
    if (VERBOSE > 0) std::cout << " keeped " << keeped << '/' << node_count << std::endl;
    // Get rid of non-keeped nodes
    nodes.resize(keeped);
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::recreate_BFS_indexes(){
    if (VERBOSE > 0) std::cout << "Recreating BFS indexes...";

    BFS_indexes.clear();
    BFS_indexes.push_back(0);

    std::queue<size_t_max> queue;
    for (size_t_max child_index : nodes[0].children) queue.push(child_index);

    while (!queue.empty()){
        size_t_max index = queue.front();
        queue.pop();

        BFS_indexes.push_back(index);

        for (size_t_max child_index : nodes[index].children){
            queue.push(child_index);
        }
    }

    if (VERBOSE > 0) std::cout << " done" << std::endl;
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::construct_AC()
{
    if (!nodes.empty()) throw std::invalid_argument("AC has already been constructed");

    // Add nodes and forward edges
    //nodes.resize(n); // TODO better guess
    nodes.emplace_back(0, 0, /*0,*/ 0);
    size_t_max node_count = 1;

    leaves.clear();

    std::vector<bool> has_basic_child(4, false);

    for (size_t_max kmer_index = 0; kmer_index < n; ++kmer_index){
        size_t_max node_index = 0;
        for (size_t_max depth = 0; depth < k; ++depth){
            uint8_t nucleotide_index = NucleotideIndexAtIndex(kMers[kmer_index], k, depth);

            if (has_basic_child[node_index * 4 + nucleotide_index]){
                for (size_t_max child : nodes[node_index].children){
                    if (get_edge_nucleotide(child) == nucleotide_index){
                        node_index = child;
                        break;
                    }
                }
            }
            else {
                // Index of new node = node count
                nodes.emplace_back(kmer_index, depth + 1, /*node_count,*/ node_index);
                for (auto ii = 0; ii < 4; ++ii) has_basic_child.push_back(false);
                // Add new child to parent
                nodes[node_index].children.push_back(node_count);
                has_basic_child[node_index * 4 + get_edge_nucleotide(node_count)] = true;
                // Increase node count
                node_count++;

                // Next step
                node_index = node_count - 1;
            }
        }
        leaves.push_back(node_index);
    }

    // Add fake goto links to root
    for (uint8_t i = 0; i < 4; ++i) has_basic_child[i] = true;

    BFS_indexes.push_back(0);

    // Add failure links and create BFS ordering
    std::queue<size_t_max> queue;
    for (size_t_max child_index : nodes[0].children) queue.push(child_index);

    while (!queue.empty()){
        size_t_max index = queue.front();
        queue.pop();

        BFS_indexes.push_back(index);

        // Add FL to all children of nodes[index]
        for (size_t_max child_index : nodes[index].children){
            uint8_t goto_char = get_edge_nucleotide(child_index);

            // Search for first node with valid goto edge
            size_t_max failure_index = nodes[index].failure;
            while(!has_basic_child[failure_index * 4 + goto_char])
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

            queue.push(child_index);
        }
    }
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::convert_AC_to_EHOG() {
    if (nodes.empty()) throw std::invalid_argument("Cannot convert empty AC to EHOG");

    // Mark nodes to be keeped
    size_t_max node_count = nodes.size();
    std::vector<bool> keep(node_count, false);
    std::vector<bool> moved(node_count, false);

    keep[0] = true;
    for (size_t_max leaf_index : leaves){
        keep[leaf_index] = true;
        size_t_max failure_index = nodes[leaf_index].failure;
        while (failure_index != 0){
            keep[failure_index] = true;
            failure_index = nodes[failure_index].failure;
        }
    }

    remove_redundant_edges(keep);
    recreate_BFS_indexes();
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::convert_EHOG_to_HOG() {
    if (nodes.empty()) throw std::invalid_argument("Cannot convert empty EHOG to HOG");

    // Initialize resulting array
    size_t_max node_count = nodes.size();
    std::vector<bool> keep(node_count, false);

    // Initialize helper arrays
    std::vector<bool> white(node_count, false);
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
    for (size_t_max i = node_count; i > 0; --i){
        size_t_max descendant = BFS_indexes[i - 1];
        if (child_count[descendant] == 1)
            descendant = favorite_descendant[nodes[descendant].children[0]];
        favorite_descendant[BFS_indexes[i - 1]] = descendant;
    }
    // Mark all leaves white
    for (size_t_max i = 0; i < n; ++i) white[leaves[i]] = true;

    keep[0] = true;
    for (size_t_max x : leaves){
        keep[x] = true;
        size_t_max v = nodes[x].failure;
        modified_vertices.clear();
        while (v != 0){
            size_t_max u = favorite_descendant[v];
            if (child_count[u] != 0){
                child_count[u] = 0;
                keep[v] = true;

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
    recreate_BFS_indexes();
}

template <typename kmer_t, typename size_t_max>
inline std::vector<size_t_max> HOGConstructer<kmer_t, size_t_max>::compute_ordering(size_t_max new_run_score, size_t_max base_score) {
    size_t_max node_count = nodes.size();
    
    std::vector<size_t_max> indexes;
    std::vector<bool> was_used(node_count, false);
    size_t_max used = 0;

    std::vector<std::pair<size_t_max, size_t_max>> heap;
    
    size_t_max last_index = 0;
    indexes.push_back(last_index);
    was_used[last_index] = true;
    ++used;

    while (used < n){
        heap.clear();
        size_t_max failure_index = nodes[last_index].failure;
        heap.push_back(std::make_pair(0, failure_index));
        bool found_next = false;

        while (!heap.empty()){
            size_t_max current = heap[0].second;
            size_t_max score = heap[0].first;
            std::pop_heap(heap); heap.pop_back();
            
            size_t_max current_depth = nodes[current].depth;
            failure_index = nodes[current].failure;
            size_t_max failure_depth = nodes[failure_index].depth;
            //bool failure_disrupts_run = nodes[next].depth == k - 1;
            heap.push_back(std::make_pair(score + base_score * (current_depth - failure_depth) +
                                          current_depth == k - 1 ? new_run_score : 0),
                                          failure_index); std::push_heap(heap);

            for (size_t_max child_index : nodes[current].children){
                size_t_max child_depth = nodes[child_index].depth;
                if (child_depth == k){
                    if (!was_used[child_index]){
                        last_index = child_index;
                        found_next = true;
                        break;
                    }

                    failure_index = nodes[child_index].failure;
                    failure_depth = nodes[failure_index].depth;
                    heap.push_back(std::make_pair(score + base_score * (k - failure_depth) +
                                                  failure_depth < k - 1 ? new_run_score : 0),
                                                  failure_index); std::push_heap(heap);
                }
                else {
                    heap.push_back(std::make_pair(score, child_index));
                }
            }
            
        }
    }

    was_used[last] = true;
    ++used;

    /*std::vector<size_t_max> successors(node_count, 0);

    std::vector<bool> has_in_edge(n, false), has_out_edge(n, false);
    std::vector<std::vector<std::pair<int, size_t_max>>> outs(node_count);
    std::vector<std::vector<size_t_max>> ins(node_count);

    for (size_t_max i = node_count - 1; i >= 0; ++i){
        size_t_max node_index = BFS_indexes[i];

        size_t_max parent_index = nodes[node_index].parent;
        //int depth_difference_parent = nodes[node_index].depth - nodes[parent_index].depth;
        size_t_max failure_index = nodes[node_index].failure;
        //int depth_difference_failure = nodes[node_index].depth - nodes[failure_index].depth;

        if (nodes[node_index].depth == k){            
            ins[failure_index].push_back(node_index);
            continue;
        }

        // Iterate through children and their outs, preserve sorted order

        /*size_t_max out_count = outs[node_index].size();
        size_t_max out = 0;
        while (out < out_count) {
            if (!has_in_edge[outs[node_index][out].second]) break;
            ++out;
        }//
        size_t_max in_count = ins[node_index].size();
        size_t_max in = 0;
        while (in < in_count) {
            if (!has_out_edge[ins[node_index][in]]) break;
            ++in;
        }

        if (in != in_count){
            size_t_max best_out = outs[node_index][out].second;
            size_t_max best_in = ins[node_index][in];
            successors[best_out] = best_in;

            has_out_edge[best_out] = true;
            has_in_edge[best_in] = true;
        }

        for (size_t_max j = in + 1; j < in_count; ++j){
            ins[failure_index].push_back(ins[node_index][j]);
        }
    }*/
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::print_stats(std::ostream &os)
{
    std::vector<size_t_max> depth_count(k + 1, 0);
    std::vector<size_t_max> real_depth_count(k + 1, 0);
    std::vector<size_t_max> real_depths(nodes.size(), 0);
    for (size_t_max i : BFS_indexes){
        depth_count[nodes[i].depth]++;
        if (i != 0) real_depths[i] = real_depths[nodes[i].parent] + 1;
        real_depth_count[real_depths[i]]++;
    }
    for (size_t_max i=0; i <= k; i++){
        std::cout << i << ":\t" << depth_count[i] << std::endl;
    }
    for (size_t_max i=0; i <= k; i++){
        std::cout << i << ":\t" << real_depth_count[i] << std::endl;
    }
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::print(std::ostream& os) {
    for (size_t_max i : BFS_indexes){
        os << i << ": ";
        nodes[i].print(os);
        os << std::endl;
    }
}

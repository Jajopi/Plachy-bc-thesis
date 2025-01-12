#pragma once
/*
#define kmer_index unsigned int
#define node_index unsigned int
#define kmar_max_len unsigned int
*/
#include <vector>
#include <queue>
#include <algorithm>
#include <limits>

#include "../kmers.h"

template<typename size_t_max>
constexpr size_t_max INVALID_NODE_TEMPLATE = std::numeric_limits<size_t_max>::max();

//in-place radix sort as in https://stackoverflow.com/questions/463105/in-place-radix-sort
template <typename kmer_t, typename size_t_max>
void radix_sort(std::vector<kmer_t>& kMers, size_t_max k, size_t_max begin, size_t_max end, size_t_max base = 0){
    if (begin >= end) return;

    size_t_max AEnd = begin, TBegin = end;
    size_t_max i = begin;
    while (i < TBegin){
        char nucleotide = NucleotideAtIndex(kMers[i], k, base);
        if (nucleotide == 'A') std::swap(kMers[i++], kMers[AEnd++]);
        else if (nucleotide == 'T') std::swap(kMers[i], kMers[--TBegin]);
        else ++i;
    }
    size_t_max CEnd = i = AEnd;
    while (i < TBegin){
        char nucleotide = NucleotideAtIndex(kMers[i], k, base);
        if (nucleotide == 'C') std::swap(kMers[i], kMers[CEnd++]);
        ++i;
    }

    if (base < k){
        radix_sort(kMers, k, begin, AEnd, (size_t_max)(base + 1));  // A
        radix_sort(kMers, k, AEnd, CEnd, (size_t_max)(base + 1));   // C
        radix_sort(kMers, k, CEnd, TBegin, (size_t_max)(base + 1)); // G
        radix_sort(kMers, k, TBegin, end, (size_t_max)(base + 1));  // T
    }
}

template <typename kmer_t, typename size_t_max>
void radix_sort(std::vector<kmer_t>& kMers, size_t_max k){
    radix_sort(kMers, k, size_t_max(0), size_t_max(kMers.size()), size_t_max(0));
}

template<typename size_t_max>
struct Node {
    size_t_max kmer_index;
    size_t_max depth;
    size_t_max parent;
    size_t_max failure;
    size_t_max inverse_failure;
    size_t_max child_range_begin;
    size_t_max child_range_end;

    using INVALID_NODE = INVALID_NODE_TEMPLATE<size_t_max>;

    Node(size_t_max kmer_index, size_t_max depth, size_t_max parent) :
        kmer_index(kmer_index), depth(depth), parent(parent),
        failure(INVALID_NODE), inverse_failure(INVALID_NODE),
        child_range_begin(INVALID_NODE), child_range_end(INVALID_NODE) {};

    Node(size_t_max kmer_index, size_t_max depth, size_t_max first_child, size_t_max last_child) :
        kmer_index(kmer_index), depth(depth), parent(INVALID_NODE),
        failure(INVALID_NODE), inverse_failure(INVALID_NODE),
        child_range_begin(first_child), child_range_end(last_child) {};
    
    Node() : kmer_index(INVALID_NODE) {}; // Special invalid type of node

    inline bool is_invalid () const {
        return kmer_index == INVALID_NODE;
    };

    void print(std::ostream& os) const;
};

template <typename size_t_max>
inline void Node<size_t_max>::print(std::ostream &os) const {
    os << depth << '/' << kmer_index << ", P: " << parent << ", F: " << failure << ", CH: [ ";
    //for (auto child : children) os << child << ' ';
    os << child_range_begin << " -- " << child_range_end << " ]";
}

template <typename kmer_t, typename size_t_max>
class HOGConstructer {
    std::vector<kmer_t>& kMers;
    std::vector<Node<size_t_max>> nodes;
    //std::vector<size_t_max> leaves; // TODO remove (not needed: BFS_indexes[-n:]) ... Or maybe not
    //std::vector<size_t_max> BFS_indexes; // TODO sort nodes in this order effectively and get rid of this
    std::vector<std::pair<size_t_max, size_t_max>> leaf_intervals;
    
    size_t_max k; // kmer-length
    size_t_max n; // number of kmers, number of leaves
    size_t_max full_depth = 0; // depth at which there are still all possible AC nodes;
    size_t_max INVALID_NODE = INVALID_NODE_TEMPLATE<size_t_max>;

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

    void precompute_full_depth();
    void construct_EHOG();
    void convert_EHOG_to_HOG();
    
    inline void create() {
        if (VERBOSE > 0) std::cout << "Constructing AC..." << std::endl;
        //construct_AC();
        if (VERBOSE > 1) print();
        if (VERBOSE > 0) std::cout << "Converting to EHOG..." << std::endl;
        //convert_AC_to_EHOG();
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
    /*size_t_max node_count = nodes.size();
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
    nodes.resize(keeped);*/
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::recreate_BFS_indexes(){
    /*if (VERBOSE > 0) std::cout << "Recreating BFS indexes...";

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

    if (VERBOSE > 0) std::cout << " done" << std::endl;*/
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::precompute_full_depth(){
    full_depth = k - 1;
    kmer_t max_diff = 1 << 2; //1 << (2 * (k - full_depth));

    kmer_t last_kmer = kMers[0];
    for (kmer_t actual_kmer : kMers){
        while (actual_kmer - last_kmer >= max_diff){
            --full_depth;
            max_diff <<= 2;
        }
        last_kmer = actual_kmer;
    }
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::construct_EHOG() {
    if (!nodes.empty()) throw std::invalid_argument("AC has already been constructed");

    radix_sort(kMers, k);
    precompute_full_depth();

    nodes.reserve((k - full_depth) * kMers.size()); // TODO decrease

    /*// Add leaves
    for (size_t_max l = kMers.size() - 1; l >= 0; --l){
        nodes.push_back(Node(l, k, INVALID_NODE, INVALID_NODE));
    }

    // Create other nodes, set children and parents
    //size_t_max node_count = kMers.size();
    size_t_max actual_depth = k;
    kmer_t max_diff = 1;

    for (size_t_max i = 0; i < node_count; ++i){
        if (nodes[i].depth != actual_depth){
            --actual_depth;
            max_diff <<= 2;

            if (actual_depth == full_depth){
                // Set parents of top nodes (in depth = full-depth)
                while (i < node_count) nodes[i++].parent = INVALID_NODE;
                break;
            };
        }

        size_t_max j = i;
        while (nodes[j].depth == actual_depth && // TODO can first part be removed?
               kMers[nodes[i].kmer_index] - kMers[nodes[j].kmer_index] < max_diff){
                    nodes[j].parent = node_count;
                    ++j;
                };
        
        nodes.push_back(Node(nodes[i].kmer_index, (size_t_max)(actual_depth - 1), i, j));
        ++node_count;
        i = j - 1;
    }

    nodes.shrink_to_fit();*/


    std::vector<Node<size_t_max>> current_nodes, last_nodes;
    current_nodes.reserve(n); // Might not be enought? Or?
    last_nodes.reserve(n); // Might not be enought? Or?
    
    std::vector<kmer_t> failures; failures.reserve(n);
    std::vector<size_t_max> last_failure_indexes; last_failure_indexes.reserve(n);
    
    size_t_max new_node_index = 0;

    // Add leaves
    for (size_t_max i = kMers.size() - 1; i >= 0; --i){
        last_nodes.push_back(Node(i, k, INVALID_NODE, INVALID_NODE));
    }
    new_node_index = kMers.size();
    failures = kMers; // Copy
    last_failure_indexes.resize(n); for (size_t_max i = 0; i < n; ++i) last_failure_indexes[i] = i;
    
    // Add other nodes
    kmer_t max_diff = 1 << 2;

    for (size_t_max depth = k - 1; depth >= full_depth; --depth){
        size_t_max saved_node_count = nodes.size();
        size_t_max last_node_count = last_nodes.size();
        size_t_max shift_from_removed = 0;

        for (size_t_max i = 0; i < n; ++i) failures[i] = BitSuffix(failures[i], depth * 2); // Is it really * 2 ?
        radix_sort(failures, k, size_t_max(0), n, size_t_max(k - depth));

        size_t_max failure_index = n - 1;
        for (size_t_max i = 0; i < last_node_count; ++i){
            /*if (inverse_failure_indexes[i] != INVALID_NODE){
                nodes[inverse_failure_indexes[i]].failure -= shift_from_removed;
            }*/

            size_t_max j = i;
            while (j < last_node_count &&
                    kMers[last_nodes[i].kmer_index] - kMers[last_nodes[j].kmer_index] < max_diff){
                nodes[j].parent = new_node_index;
                ++j;
            }

            kmer_t current_prefix = BitPrefix(kMers[i], k * 2, depth * 2); // Is it really * 2 ?
            while (failure_index >= 0 && current_prefix < failures[failure_index]) --failure_index;
            
            if (current_prefix == failures[failure_index]){
                current_nodes.push_back(Node(last_nodes[i].kmer_index,
                                                size_t_max(depth - 1),
                                                size_t_max(i - shift_from_removed),
                                                size_t_max(j - shift_from_removed)));
                
                size_t_max last_failure_index = last_failure_indexes[failure_index];
                if (last_failure_index < saved_node_count){
                    nodes[last_failure_index].failure = new_node_index;
                    
                }
                current_nodes.back().inverse_failure = last_failure_index;
                last_failure_indexes[failure_index] = new_node_index;
                //nodes[last_failure_indexes[failure_index]].failure = new_node_index;
                ++new_node_index;
            }
            else {
                for (size_t_max swap_index = i; swap_index < j; ++swap_index){
                    current_nodes.push_back(Node<size_t_max>());
                    std::swap(current_nodes.back(), last_nodes[swap_index]);
                }
                shift_from_removed += j - i;
                new_node_index += j - i;
            }

            i = j - 1;
        }

        
        for (Node node: last_nodes){
            if (!node.is_invalid()){
                node.parent -= shift_from_removed;
                //if (node.failure != INVALID_NODE) node.failure -= shift_from_removed;
                nodes.push_back(node);
            }
        }

        max_diff <<= 2;
        last_nodes = std::move(current_nodes);
    }

    for (Node node: last_nodes) nodes.push_back(node);

    // Create failure links

    size_t_max node_index = 0;
    failures = kMers;
    for (size_t_max depth = k; depth > full_depth; --depth){
        if (depth != k){
            for (size_t_max i = 0; i < n; ++i) failures[i] = BitSuffix(failures[i], depth * 2); // Is it really * 2 ?
            radix_sort(failures, k, size_t_max(0), n, size_t_max(k - depth));
        }

        for (size_t_max i = 0; i < n; ++i){
            if (i < n - 1 && failures[i] == failures[i + 1]) continue;

            nodes[node_index].failure;
        }

        
    }
    
    /*if (!nodes.empty()) throw std::invalid_argument("AC has already been constructed");

    // Add nodes and forward edges
    //nodes.resize(n); // TODO better guess
    nodes.emplace_back(0, 0, 0);
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
                nodes.emplace_back(kmer_index, depth + 1, node_index);
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
    }*/
}

//template <typename kmer_t, typename size_t_max>
//inline void HOGConstructer<kmer_t, size_t_max>::convert_AC_to_EHOG() {
    /*if (nodes.empty()) throw std::invalid_argument("Cannot convert empty AC to EHOG");

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
    recreate_BFS_indexes();*/
//}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::convert_EHOG_to_HOG() {
    /*if (nodes.empty()) throw std::invalid_argument("Cannot convert empty EHOG to HOG");

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
    recreate_BFS_indexes();*/
}

template <typename kmer_t, typename size_t_max>
inline std::vector<size_t_max> HOGConstructer<kmer_t, size_t_max>::compute_ordering(size_t_max new_run_score, size_t_max base_score) {
    return std::vector<size_t_max>();
    /*size_t_max node_count = nodes.size();
    
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
                    heap.push_back(std::make_pair(score, child_index))
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
    /*std::vector<size_t_max> depth_count(k + 1, 0);
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
    }*/
}

template <typename kmer_t, typename size_t_max>
inline void HOGConstructer<kmer_t, size_t_max>::print(std::ostream& os) {
    /*for (size_t_max i : BFS_indexes){
        os << i << ": ";
        nodes[i].print(os);
        os << std::endl;
    }*/
}

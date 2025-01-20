#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>

#include "../kmers.h"

template<typename size_t_max>
class UnionFind {
    std::vector<size_t_max> roots;
    size_t_max component_count;

public:
    UnionFind(size_t_max size) : roots(size), component_count(size) { for (size_t_max i = 0; i < size; ++i) roots[i] = i; };
    
    inline size_t_max find(size_t_max x) {
        size_t_max root = roots[x];
        if (roots[root] == root) return root;

        while (roots[root] != root) root = roots[root];
        
        while (x != root){
            size_t_max new_x = roots[x];
            roots[x] = root;
            x = new_x;
        }

        return root;
    }

    inline bool are_connected(size_t_max x, size_t_max y){
        return find(x) == find(y);
    }

    inline void connect(size_t_max from, size_t_max to){
        if (are_connected(from, to)) return;
        roots[from] = to;
        --component_count;
    }

    inline size_t_max count() const { return component_count; };
};


constexpr size_t K_SIZE_CONSTANT = 8;

template<typename size_t_max, size_t_max K_SIZE> // TODO replace with better suited value?
struct CS_AC_Node {
    size_t_max depth_and_kmer_index;            // Construction + searching
    size_t_max failure;                         // Construction + searching
    union {
        size_t_max child_range_begin;           // Construction
        size_t_max leaf_range_begin; // changes // Searching
    };
    union {
        size_t_max child_range_end;             // Construction
        size_t_max leaf_range_end;              // Searching
    };

    static inline size_t_max INVALID_NODE() { return std::numeric_limits<size_t_max>::max() >> K_SIZE; };
    
    inline size_t_max kmer_index() const { return depth_and_kmer_index >> K_SIZE; };
    inline size_t_max depth() const { return depth_and_kmer_index & ((1 << K_SIZE) - 1); };
    /*inline void reverse_kmer_index(size_t_max n) {
        depth_and_kmer_index = depth() + ((n - 1 - kmer_index()) << K_SIZE);
    }*/

    inline size_t_max original_leaf_range_begin(size_t_max N) const { return N - 1 - (depth_and_kmer_index >> K_SIZE); };

    CS_AC_Node(size_t_max kmer_index, size_t_max depth, size_t_max first_child, size_t_max after_last_child) :
        depth_and_kmer_index((kmer_index << K_SIZE) + depth), failure(INVALID_NODE()),
        child_range_begin(first_child), child_range_end(after_last_child) {};

    void print(std::ostream& os) const;
};

template<typename size_t_max, size_t_max K_SIZE>
inline void CS_AC_Node<size_t_max, K_SIZE>::print(std::ostream &os) const {
    os << depth() << '/';
    if (kmer_index() == INVALID_NODE()) os << "INV"; else os << kmer_index();
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
    using Node = CS_AC_Node<size_t_max, K_SIZE_CONSTANT>;

    std::vector<kmer_t> kMers;
    std::vector<Node> nodes;
    
    size_t_max K;            // kmer-length
    size_t_max N;            // number of kmers, number of leaves
    size_t_max DEPTH_CUTOFF; // first not-to-be-reached depth
    bool COMPLEMENTS;        // whether or not complements are used

    bool CONVERTED_TO_SEARCHABLE = false;
    size_t_max NEW_RUN_SCORE = 0;
    size_t_max BASE_EXTENSION_SCORE = 0;
    
    static inline size_t_max INVALID_NODE() { return Node::INVALID_NODE(); };

    void sort(std::vector<kmer_t>& kMers);

    void resort_and_shorten_failures(std::vector<std::pair<size_t_max, size_t_max>> &failures, size_t_max depth);

    void push_failure_of_node_into_hq(std::vector<std::tuple<size_t_max, size_t_max, size_t_max>>& hq,
                                      size_t_max node_index, size_t_max origin_leaf_index, size_t_max current_priority);
    
public:
    CuttedSortedAC(const std::vector<kmer_t>& kMers, size_t_max K, size_t_max DEPTH_CUTOFF = 0, bool complements = false) :
        kMers(kMers), K(K), N(kMers.size()), DEPTH_CUTOFF(DEPTH_CUTOFF), COMPLEMENTS(complements) {};

    CuttedSortedAC(std::vector<kmer_t>&& kMers, size_t_max K, size_t_max DEPTH_CUTOFF = 0, bool complements = false) :
        kMers(std::move(kMers)), K(K), N(kMers.size()), DEPTH_CUTOFF(DEPTH_CUTOFF), COMPLEMENTS(complements) {};

    void construct_graph();

    void convert_to_searchable_representation();
    void set_search_parameters(size_t_max new_run_score, size_t_max base_extension_score = 1);
    std::vector<size_t_max> compute_indexes(bool print = false);

    std::ostream& LOG_STREAM = std::cerr;
    void print_stats(std::ostream& os = std::cout);
    void print_sorted(std::ostream& os = std::cout);
    void print_topological(std::ostream& os = std::cout);
    void print_topological(std::ostream& os, size_t_max root, size_t_max depth);
};


template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::construct_graph() {
    if (!nodes.empty()){
        throw std::invalid_argument("Graph has already been constructed.");
    }

    LOG_STREAM << "Preparing..." << std::endl;
    sort(kMers);

    nodes.reserve((K - DEPTH_CUTOFF) * N);

    std::vector<Node> current_nodes; current_nodes.reserve(N);
    std::vector<std::pair<size_t_max, size_t_max>> failures(N); // kmer index, last index
    
    // Add leaves
    for (size_t_max i = 0; i < N; ++i){
        current_nodes.emplace_back(size_t_max(N - i - 1), K, i, i + 1);
        failures[i] = std::make_pair(i, N - i - 1);
    }

    LOG_STREAM << K << " --> " << DEPTH_CUTOFF << ": " << std::setfill(' ') << std::setw(3) << K; LOG_STREAM.flush();

    // Add other nodes
    for (size_t_max depth = K - 1; depth > DEPTH_CUTOFF; --depth){
        LOG_STREAM << "\b\b\b" << std::setfill(' ') << std::setw(3) << depth; LOG_STREAM.flush();
        
        size_t_max last_node_count = nodes.size();
        for (Node node: current_nodes) nodes.push_back(node);
        current_nodes.clear();
        size_t_max node_count = nodes.size();
        size_t_max new_node_index = node_count;

        resort_and_shorten_failures(failures, depth);
        size_t_max failure_index = failures.size() - 1;

        for (size_t_max i = last_node_count; i < node_count; ++i){
            kmer_t current_prefix = BitPrefix(kMers[nodes[i].kmer_index()], K, depth);
            while (failure_index > 0 && current_prefix < BitSuffix(kMers[failures[failure_index].first], depth)) --failure_index;

            bool new_node_on_failure_path = (current_prefix == BitSuffix(kMers[failures[failure_index].first], depth));

            size_t_max j = i;
            while (j < node_count && current_prefix == BitPrefix(kMers[nodes[j].kmer_index()], K, depth)) ++j;

            current_nodes.emplace_back(nodes[i].kmer_index(),
                                       size_t_max(depth),
                                       size_t_max(i),
                                       size_t_max(j));
            
            if (new_node_on_failure_path){
                nodes[failures[failure_index].second].failure = new_node_index;
                failures[failure_index].second = new_node_index;

                while (failure_index != 0 &&
                       current_prefix == BitSuffix(kMers[failures[failure_index - 1].first], depth)){
                    --failure_index;
                    
                    nodes[failures[failure_index].second].failure = new_node_index;
                    failures[failure_index] = std::make_pair(INVALID_NODE(), INVALID_NODE()); // Will be last, lost
                }
                --failure_index;
            }
            
            ++new_node_index;
            i = j - 1;
        }
    }
    LOG_STREAM << "\b\b\b" << std::setfill(' ') << std::setw(3) << DEPTH_CUTOFF << std::endl;

    for (Node node: current_nodes) nodes.push_back(node);

    nodes.emplace_back(0, 0, nodes.size() - current_nodes.size(), nodes.size());
    size_t_max root_node = nodes.size() - 1;
    
    for (size_t_max i = 0; i < root_node; ++i){
        if (nodes[i].failure == INVALID_NODE()) nodes[i].failure = root_node;
    }

    LOG_STREAM << "Graph construction finished." << std::endl;
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::convert_to_searchable_representation() {
    if (nodes.empty()){
        throw std::invalid_argument("Graph has not been constructed yet.");
    }
    if (CONVERTED_TO_SEARCHABLE){
        throw std::invalid_argument("Graph has already been converted to searchable.");
    }

    LOG_STREAM << "Converting graph to searchable..." << std::endl;

    size_t_max node_count = nodes.size();
    for (size_t_max i = N; i < node_count; ++i){
        nodes[i].leaf_range_begin = nodes[nodes[i].child_range_begin].child_range_begin;
        nodes[i].leaf_range_end = nodes[nodes[i].child_range_end - 1].child_range_end;
    }

    CONVERTED_TO_SEARCHABLE = true;
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::set_search_parameters(size_t_max new_run_score, size_t_max base_extension_score) {
    NEW_RUN_SCORE = new_run_score;
    BASE_EXTENSION_SCORE = base_extension_score;
}

template <typename kmer_t, typename size_t_max>
inline std::vector<size_t_max> CuttedSortedAC<kmer_t, size_t_max>::compute_indexes(bool print) {
    if (!CONVERTED_TO_SEARCHABLE){
        throw std::invalid_argument("Graph has not been converted to searchable yet.");
    }
    if (NEW_RUN_SCORE == 0){
        throw std::invalid_argument("Mandatory parameter NEW_RUN_SCORE was not set.");
    }
    if (BASE_EXTENSION_SCORE == 0){
        throw std::invalid_argument("Mandatory parameter BASE_EXTENSION_SCORE was not set.");
    }

    LOG_STREAM << "Computing indexes..." << std::endl;

    std::vector<std::tuple<size_t_max, size_t_max, size_t_max>> hq; // priority, current node_index, origin_leaf_index
    hq.reserve(2 * N); // TODO better guess
    //for (size_t_max i = 0; i < N; ++i) hq.emplace_back(std::numeric_limits<size_t_max>::max(), i, i);
    //std::make_heap(hq.begin(), hq.end()); // Already is a heap
    for (size_t_max i = 0; i < N; ++i) push_failure_of_node_into_hq(hq, i, i, std::numeric_limits<size_t_max>::max());

    size_t_max root_node_index = nodes.size() - 1;
    size_t_max resigned_count = 0;
    size_t_max last_origin_leaf = INVALID_NODE(), last_node = INVALID_NODE();

    std::vector<size_t_max> previous(N, INVALID_NODE()); // leaf number + N if new run was started, else leaf number
    UnionFind components(N);

    while (!hq.empty() && components.count() - resigned_count > 1){
        auto p = hq.front(); std::pop_heap(hq.begin(), hq.end()); hq.pop_back();
        // LOG_STREAM << std::get<0>(p) << ' ' << std::get<1>(p) << ' ' << std::get<2>(p) << std::endl;

        size_t_max origin_leaf_index = std::get<2>(p);
        Node& origin_leaf_node = nodes[origin_leaf_index];
        if (origin_leaf_node.leaf_range_begin == origin_leaf_node.leaf_range_end) continue; // Already found next or resigned for this leaf

        size_t_max priority = std::get<0>(p), node_index = std::get<1>(p);

        if (last_origin_leaf == origin_leaf_index && last_node == node_index) continue; // Skipping duplicates, TODO better?
        last_origin_leaf = origin_leaf_index; last_node = node_index;

        if (node_index == root_node_index){
            --origin_leaf_node.leaf_range_end; // We resign for current leaf
            ++resigned_count;
            continue;
        }

        Node& node = nodes[node_index];
        while (node.leaf_range_begin != node.leaf_range_end &&
                (previous[node.leaf_range_begin] != INVALID_NODE() || // Todo make efficient
                components.are_connected(node.leaf_range_begin, origin_leaf_index)))
            ++node.leaf_range_begin;
        
        // LOG_STREAM << node.original_leaf_range_begin(N) << ' ' << node.leaf_range_begin << ' ' << node.leaf_range_end << std::endl;

        if (node.leaf_range_begin != node.leaf_range_end){
            previous[node.leaf_range_begin] = origin_leaf_index; // Found previous for this leaf
            ++nodes[origin_leaf_index].leaf_range_begin; // Found next for origin leaf
            // LOG_STREAM << origin_leaf_index << " -> " << node.leaf_range_begin << std::endl;

            components.connect(origin_leaf_index, node.leaf_range_begin); // First one pointing to the second one = always point to last in a chain
            continue;
        }

        for (size_t_max i = node.original_leaf_range_begin(N); i < node.leaf_range_end; ++i){
            if (components.are_connected(origin_leaf_index, i) || nodes[i].leaf_range_end == i) continue; // Is already connected to origin leaf or was resigned for
            push_failure_of_node_into_hq(hq, i, origin_leaf_index, priority); // Add failure of that leaf
        }

        push_failure_of_node_into_hq(hq, node_index, origin_leaf_index, priority); // Add failure of current node
    }

    hq.clear(); hq.shrink_to_fit(); // Should deallocate, ?


    std::vector<size_t_max> indexes;
    indexes.resize(N);
    
    size_t_max current_final_index = N - 1;
    for (size_t_max i = 0; i < N; ++i){
        if (components.find(i) == i){ // Begining of a chain
            size_t_max actual = i;
            while (actual != INVALID_NODE()){
                indexes[current_final_index--] = N - actual - 1;
                // LOG_STREAM << actual << std::endl;
                actual = previous[actual];
            }
            // LOG_STREAM << std::endl;
        }
    }
    LOG_STREAM << "Resigned leaves: " << resigned_count << std::endl;

    return indexes;
}

// Sorting

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::sort(std::vector<kmer_t> &kMers) {
    for (size_t_max i = 1; i <= N; ++i){
        if (i == N){
            LOG_STREAM << "Already sorted.";
            return;
        }
        if (kMers[i] < kMers[i - 1]) break;
    }

    LOG_STREAM << "Using std::sort..." << std::endl;
    std::sort(kMers.begin(), kMers.end());
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::resort_and_shorten_failures(std::vector<std::pair<size_t_max, size_t_max>> &failures, size_t_max depth) {
    if (CONVERTED_TO_SEARCHABLE){
        throw std::invalid_argument("Graph has already been converted to searchable.");
    }

    size_t_max count = failures.size();

    size_t_max skipped = 0;
    for (size_t_max i = 0; i < count; ++i){
        if (failures[i].second == INVALID_NODE()) ++skipped;
        else failures[i - skipped] = failures[i];
    }

    count -= skipped;
    failures.resize(count);
    std::vector<std::pair<size_t_max, size_t_max>> sorted_failures(count);

    size_t_max end_indexes[4] = {0, 0, 0, 0}; // Starts of: A, C, G, T
    for (size_t_max i = 0; i < count; ++i){
        uint8_t base_index = NucleotideIndexAtIndex(kMers[failures[i].first], K, K - depth - 1);
        ++end_indexes[base_index];
    }

    size_t_max start_indexes[4] = {0, 0, 0, 0};
    for (uint8_t i = 1; i < 4; ++i){
        end_indexes[i] += end_indexes[i - 1];
        start_indexes[i] = end_indexes[i - 1];
    }

    for (size_t_max s = 0; s < count; ++s){
        uint8_t best_i = 4;
        kmer_t best_suffix;
        for (uint8_t i = 0; i < 4; ++i){
            if (start_indexes[i] == end_indexes[i]) continue;
            if (best_i == 4 || best_suffix > BitSuffix(kMers[failures[start_indexes[i]].first], depth)){
                best_i = i;
                best_suffix = BitSuffix(kMers[failures[start_indexes[best_i]].first], depth);
            }
        }
        sorted_failures[s] = failures[start_indexes[best_i]++];
    }

    failures = std::move(sorted_failures);
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::push_failure_of_node_into_hq(std::vector<std::tuple<size_t_max, size_t_max, size_t_max>> &hq,
                                                                             size_t_max node_index, size_t_max origin_leaf_index, size_t_max current_priority) {
    Node& node = nodes[node_index];
    Node& failure_node = nodes[node.failure];

    // LOG_STREAM << node_index << ' ' << node.failure << ' ' << std::endl;

    size_t_max failure_priority = current_priority - (node.depth() - failure_node.depth()) * BASE_EXTENSION_SCORE;
    if (node.depth() == K - 1 || (node.depth() == K && failure_node.depth() < K - 1)) failure_priority -= NEW_RUN_SCORE;

    hq.push_back(std::make_tuple(failure_priority, node.failure, origin_leaf_index)); std::push_heap(hq.begin(), hq.end());
}

// Printing

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::print_stats(std::ostream &os)
{
    os << "Computing stats..." << std::endl;
    std::vector<size_t_max> depth_count(K + 1, 0);
    
    size_t_max node_count = nodes.size();
    for (size_t_max i = 0; i < node_count; ++i){
        depth_count[nodes[i].depth()]++;
    }

    os << "Depths:" << std::endl;
    for (size_t_max i = 0; i <= K; i++){
        if (depth_count[i] == 0) continue;
        os << i << ":\t" << depth_count[i] << std::endl;
    }
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::print_sorted(std::ostream& os) {
    size_t_max i = 0;
    for (Node node : nodes){
        os << i++ << ":\t";
        node.print(os);
        os << ":\t";
        for (size_t_max c = 0; c < node.depth(); ++c) os << NucleotideAtIndex(kMers[node.kmer_index], K, c);
        os << std::endl;
    }
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::print_topological(std::ostream& os) {
    if (CONVERTED_TO_SEARCHABLE){
        throw std::invalid_argument("Graph has already been converted to searchable.");
    }
    
    print_topological(os, nodes.size() - 1, 0);
}

template <typename kmer_t, typename size_t_max>
inline void CuttedSortedAC<kmer_t, size_t_max>::print_topological(std::ostream& os, size_t_max root, size_t_max depth) {
    Node node = nodes[root];
    
    for (size_t_max i = 0; i < depth; ++i) os << "|  ";
    os << "+->";
    
    os << root << ":\t";
    node.print(os);
    os << ":\t";
    for (size_t_max c = 0; c < node.depth(); ++c) os << NucleotideAtIndex(kMers[node.kmer_index()], K, c);
    
    os << std::endl;

    if (root < N) return;
    
    for (size_t_max i = node.child_range_begin; i < node.child_range_end; ++i){
        if (i == INVALID_NODE()) break;
        print_topological(os, i, depth + 1);
    }
}

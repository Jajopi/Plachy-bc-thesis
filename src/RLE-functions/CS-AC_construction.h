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
    UnionFind(size_t_max size) : roots(size), component_count(size) {
        for (size_t_max i = 0; i < size; ++i) roots[i] = i;
    };
    UnionFind() = default;
    
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

    inline void connect(size_t_max to, size_t_max from){ // Second one points to the first one - points to the begining of a chain
        if (are_connected(from, to)) return;
        roots[from] = to;
        --component_count;
    }

    inline size_t_max count() const { return component_count; };
};


template<typename size_t_max, size_t_max K_BIT_SIZE> // K_BIT_SIZE must be >= 3
struct CS_AC_Node {
    union {
        size_t_max depth_and_kmer_index;        // Construction + searching - original_leaf_range_begin comes from constructing
        size_t_max complement_index;            // Leaves while searching
    };
    size_t_max failure;                         // Construction + searching + leaves while searching
    union {
        size_t_max leaf_range_begin; // changes // Construction, searching
        size_t_max bitmask_and_next;            // Leaves while searching
    };

    static inline size_t_max INVALID_LEAF() { return std::numeric_limits<size_t_max>::max() >> K_BIT_SIZE; };
    static inline size_t_max INVALID_NODE() { return std::numeric_limits<size_t_max>::max(); };
    
    CS_AC_Node(size_t_max kmer_index, size_t_max depth, size_t_max first_child) :
        depth_and_kmer_index((kmer_index << K_BIT_SIZE) + depth), failure(INVALID_NODE()),
        leaf_range_begin(first_child) {};
    
    void print(std::ostream& os) const; // Construction

    inline size_t_max kmer_index() const { return depth_and_kmer_index >> K_BIT_SIZE; }; // Construction
    inline size_t_max depth() const { return depth_and_kmer_index & ((1 << K_BIT_SIZE) - 1); }; // Construction + searching

    inline size_t_max original_leaf_range_begin() const { return kmer_index(); }; // Searching, does not change through search

    // ... -> node = was used, node -> ... = was completed
    // Leaves while searching - bitmask:
    // 1 - was used
    // 2 - complement was chosen (used / completed)
    inline bool used() const { return (bitmask_and_next & size_t_max(1)) != 0; };
    inline void set_used() { bitmask_and_next |= size_t_max(1); };
    inline bool complement_chosen() const { return (bitmask_and_next & size_t_max(2)) != 0; };
    inline void set_complement_chosen() { bitmask_and_next |= size_t_max(2); };

    inline size_t_max next() const { return bitmask_and_next >> K_BIT_SIZE; };
    inline void set_next(size_t_max next) {
        bitmask_and_next &= size_t_max(3);
        bitmask_and_next |= (next << K_BIT_SIZE);
    };
    inline bool completed() const { return next() != INVALID_LEAF(); };
    inline void reset_bitmask_and_next() { bitmask_and_next = (INVALID_LEAF() << K_BIT_SIZE); };
};

template<typename size_t_max, size_t_max K_BIT_SIZE>
inline void CS_AC_Node<size_t_max, K_BIT_SIZE>::print(std::ostream &os) const {
    os << depth() << '/';
    if (kmer_index() == INVALID_LEAF()) os << "INV"; else os << kmer_index();
    os << ",\tF: ";
    if (failure == INVALID_NODE()) os << "INV"; else os << failure;
    os << ",\tCH: [ ";
    if (leaf_range_begin == INVALID_NODE()) os << "INV"; else os << leaf_range_begin;
    os << " - ";
}


template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
class CuttedSortedAC {
    using Node = CS_AC_Node<size_t_max, K_BIT_SIZE>;
    static inline size_t_max INVALID_LEAF() { return Node::INVALID_LEAF(); };
    static inline size_t_max INVALID_NODE() { return Node::INVALID_NODE(); };

    std::vector<kmer_t> kMers;  // sorted
    std::vector<Node> nodes;    // sorted by depth from K to DEPTH_CUTOFF, lexicographically
    
    size_t_max K;               // kmer-length
    size_t_max N;               // number of kmers, number of leaves
    size_t_max DEPTH_CUTOFF;    // first not-to-be-reached depth
    size_t_max PRACTICAL_DEPTH; // used to resize nodes
    size_t_max SEARCH_CUTOFF;   // fraction of uncompleted leaves to switch to fast mode
    bool COMPLEMENTS;           // whether or not complements are used

    bool CONVERTED_TO_SEARCHABLE = false;
    bool COMPUTED_RESULT = false;
    size_t_max RUN_PENALTY = INVALID_NODE();
    size_t_max EXTENSION_PENALTY = INVALID_NODE();

    std::vector<std::tuple<size_t_max, size_t_max, size_t_max>> stack;
    UnionFind<size_t_max> components;
    std::vector<size_t_max> backtracks;
    std::vector<size_t_max> backtrack_indexes;
    std::vector<size_t_max> previous;
    std::vector<size_t_max> remaining_priorities;

    void sort_and_remove_duplicate_kmers();

    void resort_and_shorten_failures(std::vector<std::pair<size_t_max, size_t_max>> &failures, size_t_max depth);

    bool try_complete_leaf(size_t_max leaf_index, size_t_max priority_drop_limit);
    void push_failure_of_node_into_stack(size_t_max riority, size_t_max node_index, size_t_max priority_drop_limit, size_t_max last_leaf);
    void squeeze_uncompleted_leaves(std::vector<size_t_max>& unclompleted_leaves);
    void set_backtrack_path_for_leaf(size_t_max origin_leaf, size_t_max next_leaf);
    
    size_t_max find_complement_kmer_index(size_t_max kmer_index);
public:
    CuttedSortedAC(const std::vector<kmer_t>& kmers, size_t_max K, size_t_max DEPTH_CUTOFF = 0, bool complements = false) :
        kMers(kmers), K(K), N(kMers.size()), DEPTH_CUTOFF(DEPTH_CUTOFF), PRACTICAL_DEPTH(K - DEPTH_CUTOFF),
        COMPLEMENTS(complements) { sort_and_remove_duplicate_kmers(); };
    CuttedSortedAC(std::vector<kmer_t>&& kmers, size_t_max K, size_t_max DEPTH_CUTOFF = 0, bool complements = false) :
        kMers(std::move(kmers)), K(K), N(kMers.size()), DEPTH_CUTOFF(DEPTH_CUTOFF), PRACTICAL_DEPTH(K - DEPTH_CUTOFF),
        COMPLEMENTS(complements) { sort_and_remove_duplicate_kmers(); };

    CuttedSortedAC(const std::vector<kmer_t>& kmers, size_t_max K, size_t_max DEPTH_CUTOFF = 0, size_t_max PRACTICAL_DEPTH = 0, bool complements = false) :
        kMers(kmers), K(K), N(kMers.size()), DEPTH_CUTOFF(DEPTH_CUTOFF), PRACTICAL_DEPTH(PRACTICAL_DEPTH),
        COMPLEMENTS(complements) { sort_and_remove_duplicate_kmers(); };
    CuttedSortedAC(std::vector<kmer_t>&& kmers, size_t_max K, size_t_max DEPTH_CUTOFF = 0, size_t_max PRACTICAL_DEPTH = 0, bool complements = false) :
        kMers(std::move(kmers)), K(K), N(kMers.size()), DEPTH_CUTOFF(DEPTH_CUTOFF), PRACTICAL_DEPTH(PRACTICAL_DEPTH),
        COMPLEMENTS(complements) { sort_and_remove_duplicate_kmers(); };

    void construct_graph();

    void convert_to_searchable_representation();
    void set_search_parameters(size_t_max run_penalty,
                               size_t_max extension_penalty = 1,
                               size_t_max precision = 10);
    void compute_result();
    size_t_max print_result(std::ostream& os);

    std::ostream& LOG_STREAM = std::cerr;
    void print_stats(std::ostream& os = std::cerr);
    void print_sorted(std::ostream& os = std::cerr);
    void print_topological(std::ostream& os = std::cerr);
    void print_topological(std::ostream& os, size_t_max root, size_t_max depth);
};

// Constructing

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::construct_graph() {
    if (!nodes.empty()){
        throw std::invalid_argument("Graph has already been constructed.");
    }

    // LOG_STREAM << "Constructing the graph..." << std::endl;

    nodes.reserve(PRACTICAL_DEPTH * N);

    std::vector<std::pair<size_t_max, size_t_max>> failures(N); // kmer index, last index
    
    // Add leaves
    for (size_t_max i = 0; i < N; ++i){
        nodes.emplace_back(i, K, i);
        failures[i] = std::make_pair(i, i);
    }
    size_t_max node_count = 0;

    LOG_STREAM << K << " --> " << DEPTH_CUTOFF << ": " << std::setfill(' ') << std::setw(3) << K; LOG_STREAM.flush();

    // Add other nodes
    // size_t_max remaining_depth_steps = K - DEPTH_CUTOFF - 1;
    // for (size_t_max depth = K - 1; depth > 0; depth -= depth / remaining_depth_steps--){
    for (size_t_max depth = K - 1; depth > DEPTH_CUTOFF; --depth){
        LOG_STREAM << "\b\b\b" << std::setfill(' ') << std::setw(3) << depth; LOG_STREAM.flush();
        
        size_t_max last_node_count = node_count;
        node_count = nodes.size();
        size_t_max new_node_index = node_count;

        resort_and_shorten_failures(failures, depth);
        size_t_max failure_count = failures.size();
        size_t_max failure_index = 0;

        for (size_t_max i = last_node_count; i < node_count; ++i){
            kmer_t current_prefix = BitPrefix(kMers[nodes[i].kmer_index()], K, depth);
            while (failure_index < failure_count &&
                   current_prefix > BitSuffix(kMers[failures[failure_index].first], depth))
                ++failure_index;

            bool new_node_on_failure_path = (failure_index < failure_count &&
                    current_prefix == BitSuffix(kMers[failures[failure_index].first], depth));

            size_t_max j = i;
            while (j < node_count && current_prefix == BitPrefix(kMers[nodes[j].kmer_index()], K, depth)) ++j;

            nodes.emplace_back(nodes[i].kmer_index(), depth, nodes[i].kmer_index());
            
            if (new_node_on_failure_path){
                nodes[failures[failure_index].second].failure = new_node_index;
                failures[failure_index].second = new_node_index;

                ++failure_index;
                while (failure_index < failure_count &&
                       current_prefix == BitSuffix(kMers[failures[failure_index].first], depth)){
                    
                    nodes[failures[failure_index].second].failure = new_node_index;
                    failures[failure_index] = std::make_pair(INVALID_LEAF(), INVALID_NODE()); // Will be last, lost
                    
                    ++failure_index;
                }
            }
            
            ++new_node_index;
            i = j - 1;
        }
    }
    LOG_STREAM << "\b\b\b" << std::setw(3) << DEPTH_CUTOFF << std::endl;

    nodes.emplace_back(0, 0, 0);
    
    size_t_max root_node = nodes.size() - 1;
    for (size_t_max i = 0; i < root_node; ++i){
        if (nodes[i].failure == INVALID_NODE()) nodes[i].failure = root_node; // Set all not-yet-set failures
    }

    // LOG_STREAM << "Graph construction finished." << std::endl;
}

// Converting

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::convert_to_searchable_representation() {
    if (nodes.empty()){
        throw std::invalid_argument("Graph has not been constructed yet.");
    }
    if (CONVERTED_TO_SEARCHABLE){
        throw std::invalid_argument("Graph has already been converted to searchable.");
    }

    // LOG_STREAM << "Converting graph to searchable..." << std::endl;

    for (size_t_max i = 0; i < N; ++i){ // Convert leaves
        Node& node = nodes[i];
        node.reset_bitmask_and_next();
        if (COMPLEMENTS) node.complement_index = find_complement_kmer_index(node.kmer_index());
        else node.complement_index = INVALID_LEAF(); // Let that crash if used somewhere
    }

    CONVERTED_TO_SEARCHABLE = true;
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::set_search_parameters(
            size_t_max run_penalty, size_t_max extension_penalty, size_t_max precision) {
    RUN_PENALTY = run_penalty;
    EXTENSION_PENALTY = extension_penalty;
    SEARCH_CUTOFF = N / std::max(2, 1 << precision);

    LOG_STREAM << "Run penalty: " << run_penalty << std::endl << "Precision: " << precision << std::endl;
}

// Sorting

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::sort_and_remove_duplicate_kmers() {
    // LOG_STREAM << "Sorting kmers..." << std::endl;
    std::sort(kMers.begin(), kMers.end());
    
    if (COMPLEMENTS){
        LOG_STREAM << "Removing duplicate kmers... " << N;

        size_t_max shift = 0;
        for (size_t_max i = 1; i < N; ++i){
            if (kMers[i] == kMers[i - 1]) ++shift;
            else kMers[i - shift] = kMers[i];
        }

        N -= shift;
        kMers.resize(N);

        LOG_STREAM << " -> " << N << std::endl;
    }
}

// Internal functions

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::resort_and_shorten_failures(
            std::vector<std::pair<size_t_max, size_t_max>> &failures, size_t_max depth) {
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

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline size_t_max CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::find_complement_kmer_index(size_t_max kmer_index) {
    kmer_t complement = ReverseComplement(kMers[kmer_index], K);
    
    size_t_max begin = 0, end = N - 1;
    while (begin < end){
        size_t_max middle = (begin + end) / 2;
        if (kMers[middle] == complement) return middle;
        else if (kMers[middle] < complement) begin = middle + 1;
        else end = middle - 1;
    }
    // print_kmer(kMers[kmer_index], K, LOG_STREAM); LOG_STREAM << ' ';
    // print_kmer(kMers[begin], K, LOG_STREAM); LOG_STREAM << std::endl;
    return begin; // Complement should always be present
}

// Debug printing

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::print_stats(std::ostream &os)
{
    if (CONVERTED_TO_SEARCHABLE){
        throw std::invalid_argument("Graph has already been converted to searchable and stats cannot be computed.");
    }

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

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::print_sorted(std::ostream& os) {
    if (CONVERTED_TO_SEARCHABLE){
        throw std::invalid_argument("Graph has already been converted to searchable and cannot be printed.");
    }
    
    size_t_max i = 0;
    for (Node node : nodes){
        os << i++ << ":\t";
        node.print(os);
        os << ":\t";
        for (size_t_max c = 0; c < node.depth(); ++c) os << NucleotideAtIndex(kMers[node.kmer_index()], K, c);
        os << std::endl;
    }
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::print_topological(std::ostream& os) {
    if (CONVERTED_TO_SEARCHABLE){
        throw std::invalid_argument("Graph has already been converted to searchable and cannot be printed.");
    }
    
    print_topological(os, nodes.size() - 1, 0);
}

template <typename kmer_t, typename size_t_max, size_t_max K_BIT_SIZE>
inline void CuttedSortedAC<kmer_t, size_t_max, K_BIT_SIZE>::print_topological(std::ostream& os, size_t_max root, size_t_max depth) {
    Node node = nodes[root];
    
    for (size_t_max i = 0; i < depth; ++i) os << "|  ";
    os << "+->";
    
    os << root << ":\t";
    node.print(os);
    os << ":\t";
    for (size_t_max c = 0; c < node.depth(); ++c) os << NucleotideAtIndex(kMers[node.kmer_index()], K, c);
    
    os << std::endl;

    if (root < N) return;
    
    /*for (size_t_max i = node.leaf_range_begin; i < node.child_range_end; ++i){
        if (i == INVALID_NODE()) break;
        print_topological(os, i, depth + 1);
    }*/
}

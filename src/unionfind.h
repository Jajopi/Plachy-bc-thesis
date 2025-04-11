#pragma once

#include <vector>

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

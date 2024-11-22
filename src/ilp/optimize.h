
#define GUROBI_DIR "/opt/gurobi1103/linux64/include/gurobi_c++.h"

#include GUROBI_DIR
#include <vector>
#include <sstream>
#include <math.h>

template <typename kmer_t>
std::vector<size_t> optimize_indexes(const std::vector<kmer_t>& kmers, int (*distance)(const std::vector<kmer_t>&, size_t, size_t, size_t), size_t k){
    try {
        const size_t n = kmers.size();

        std::cout << n << std::endl;

        // n k-mers + 1 s_0
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        std::cout << "Gurobi env + model created" << std::endl;

        // Edge variables
        GRBVar edges[n + 1][n + 1];
        for (size_t i = 0; i < n + 1; i++){
            for (size_t j = 0; j < n + 1; j++){
                if (i == j) continue;

                std::string name = "Edge_" + std::to_string(i) + "_" + std::to_string(j);
                if (i == n || j == n)
                    edges[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                else
                    edges[i][j] = model.addVar(0.0, 1.0, distance(kmers, k, i, j), GRB_BINARY, name);
            }
        }
        // Edge constraints
        // Out-degrees of i
        for (size_t i = 0; i < n + 1; i++){
            GRBLinExpr outdegree = 0;
            for (size_t j = 0; j < n + 1; j++){
                if (i == j) continue;

                outdegree += edges[i][j];
            }
            std::string name = "OutDegree_" + std::to_string(i);
            model.addConstr(outdegree == 1.0, name);
        }
        // In-degrees of j
        for (size_t j = 0; j < n + 1; j++){
            GRBLinExpr indegree = 0;
            for (size_t i = 0; i < n + 1; i++){
                if (i == j) continue;

                indegree += edges[i][j];
            }
            std::string name = "InDegree_" + std::to_string(j);
            model.addConstr(indegree == 1.0, name);
        }
        
        std::cout << "Edge constraints created" << std::endl;

        // Index variables (u)
        GRBVar indexes[n + 1];
        for (size_t i = 0; i < n + 1; i++){
            std::string name = "Index_" + std::to_string(i);
            indexes[i] = model.addVar(0.0, n, 0.0, GRB_INTEGER, name);
        }
        // Index constraints
        // Indexes ascend
        for (size_t i = 0; i < n; i++){
            for (size_t j = 0; j < n; j++){
                if (i == j) continue;

                GRBLinExpr difference = indexes[i] - indexes[j] + 1 - n + (n * edges[i][j]);
                std::string name = "IndexDiff_" + std::to_string(i) + "_" + std::to_string(j);
                model.addConstr(difference <= 0.0, name);
            }
        }

        std::cout << "Index constraints created" << std::endl;

        // Does it minimize?
        model.optimize();

        std::cout << "Optimization performed" << std::endl;

        std::vector<size_t> optimized_indexes(n);
        for (size_t i = 0; i < n; i++){
            optimized_indexes[i] = round(indexes[i].get(GRB_DoubleAttr_X));
        }
        
        std::cout << "Indexes retreived" << std::endl;

        return optimized_indexes;
    }
    catch(GRBException e){
        std::cout << e.getErrorCode() << std::endl << e.getMessage() << std::endl;
        return std::vector<size_t>();
    }
}

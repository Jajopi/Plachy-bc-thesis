
#define GUROBI_DIR "/opt/gurobi1103/linux64/include/gurobi_c++.h"
#define MEMORY_LIMIT_GB 16
#define THREAD_COUNT 1

#include GUROBI_DIR
#include <vector>
#include <sstream>
#include <math.h>

template <typename kmer_t>
std::vector<size_t> optimize_indexes_with_complements(const std::vector<kmer_t>& kMers,
        int (*distance)(const std::vector<kmer_t>&, size_t, size_t, size_t),
        size_t k){
    try {
        // 2 n k-mers + 1 s_0 node = 2n + 1 nodes
        // each node has its reverse complement on position + n
        const size_t nn = kMers.size();
        const size_t n = nn / 2;

        // Creating model
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        // Setting model parameters
        model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_OFF); // seems to not help in this case
        model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER); // consumes less memory than concurent, barrier seems to be the best algortihm for this case
        //model.set(GRB_DoubleParam_MemLimit, MEMORY_LIMIT_GB); // prevents the computer from freezing, maybe?
        //model.set(GRB_IntParam_Threads, THREAD_COUNT); // less threads use less memory, barrier seems to never use more than one anyway

        // Edge variables
        GRBVar** edges = new GRBVar*[nn + 1];
        for (size_t i = 0; i < nn + 1; i++){
            edges[i] = new GRBVar[nn + 1];
            for (size_t j = 0; j < nn + 1; j++){
                if ((i < n ? i : i - n) == (j < n ? j : j - n)) continue;

                std::string name = "Edge_" + std::to_string(i) + "_" + std::to_string(j);
                if (i == nn || j == nn)
                    edges[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                else
                    edges[i][j] = model.addVar(0.0, 1.0, distance(kMers, k, i, j), GRB_BINARY, name);
            }
        }
        // Edge constraints
        // Outdegrees and Indegrees, first half also with complements
        for (size_t i = 0; i < nn + 1; i++){
            GRBLinExpr outdegree = 0;
            GRBLinExpr indegree = 0;
            GRBLinExpr complementOutdegree = 0;
            GRBLinExpr complementIndegree = 0;
            for (size_t j = 0; j < nn + 1; j++){
                if ((i < n ? i : i - n) == (j < n ? j : j - n)) continue;

                outdegree += edges[i][j];
                indegree += edges[j][i];
                if (i < n){
                    complementOutdegree += edges[i + n][j];
                    complementIndegree += edges[j][i + n];
                }
                
            }
            std::string name;
            /*name = "OutDegree_" + std::to_string(i);
            model.addConstr(outdegree <= 1.0, name);*/

            /*name = "InDegree_" + std::to_string(i);
            model.addConstr(indegree <= 1.0, name);*/

            name = "InOutSameComplement_" + std::to_string(i);
            model.addConstr(outdegree - indegree == 0.0, name);

            if (i < n || i == nn){
                name = "TotalOutDegree_" + std::to_string(i);
                model.addConstr(outdegree + complementOutdegree == 1.0, name);

                name = "TotalInDegree_" + std::to_string(i);
                model.addConstr(indegree + complementIndegree == 1.0, name);
            }
        }
        /*// In-degrees of j, first half also with complements
        for (size_t j = 0; j < nn + 1; j++){
            GRBLinExpr indegree = 0;
            GRBLinExpr complementIndegree = 0;
            for (size_t i = 0; i < nn + 1; i++){
                if ((i < n ? i : i - n) == (j < n ? j : j - n)) continue;

                indegree += edges[i][j];
                if (i != nn && j < n) complementIndegree += edges[i][j + n];
            }
            std::string name = "InDegree_" + std::to_string(j);
            model.addConstr(indegree <= 1.0, name);

            if (j < n){
                name = "TotalInDegree_" + std::to_string(j);
                model.addConstr(indegree + complementIndegree == 1.0, name);
            }
        }*/
        
        // Index variables (u)
        GRBVar* indexes = new GRBVar[n + 1];
        for (size_t i = 0; i < n; i++){
            std::string name = "Index_" + std::to_string(i);
            indexes[i] = model.addVar(0.0, n, 0.0, GRB_INTEGER, name);
        }
        // Index constraints
        // Indexes ascend
        for (size_t i = 0; i < n; i++){
            for (size_t j = 0; j < n; j++){
                if (i == j) continue;

                GRBLinExpr difference = indexes[i] - indexes[j] +
                                        1 - n + n * (
                                            edges[i][j] + edges[i + n][j + n] +
                                            edges[i + n][j] + edges[i][j + n]);
                std::string name = "IndexDiff_" + std::to_string(i) + "_" + std::to_string(j);
                model.addConstr(difference <= 0.0, name);
            }
        }


        // Run optimization
        model.optimize();


        // Get results
        std::vector<size_t> optimized_indexes(n);
        for (size_t i = 0; i < n; i++){
            optimized_indexes[round(indexes[i].get(GRB_DoubleAttr_X))] = i;
        }

        for (size_t i = 1; i < n; i++){
            size_t edge_begin = optimized_indexes[i - 1];
            size_t edge_end = optimized_indexes[i];
            //std::cout << edge_begin << " " << edge_end << ": ";
            //std::cout << edges[edge_begin][edge_end].get(GRB_DoubleAttr_X) << " " <<
            //    edges[edge_begin][edge_end + n].get(GRB_DoubleAttr_X) << ", ";
            //std::cout << edges[edge_begin + n][edge_end].get(GRB_DoubleAttr_X) << " " <<
            //    edges[edge_begin + n][edge_end + n].get(GRB_DoubleAttr_X) << ": ";
            if (edges[edge_begin][edge_end].get(GRB_DoubleAttr_X) +
                edges[edge_begin][edge_end + n].get(GRB_DoubleAttr_X) == 0){
                    optimized_indexes[i - 1] += n;
                    //std::cout << "R" << std::endl;
                }
            //else std::cout << "N" << std::endl;
        }
        size_t edge_begin = optimized_indexes[n - 2];
        if (edge_begin > n) edge_begin -= n;
        size_t edge_end = optimized_indexes[n - 1];
        if (edges[edge_begin][edge_end].get(GRB_DoubleAttr_X) +
            edges[edge_begin + n][edge_end].get(GRB_DoubleAttr_X) == 0){
                optimized_indexes[n - 1] += n;
        }
        

        // Free memory from edges
        for (size_t i = 0; i < n + 1; ++i){
            delete [] edges[i];
        }
        delete [] edges;
        
        // Free memory from indexes
        delete [] indexes;

        return optimized_indexes;
    }
    catch(GRBException e){
        std::cout << e.getErrorCode() << std::endl << e.getMessage() << std::endl;
        return std::vector<size_t>();
    }
}


template <typename kmer_t>
std::vector<size_t> optimize_indexes(const std::vector<kmer_t>& kMers,
        int (*distance)(const std::vector<kmer_t>&, size_t, size_t, size_t),
        size_t k, bool complements = false){
    if (complements){
        return optimize_indexes_with_complements(kMers, distance, k);
    }
    try {
        // n k-mers + 1 s_0 node = n + 1 nodes
        const size_t n = kMers.size();

        // Creating model
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        // Setting model parameters
        model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_OFF); // seems to not help in this case
        model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER); // consumes less memory than concurent, barrier seems to be the best algortihm for this case
        //model.set(GRB_DoubleParam_MemLimit, MEMORY_LIMIT_GB); // prevents the computer from freezing, maybe?
        //model.set(GRB_IntParam_Threads, THREAD_COUNT); // less threads use less memory, barrier seems to never use more than one anyway

        // Edge variables
        GRBVar** edges = new GRBVar*[n + 1];
        for (size_t i = 0; i < n + 1; i++){
            edges[i] = new GRBVar[n + 1];
            for (size_t j = 0; j < n + 1; j++){
                if (i == j) continue;

                std::string name = "Edge_" + std::to_string(i) + "_" + std::to_string(j);
                if (i == n || j == n)
                    edges[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                else
                    edges[i][j] = model.addVar(0.0, 1.0, distance(kMers, k, i, j), GRB_BINARY, name);
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
        
        // Index variables (u)
        GRBVar* indexes = new GRBVar[n + 1];
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

        // Free memory from edges
        for (size_t i = 0; i < n + 1; ++i){
            delete [] edges[i];
        }
        delete [] edges;

        model.optimize();

        std::vector<size_t> optimized_indexes(n);
        for (size_t i = 0; i < n; i++){
            optimized_indexes[round(indexes[i].get(GRB_DoubleAttr_X))] = i;
        }
        
        // Free memory from indexes
        delete [] indexes;

        return optimized_indexes;
    }
    catch(GRBException e){
        std::cout << e.getErrorCode() << std::endl << e.getMessage() << std::endl;
        return std::vector<size_t>();
    }
}

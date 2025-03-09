
#define GUROBI_DIR "/opt/gurobi1103/linux64/include/gurobi_c++.h"
#define MEMORY_LIMIT_GB 16
#define THREAD_COUNT 1

#include GUROBI_DIR
#include <vector>
#include <sstream>
#include <math.h>

#include "distance_functions.h"

template <typename kmer_t>
std::vector<size_t> compute_indexes_ilp(const std::vector<kmer_t>& kMers,
        size_t K, bool complements = false, size_t lower_bound = 0){
    try {
        // n K-mers + 1 s_0 node = n + 1 nodes
        const size_t N = kMers.size();

        // Creating model
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        // Setting model parameters
        // model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_OFF); // seems to not help in this case
        // model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER); // consumes less memory than concurent, barrier seems to be the best algortihm for this case
        //model.set(GRB_DoubleParam_MemLimit, MEMORY_LIMIT_GB); // prevents the computer from freezing, maybe?
        //model.set(GRB_IntParam_Threads, THREAD_COUNT); // less threads use less memory, barrier seems to never use more than one anyway

        // Kmer encoding
        GRBVar** kPos = new GRBVar*[N + 1];
        for (size_t i = 0; i < N + 1; ++i){
            kPos[i] = new GRBVar[2 * K];
            for (size_t j = 0; j < K; ++j){
                std::string name = "KPos_" + std::to_string(i) + "_" + std::to_string(j) + "_f";
                kPos[i][2 * j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                
                name = "KPos_" + std::to_string(i) + "_" + std::to_string(j) + "_s";
                kPos[i][2 * j + 1] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name)
            }
        }

        /*// KMer indexes (u)
        GRBVar* indexes = new GRBVar[N + 1];
        for (size_t i = 0; i < N + 1; i++){
            std::string name = "Index_" + std::to_string(i);
            indexes[i] = model.addVar(0.0, N, 0.0, GRB_INTEGER, name);
        }
        // Index constraints -- indexes ascend
        for (size_t i = 0; i < N; i++){
            for (size_t j = 0; j < N + 1; j++){
                if (i == j) continue;

                GRBLinExpr difference = indexes[i] - indexes[j] + 1 - N + (N * edges[i][j]);
                std::string name = "IndexDiff_" + std::to_string(i) + "_" + std::to_string(j);
                model.addConstr(difference <= 0.0, name);
            }
        }*/

        // Overlaps at positions
        GRBVar** pOv = new GRBVar*[N];
        for (size_t i = 0; i < N; ++i){
            pOv[i] = new GRBVar[K];
            for (size_t j = 0; j < K; ++j){
                std::string name = "pOv_" + std::to_string(i) + "_" + std::to_string(j);
                pOv[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                // pOv[i] means overlap between i - 1 and i, at pos 0 is beginning
            }

            GRBLinExpr ovSpecificity = 0;
            for (size_t j = 0; j < K; ++j){
                ovSpecificity += pOv[i][j];
            }
            std::string name = "Ov_" + std::to_string(i) + "_" + std::to_string(j);
            model.addConstr(ovSpecificity == 1.0, name);
        }

        


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
                    edges[i][j] = model.addVar(0.0, 1.0, trivial_distance(kMers, K, i, j), GRB_BINARY, name);
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

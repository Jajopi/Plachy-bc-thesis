
#define GUROBI_DIR "/opt/gurobi1103/linux64/include/gurobi_c++.h"
#define MEMORY_LIMIT_GB 16
#define THREAD_COUNT 1

#include GUROBI_DIR
#include <vector>
#include <sstream>
#include <math.h>
#include <map>
#include <set>

#include "../kmers.h"

template <typename kmer_t>
class CycleProhibiter : public GRBCallback
{
    size_t K;
    std::vector<kmer_t> const &kMers;
    std::vector<std::vector<GRBVar>> const &in_edges, &out_edges;
    std::vector<GRBVar> const &starting, &ending;
    std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> const &in_overlaps, &out_overlaps;
public:
    CycleProhibiter(size_t k,
        std::vector<kmer_t> const &kmers,
        std::vector<std::vector<GRBVar>> const &in_es,
        std::vector<std::vector<GRBVar>> const &out_es,
        std::vector<GRBVar> const &start,
        std::vector<GRBVar> const &end,
        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> const &in_ovs,
        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> const &out_ovs) :
            K(k), kMers(kmers),
            in_edges(in_es), out_edges(out_es),
            starting(start), ending(end),
            in_overlaps(in_ovs), out_overlaps(out_ovs) {};
protected:
    void callback() {
        try {
            if (where != GRB_CB_MIPSOL) return;
            
            size_t N = in_edges.size();

            std::vector<bool> visited(N, false);

            // Set begin-end path as visited
            for (size_t starting_index = 0; starting_index < N; ++starting_index){
                if (getSolution(starting[starting_index]) != 1) continue;
                
                size_t index = starting_index, length = 0;
                while (getSolution(ending[index]) != 1){
                    // Find outgoing overlap length
                    size_t depth = 0;
                    for (size_t d = 0; d < K; ++d){
                        if (getSolution(out_edges[index][d]) == 1){
                            depth = d;
                            break;
                        }
                    }
                    visited[index] = true;
                    ++length;

                    kmer_t next_overlap = BitSuffix(kMers[index], depth);
                    for (size_t i = 0; i < N; ++i){
                        if (BitPrefix(kMers[i], K, depth) == next_overlap &&
                                getSolution(in_edges[i][depth]) == 1){
                            index = i;
                            if (!visited[index]) break;
                        }
                    }
                }
                std::cerr << "Found path of length " << length << std::endl;
                break; // Only one this path exists
            }

            // Find all cycles and add constraints
            for (size_t starting_index = 0; starting_index < N; ++starting_index){
                if (visited[starting_index]) continue;

                std::vector<std::pair<size_t, size_t>> path; // index, depth

                size_t index = starting_index;
                while (getSolution(ending[index]) != 1){
                    // Find outgoing overlap length
                    size_t depth = 0;
                    for (size_t d = 0; d < K; ++d){
                        if (getSolution(out_edges[index][d]) == 1){
                            depth = d;
                            break;
                        }
                    }
                    path.emplace_back(index, depth);

                    if (visited[index]){ // Cycle detected, add lazy constraint and exit
                        size_t path_start = 0;
                        while (path[path_start].first != index) ++path_start;
                        
                        GRBLinExpr cycle = 0;
                        size_t path_index, length = 0;
                        
                        path_index = path_start;
                        while (path_index < path.size() - 1){
                            cycle += out_edges[path[path_index].first][path[path_index].second];
                            cycle += in_edges[path[path_index + 1].first][path[path_index].second];
                            ++length;
                            ++path_index;
                        }
                        cycle += out_edges[path[path_index].first][path[path_index].second];
                        cycle += in_edges[path[path_start].first][path[path_index].second];
                        ++length;
                        if (length == 1) break;

                        std::cerr << "Adding lazy constraint to cycle of length " << length << std::endl;
                        addLazy(cycle <= 2 * length - 1);
                        break;
                    }
                    visited[index] = true;

                    kmer_t next_overlap = BitSuffix(kMers[index], depth);
                    for (size_t i = 0; i < N; ++i){
                        if (BitPrefix(kMers[i], K, depth) == next_overlap &&
                                getSolution(in_edges[i][depth]) == 1){
                            index = i;
                            if (!visited[index]) break;
                        }
                    }
                }
            }
        } catch (GRBException const & e) {
            std::cerr << "Error during callback: " << e.getErrorCode() << std::endl << e.getMessage() << std::endl;
        } catch (std::exception const & e) {
            std::cerr << "Error during callback:" << e.what() << std::endl;
        }
    }
};

template <typename kmer_t>
std::vector<size_t> compute_indexes(const std::vector<kmer_t>& kMers,
        size_t K, bool complements = false){
    try {
        const size_t N = kMers.size();
        
        // Create model and set parameters

        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_LazyConstraints, 1);

        // Set variables and constraints

        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> in_overlaps; // kmer prefix / suffix , depth
        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> out_overlaps;
        
        std::vector<GRBVar> starting(N), ending(N);
        GRBLinExpr start_sum, end_sum;

        std::vector<std::vector<GRBVar>> in_edges(N);
        std::vector<std::vector<GRBVar>> out_edges(N);
        for (size_t i = 0; i < N; ++i){
            kmer_t kmer = kMers[i];
            
            GRBLinExpr indegree = 0;
            in_edges[i].resize(K);
            for (size_t depth = 0; depth < K; ++depth){
                in_edges[i][depth] = model.addVar(0, 1, 0, GRB_BINARY);
                
                indegree += in_edges[i][depth];
                
                auto key = std::make_tuple(BitPrefix(kmer, K, depth), depth);
                if (out_overlaps.find(key) == out_overlaps.end()){
                    out_overlaps[key] = GRBLinExpr();
                }
                out_overlaps[key] += in_edges[i][depth];
                
                if (depth == 0){
                    starting[i] = model.addVar(0, 1, 0, GRB_BINARY);
                    start_sum += starting[i];

                    indegree += starting[i];
                }
            }
            model.addConstr(indegree == 1.0);
            
            GRBLinExpr outdegree = 0;
            out_edges[i].resize(K);
            for (size_t depth = 0; depth < K; ++depth){
                out_edges[i][depth] = model.addVar(0, 1, K - depth, GRB_BINARY);

                outdegree += out_edges[i][depth];

                auto key = std::make_tuple(BitSuffix(kmer, depth), depth);
                if (in_overlaps.find(key) == in_overlaps.end()){
                    in_overlaps[key] = GRBLinExpr();
                }
                in_overlaps[key] += out_edges[i][depth];

                if (depth == 0){
                    ending[i] = model.addVar(0, 1, 0, GRB_BINARY);
                    end_sum += ending[i];
                    
                    outdegree += ending[i];
                }
            }
            model.addConstr(outdegree == 1);

            /*for (size_t d = 1; d < K; ++d){
                if (BitPrefix(kmer, K, d) == BitSuffix(kmer, d)){
                    GRBLinExpr self_loop = in_edges[i][d] + out_edges[i][d]
                    model.addConstr
                }
            }*/
        }

        model.addConstr(start_sum == 1);
        model.addConstr(end_sum == 1);

        for (const auto &overlap : in_overlaps){
            if (out_overlaps.find(overlap.first) == out_overlaps.end()){
                model.addConstr(overlap.second == 0);
            }
            else {
                const auto& opp_overlap = out_overlaps.find(overlap.first);
                model.addConstr(overlap.second == opp_overlap->second);
            }
        }
        for (const auto &overlap : out_overlaps){
            if (in_overlaps.find(overlap.first) == in_overlaps.end()){
                model.addConstr(overlap.second == 0);
            }
            // else option has already been added
        }

        // Set callback

        CycleProhibiter<kmer_t> cb = CycleProhibiter<kmer_t>(
            K, kMers, in_edges, out_edges, starting, ending, in_overlaps, out_overlaps);
        model.setCallback(&cb);

        // Run optimization

        model.optimize();

        // Reconstruct the path

        // ...

        // std::vector<size_t> optimized_indexes(n);
        // for (size_t i = 0; i < n; i++){
        //     optimized_indexes[round(indexes[i].get(GRB_DoubleAttr_X))] = i;
        // }

        // return optimized_indexes;
    }
    catch(GRBException const & e){
        std::cerr << "Error during optimization:" << e.getErrorCode() << std::endl << e.getMessage() << std::endl;
    }
    catch(std::exception const & e){
        std::cerr << "Error during optimization:" << e.what() << std::endl;
    }
    return std::vector<size_t>();
}

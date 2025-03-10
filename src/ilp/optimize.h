
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
    std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> const &in_overlaps, &out_overlaps;
public:
    CycleProhibiter(size_t k,
        std::vector<kmer_t> const &kmers,
        std::vector<std::vector<GRBVar>> const &in_es,
        std::vector<std::vector<GRBVar>> const &out_es,
        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> const &in_ovs,
        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> const &out_ovs) :
            K(k), kMers(kmers),
            in_edges(in_es), out_edges(out_es),
            in_overlaps(in_ovs), out_overlaps(out_ovs) {};
protected:
    void callback() {
        try {
            size_t N = in_edges.size();

            size_t starting_index = 0;
            for (size_t i = 0; i < N; ++i){
                if (in_edges[i][0].get(GRB_DoubleAttr_X) == 1){
                    starting_index = i;
                    break;
                }
            }

            std::set<size_t> visited;
            std::vector<std::pair<size_t, size_t>> path; // index, depth

            size_t index = starting_index;
            while (out_edges[index][0].get(GRB_DoubleAttr_X) != 1){
                // find outgoing overlap length
                size_t depth = 0;
                for (size_t d = 0; d < K; ++d){
                    if (out_edges[index][d].get(GRB_DoubleAttr_X) == 1){
                        depth = d;
                        break;
                    }
                }
                path.emplace_back(index, depth);

                if (visited.find(index) != visited.end()){ // cycle detected, add lazy constraint and exit
                    size_t path_index = 0;
                    while (path[path_index].first != index) ++path_index;
                    
                    GRBLinExpr cycle;
                    cycle += in_edges[path[path_index].first][path[path.size() - 1].second];
                    while (path_index < path.size() - 1){
                        cycle += out_edges[path[path_index].first][path[path_index].second];
                        cycle += in_edges[path[path_index + 1].first][path[path_index].second];
                    }
                    addLazy(cycle < visited.size());
                    return;
                }
                visited.emplace(index);

                kmer_t next_overlap = BitSuffix(kMers[index], depth);

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
        
        // Create model

        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        // Set variables and constraints

        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> in_overlaps; // kmer prefix / suffix , depth
        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> out_overlaps;

        GRBLinExpr starting, ending;

        std::vector<std::vector<GRBVar>> in_edges(N);
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

                if (depth == 0) starting += in_edges[i][depth];
            }
            model.addConstr(indegree == 1.0);
        }
        std::vector<std::vector<GRBVar>> out_edges(N);
        for (size_t i = 0; i < N; ++i){
            kmer_t kmer = kMers[i];

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

                if (depth == 0) ending += out_edges[i][depth];
            }
            model.addConstr(outdegree == 1);
        }

        model.addConstr(starting == 1);
        model.addConstr(ending == 1);

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
            kMers, in_edges, out_edges, in_overlaps, out_overlaps);
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

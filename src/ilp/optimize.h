
#define GUROBI_DIR "./gurobi/gurobi_c++.h"

#include <vector>
#include <sstream>
#include <math.h>
#include <map>
#include <set>
#include GUROBI_DIR

#include "../kmers.h"

template <typename kmer_t>
class PathFinder : public GRBCallback
{
    size_t K;
    bool COMPLEMENTS;
    std::vector<kmer_t> const &kMers;
    std::vector<std::vector<GRBVar>> const &in_edges, &out_edges;
    std::vector<GRBVar> const &starting, &ending;
    std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> const &in_overlaps, &out_overlaps;
public:
    using Node = std::pair<kmer_t, size_t>;

    size_t complement_index(size_t index){
        if (!COMPLEMENTS) return index;

        if (index < kMers.size() / 2) return index + kMers.size() / 2;
        return index - kMers.size() / 2;
    }

    PathFinder(size_t k, bool complements,
        std::vector<kmer_t> const &kmers,
        std::vector<std::vector<GRBVar>> const &in_es,
        std::vector<std::vector<GRBVar>> const &out_es,
        std::vector<GRBVar> const &start,
        std::vector<GRBVar> const &end,
        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> const &in_ovs,
        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> const &out_ovs) :
            K(k), COMPLEMENTS(complements), kMers(kmers),
            in_edges(in_es), out_edges(out_es),
            starting(start), ending(end),
            in_overlaps(in_ovs), out_overlaps(out_ovs) {};

    std::vector<size_t> get_path();
protected:
    void callback();
};

template <typename kmer_t>
std::vector<size_t> compute_indexes(const std::vector<kmer_t>& kMers,
        size_t K, bool COMPLEMENTS = false){
    try {
        const size_t N = kMers.size();
        
        // Create model and set parameters

        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_LogToConsole, 0);
        env.set(GRB_StringParam_LogFile, "/dev/stderr");
        env.start();

        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_LazyConstraints, 1);

        // Set variables and constraints

        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> in_overlaps; // kmer prefix / suffix , depth
        std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> out_overlaps;
        
        std::vector<GRBVar> starting(N), ending(N);
        GRBLinExpr start_sum, end_sum;

        std::vector<std::vector<GRBVar>> in_edges(N);
        std::vector<std::vector<GRBVar>> out_edges(N);

        std::vector<GRBLinExpr> indegrees(N), outdegrees(N);

        for (size_t i = 0; i < N; ++i){
            kmer_t kmer = kMers[i];
            
            in_edges[i].resize(K);
            for (size_t depth = 0; depth < K; ++depth){
                in_edges[i][depth] = model.addVar(0, 1, 0, GRB_BINARY);
                
                indegrees[i] += in_edges[i][depth];
                
                auto key = std::make_tuple(BitPrefix(kmer, K, depth), depth);
                if (out_overlaps.find(key) == out_overlaps.end()){
                    out_overlaps[key] = GRBLinExpr();
                }
                out_overlaps[key] += in_edges[i][depth];
                
                if (depth == 0){
                    starting[i] = model.addVar(0, 1, 0, GRB_BINARY);
                    start_sum += starting[i];

                    indegrees[i] += starting[i];
                }
            }
            
            out_edges[i].resize(K);
            for (size_t depth = 0; depth < K; ++depth){
                out_edges[i][depth] = model.addVar(0, 1, K - depth, GRB_BINARY);
                
                outdegrees[i] += out_edges[i][depth];
                
                auto key = std::make_tuple(BitSuffix(kmer, depth), depth);
                if (in_overlaps.find(key) == in_overlaps.end()){
                    in_overlaps[key] = GRBLinExpr();
                }
                in_overlaps[key] += out_edges[i][depth];
                
                if (depth == 0){
                    ending[i] = model.addVar(0, 1, 0, GRB_BINARY);
                    end_sum += ending[i];
                    
                    outdegrees[i] += ending[i];
                }
            }
        }
        if (COMPLEMENTS){
            std::vector<GRBVar> this_complement_used(N);
            
            for (size_t i = 0; i < N; ++i){
                this_complement_used[i] = model.addVar(0, 1, 0, GRB_BINARY);
                if (i >= N / 2){
                    model.addConstr(this_complement_used[i] + this_complement_used[i - N / 2] == 1);
                }
                model.addConstr(indegrees[i] == this_complement_used[i]);
                model.addConstr(outdegrees[i] == this_complement_used[i]);
            }
        }
        else {
            for (size_t i = 0; i < N; ++i){
                model.addConstr(indegrees[i] == 1);
                model.addConstr(outdegrees[i] == 1);
            }
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

        PathFinder<kmer_t> pathfinder(
            K, COMPLEMENTS, kMers, in_edges, out_edges, starting, ending, in_overlaps, out_overlaps);
        model.setCallback(&pathfinder);

        // Run optimization

        model.optimize();

        // Reconstruct the path

        std::vector<size_t> ans = pathfinder.get_path();
        std::cerr << ans.size() << std::endl;
        return ans;
    }
    catch(GRBException const & e){
        std::cerr << "Error during optimization:" << e.getErrorCode() << ' ' << e.getMessage() << std::endl;
    }
    catch(std::exception const & e){
        std::cerr << "Error during optimization:" << e.what() << std::endl;
    }
    return std::vector<size_t>();
}

template <typename kmer_t>
inline void PathFinder<kmer_t>::callback(){
    try {
        if (where != GRB_CB_MIPSOL) return;
        // if (COMPLEMENTS) return;
        
        size_t N = in_edges.size();

        std::vector<bool> visited(N, false);
        size_t reachable_count = 0;
        std::vector<size_t> reachable_nodes;

        // Set reachable path as visited
        for (size_t starting_index = 0; starting_index < N / 2; ++starting_index){
            if (getSolution(starting[starting_index]) == 1){
                reachable_nodes.push_back(starting_index);
                visited[starting_index] = true;
                if (COMPLEMENTS) visited[complement_index(starting_index)] = true;
                break;
            }
        }

        while (reachable_nodes.size() > 0){
            size_t index = reachable_nodes.back();
            reachable_nodes.pop_back();

            ++reachable_count;
            if (COMPLEMENTS) ++reachable_count;

            if (getSolution(ending[index]) == 1) continue;
            
            size_t depth = K;
            for (size_t d = 0; d < K; ++d){
                if (getSolution(out_edges[index][d]) == 1){
                    depth = d;
                    break;
                }
            }

            kmer_t next_overlap = BitSuffix(kMers[index], depth);
            for (size_t i = 0; i < N; ++i){
                if (visited[i]) continue;
                if (BitPrefix(kMers[i], K, depth) == next_overlap &&
                        getSolution(in_edges[i][depth]) == 1){

                    reachable_nodes.push_back(i);
                    visited[i] = true;
                    if (COMPLEMENTS) visited[complement_index(i)] = true;
                }
            }
        }
            
        std::cerr << "- Callback: reachable nodes " << reachable_count << " / " << N << std::endl;

        if (reachable_count == N) return;

        // Find all cycles and add constraints
        GRBLinExpr cycles;
        size_t count = 0;

        for (size_t index = 0; index < N; ++index){
            if (visited[index]) continue; // Complement is also visited

            for (size_t d = 0; d < K; ++d){
                if (getSolution(in_edges[index][d]) == 1){
                    cycles += in_edges[index][d];
                    ++count;
                    break;
                }
            }
            for (size_t d = 0; d < K; ++d){
                if (getSolution(out_edges[index][d]) == 1){
                    cycles += out_edges[index][d];
                    ++count;
                    break;
                }
            }
        }

        addLazy(cycles <= count - 1);

    } catch (GRBException const & e) {
        std::cerr << "Error during callback: " << e.getErrorCode() << ' ' << e.getMessage() << std::endl;
    } catch (std::exception const & e) {
        std::cerr << "Error during callback:" << e.what() << std::endl;
    }
}

template <typename kmer_t>
inline std::vector<size_t> PathFinder<kmer_t>::get_path(){
    try {
        size_t N = in_edges.size();

        // Find euler walk in connected graph

        std::map<Node, std::pair<std::vector<std::pair<Node, size_t>>, size_t>> edges;
        std::vector<Node> reached;
        reached.reserve(N);

        std::vector<size_t> next(N), previous(N);
        size_t start_index = N, end_index = N;
        Node start_node, end_node;

        // Construct the graph
        for (size_t i = 0; i < N; ++i){
            if (starting[i].get(GRB_DoubleAttr_X) == 1 && start_index == N) start_index = i;
            if (ending[i].get(GRB_DoubleAttr_X) == 1 && end_index == N){
                end_index = i;
                for (size_t d = 0; d < K; ++d){
                    if (in_edges[i][d].get(GRB_DoubleAttr_X) != 1) continue;

                    end_node = Node(BitPrefix(kMers[i], K, d), d);
                    break;
                }
                continue;
            }

            for (size_t depth = 0; depth < K; ++depth){
                if (out_edges[i][depth].get(GRB_DoubleAttr_X) != 1) continue;
                
                Node out_node = Node(BitSuffix(kMers[i], depth), depth);
                if (i == start_index){
                    start_node = out_node;
                    break; // Starting kmer cannot join any overlap nodes
                }
                
                for (size_t d = 0; d < K; ++d){
                    if (in_edges[i][d].get(GRB_DoubleAttr_X) != 1) continue;

                    edges[Node(BitPrefix(kMers[i], K, d), d)].first.emplace_back(out_node, i);
                    break; // Each kmer has exactly one in-edge
                }
                break; // Each kmer has exactly one out-edge
            }
        }
        edges[end_node].first.emplace_back(Node(0, 0), end_index);

        std::vector<bool> visited(N);

        // Find walk from start to end
        Node node = start_node;
        size_t last_index = start_index;
        while (last_index != end_index){

            if (edges[node].second == 0) reached.push_back(node);

            auto p = edges[node].first[edges[node].second++];
            size_t index = p.second;
            
            visited[index] = true;

            next[last_index] = index;
            previous[index] = last_index;
            
            last_index = index;
            node = p.first;
        }

        // Extend the walk wherever possible
        for (size_t i = 0; i < reached.size(); ++i){
            Node first_node = reached[i];
            if (edges[first_node].second == edges[first_node].first.size()) continue;
            --i;

            size_t edge_index = edges[first_node].second - 1;

            size_t next_index = edges[first_node].first[edge_index].second;
            size_t last_index = previous[next_index];

            node = first_node;
            do {
                if (edges[node].second == 0) reached.push_back(node);

                auto p = edges[node].first[edges[node].second++];
                size_t index = p.second;

                next[last_index] = index;
                previous[index] = last_index;
    
                last_index = index;
                node = p.first;
            } while (node != first_node);

            previous[next_index] = last_index;
            next[last_index] = next_index;
        }
        
        std::vector<size_t> indexes(COMPLEMENTS ? N / 2 : N);
        size_t actual = start_index;
        for (size_t i = 0; i < (COMPLEMENTS ? N / 2 : N); ++i){
            indexes[i] = actual;
            actual = next[actual];
        }
        return indexes;
        
    } catch (GRBException const & e) {
        std::cerr << "Error during path reconstruction: " << e.getErrorCode() << ' ' << e.getMessage() << std::endl;
    } catch (std::exception const & e) {
        std::cerr << "Error during path reconstruction: " << e.what() << std::endl;
    }
    return std::vector<size_t>();
}

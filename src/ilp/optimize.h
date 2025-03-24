
#define GUROBI_DIR "./gurobi/gurobi_c++.h"
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
class PathFinder : public GRBCallback
{
    size_t K;
    std::vector<kmer_t> const &kMers;
    std::vector<std::vector<GRBVar>> const &in_edges, &out_edges;
    std::vector<GRBVar> const &starting, &ending;
    std::map<std::tuple<kmer_t, size_t>, GRBLinExpr> const &in_overlaps, &out_overlaps;
public:
    using Node = std::pair<kmer_t, size_t>;
<<<<<<< HEAD

    size_t complement_index(size_t index){
        if (!COMPLEMENTS) return index;

<<<<<<< Updated upstream
        return (index + kMers.size() / 2) % kMers.size();
=======
        if (index < kMers.size() / 2) return index + kMers.size() / 2;
        return index - kMers.size() / 2;
>>>>>>> Stashed changes
    }

    PathFinder(size_t k, bool complements,
=======
    PathFinder(size_t k,
>>>>>>> parent of be22c1e (Fix ILP complements)
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
    std::vector<size_t> get_path();
protected:
    void callback();
};

template <typename kmer_t>
std::vector<size_t> compute_indexes(const std::vector<kmer_t>& kMers,
        size_t K, bool complements = false){
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
        for (size_t i = 0; i < N; ++i){
            kmer_t kmer = kMers[i];
            
            GRBLinExpr indegree = 0;
            in_edges[i].resize(K);
            for (size_t depth = 0; depth < K; ++depth){
                in_edges[i][depth] = model.addVar(0, 1, 0, GRB_BINARY);

                if (COMPLEMENTS && i >= N / 2){
                    model.addConstr(in_edges[i][depth] == in_edges[i - N / 2][depth]);
                }
                
                indegree += in_edges[i][depth];
                
                auto key = std::make_tuple(BitPrefix(kmer, K, depth), depth);
                if (out_overlaps.find(key) == out_overlaps.end()){
                    out_overlaps[key] = GRBLinExpr();
                }
                out_overlaps[key] += in_edges[i][depth];
                
                if (depth == 0){
                    starting[i] = model.addVar(0, 1, 0, GRB_BINARY);
                    start_sum += starting[i];

<<<<<<< HEAD
                    indegrees[i] += starting[i];
<<<<<<< Updated upstream
=======

                    if (COMPLEMENTS && i >= N / 2){
                        model.addConstr(starting[i] == starting[i - N / 2]);
                    }
>>>>>>> Stashed changes
=======
                    indegree += starting[i];
>>>>>>> parent of be22c1e (Fix ILP complements)
                }
            }
            model.addConstr(indegree == 1.0);
            
            GRBLinExpr outdegree = 0;
            out_edges[i].resize(K);
            for (size_t depth = 0; depth < K; ++depth){
                out_edges[i][depth] = model.addVar(0, 1, K - depth, GRB_BINARY);
<<<<<<< HEAD
<<<<<<< Updated upstream
=======

                if (COMPLEMENTS && i >= N / 2){
                    model.addConstr(out_edges[i][depth] == out_edges[i - N / 2][depth]);
                }
>>>>>>> Stashed changes
                
                outdegrees[i] += out_edges[i][depth];
                
=======

                outdegree += out_edges[i][depth];

>>>>>>> parent of be22c1e (Fix ILP complements)
                auto key = std::make_tuple(BitSuffix(kmer, depth), depth);
                if (in_overlaps.find(key) == in_overlaps.end()){
                    in_overlaps[key] = GRBLinExpr();
                }
                in_overlaps[key] += out_edges[i][depth];

                if (depth == 0){
                    ending[i] = model.addVar(0, 1, 0, GRB_BINARY);
                    end_sum += ending[i];
                    
<<<<<<< HEAD
                    outdegrees[i] += ending[i];
<<<<<<< Updated upstream
                }
            }
=======
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
>>>>>>> parent of be22c1e (Fix ILP complements)
        }

        for (size_t i = 0; i < N; ++i){
            model.addConstr(indegrees[i] == 1);
            model.addConstr(outdegrees[i] == 1);
        }

        if (COMPLEMENTS){            
            for (size_t i = N / 2; i < N; ++i){
                for (size_t depth = 0; depth < K; ++depth){
                    model.addConstr(in_edges[i][depth] == out_edges[i - N / 2][depth]);
                    model.addConstr(out_edges[i][depth] == in_edges[i - N / 2][depth]);
                }    
                model.addConstr(starting[i] == ending[i - N / 2]);
                model.addConstr(ending[i] == starting[i - N / 2]);
            }    
        }    
        
        model.addConstr(start_sum == (COMPLEMENTS ? 2 : 1));
        model.addConstr(end_sum == (COMPLEMENTS ? 2 : 1));
=======

                    if (COMPLEMENTS && i >= N / 2){
                        model.addConstr(ending[i] == ending[i - N / 2]);
                    }
                }
            }
        }
        
        for (size_t i = 0; i < N; ++i){
            model.addConstr(indegrees[i] == 1);
            model.addConstr(outdegrees[i] == 1);
        }

        if (COMPLEMENTS){
            model.addConstr(start_sum == 2);
            model.addConstr(end_sum == 2);
        }
        else{
            model.addConstr(start_sum == 1);
            model.addConstr(end_sum == 1);
        }
>>>>>>> Stashed changes

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
            K, kMers, in_edges, out_edges, starting, ending, in_overlaps, out_overlaps);
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
<<<<<<< HEAD
        if (where != GRB_CB_MIPSOL) return;
<<<<<<< Updated upstream
        // if (COMPLEMENTS) return;
=======
>>>>>>> Stashed changes
=======
        if (where != GRB_CB_MIPSOL) return; // vcelku zbytecne, doporucuji smazat
>>>>>>> parent of be22c1e (Fix ILP complements)
        
        size_t N = in_edges.size();

        std::vector<bool> visited(N, false);
        size_t reachable_count = 0;
        std::vector<size_t> reachable_nodes;
        
        std::vector<std::pair<size_t, size_t>> previous; // index, depth
        if (COMPLEMENTS) previous.resize(N, std::make_pair(N, K));

        // Set reachable path as visited
<<<<<<< HEAD
<<<<<<< Updated upstream
        size_t starting_index = N;
        for (size_t i = 0; i < N / 2; ++i){
            if (getSolution(starting[i]) == 1){
                reachable_nodes.push_back(i);
                visited[i] = true;
                starting_index = i;
                // if (COMPLEMENTS) visited[complement_index(starting_index)] = true;
=======
=======
>>>>>>> parent of be22c1e (Fix ILP complements)
        for (size_t starting_index = 0; starting_index < N; ++starting_index){
            if (getSolution(starting[starting_index]) == 1){
                reachable_nodes.push_back(starting_index);
                visited[starting_index] = true;
<<<<<<< HEAD
                if (COMPLEMENTS) visited[complement_index(starting_index)] = true;
>>>>>>> Stashed changes
=======
>>>>>>> parent of be22c1e (Fix ILP complements)
                break;
            }
        }

        while (reachable_nodes.size() > 0){
            size_t index = reachable_nodes.back();
            reachable_nodes.pop_back();

            ++reachable_count;
<<<<<<< HEAD
            // if (COMPLEMENTS) ++reachable_count;
=======
>>>>>>> parent of be22c1e (Fix ILP complements)

            if (getSolution(ending[index]) == 1) continue;
            
            size_t depth = 0;
            for (size_t d = 0; d < K; ++d){
                if (getSolution(out_edges[index][d]) == 1){
                    depth = d;
                    break;
                }
            }

            kmer_t next_overlap = BitSuffix(kMers[index], depth);
            for (size_t i = 0; i < N; ++i){
<<<<<<< HEAD
<<<<<<< Updated upstream
                if (visited[i]) continue;
=======
>>>>>>> Stashed changes
                if (BitPrefix(kMers[i], K, depth) == next_overlap &&
                        getSolution(in_edges[i][depth]) == 1){

                    // if (reachable_nodes.back() != index) reachable_nodes.push_back(index);
                    
                    if (COMPLEMENTS && visited[complement_index(i)]){
                        // std::cerr << i << ' ' << complement_index(i) << std::endl;
                        // previous[i] = std::make_pair(index, depth);
                        // size_t current = i;
                        // GRBLinExpr bad = 0;
                        // size_t bad_count = 0;
                        // while (current != complement_index(i)){
                        //     auto p = previous[current];
                        //     bad += in_edges[current][p.second];
                        //     bad += out_edges[p.first][p.second];
                        //     bad_count += 2;
                        //     current = p.first;

                        //     if (current == starting_index){
                        //         bad_count = 0;
                        //         break;
                        //     }
                        // }
                        // if (bad_count > 0){
                        //     addLazy(bad <= bad_count - 1);
                        //     return;
                        // }
                        continue;
                    }

                    reachable_nodes.push_back(i);
                    visited[i] = true;
<<<<<<< Updated upstream
                    
                    std::cerr << i << ' '; print_kmer(kMers[i], K, std::cerr, K);
                    std::cerr << ' ' << (i + N / 2) % N << std::endl;

                    if (COMPLEMENTS) previous[i] = std::make_pair(index, depth);
                    // if (COMPLEMENTS) visited[complement_index(i)] = true;

                    break;
=======
                    if (COMPLEMENTS) visited[complement_index(i)] = true;
>>>>>>> Stashed changes
=======
                if (BitPrefix(kMers[i], K, depth) == next_overlap && // Dulezity radek
                        getSolution(in_edges[i][depth]) == 1){
                    if (visited[i]) continue;
                    reachable_nodes.push_back(i);
                    visited[i] = true;
>>>>>>> parent of be22c1e (Fix ILP complements)
                }
            }
        }
            
        std::cerr << "- Callback: reachable nodes " << reachable_count << " / " << N << std::endl;

        if (reachable_count == N || (COMPLEMENTS && reachable_count == N / 2)) return;

        // Find all cycles and add constraints
        GRBLinExpr cycles;

        for (size_t index = 0; index < N; ++index){
            if (visited[index]) continue;

            for (size_t d = 0; d < K; ++d){
                if (getSolution(in_edges[index][d]) == 1){
                    cycles += in_edges[index][d];
                    break;
                }
            }
            for (size_t d = 0; d < K; ++d){
                if (getSolution(out_edges[index][d]) == 1){
                    cycles += out_edges[index][d];
                    break;
                }
            }
        }

        addLazy(cycles <= 2 * (N - reachable_count) - 1);

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
        // (kmer, depth): edges ((kmer, depth), kmer index), first unused edge
        std::vector<Node> reached;
        reached.reserve(N);

        std::vector<size_t> next(N), previous(N);
        size_t start_index = N, end_index = N;
        Node start_node, end_node, complement_end_node;

        // Construct the graph
        for (size_t i = 0; i < N; ++i){
<<<<<<< HEAD
<<<<<<< Updated upstream
            if (starting[i].get(GRB_DoubleAttr_X) == 1 && start_index == N &&
                    i != end_index + N / 2) start_index = i;
            if (ending[i].get(GRB_DoubleAttr_X) == 1 && end_index == N &&
                    i != start_index + N / 2){
=======
            if (starting[i].get(GRB_DoubleAttr_X) == 1 && start_index == N) start_index = i;
            if (ending[i].get(GRB_DoubleAttr_X) == 1 && end_index == N){
>>>>>>> Stashed changes
=======
            if (starting[i].get(GRB_DoubleAttr_X) == 1) start_index = i;
            if (ending[i].get(GRB_DoubleAttr_X) == 1){
>>>>>>> parent of be22c1e (Fix ILP complements)
                end_index = i;
                // std::cerr << "E: " << end_index << std::endl;
                for (size_t d = 0; d < K; ++d){
                    if (in_edges[i][d].get(GRB_DoubleAttr_X) != 1) continue;

                    end_node = Node(BitPrefix(kMers[i], K, d), d);
<<<<<<< HEAD
<<<<<<< Updated upstream
=======
                    if (COMPLEMENTS) complement_end_node = Node(BitPrefix(kMers[complement_index(i)], K, d), d);
                    // edges[Node(BitPrefix(kMers[i], K, d), d)].first.emplace_back(Node(0, 0), i);
>>>>>>> Stashed changes
=======
                    // edges[Node(BitPrefix(kMers[i], K, d), d)].first.emplace_back(Node(0, 0), i);
>>>>>>> parent of be22c1e (Fix ILP complements)
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
        if (COMPLEMENTS){
            edges[complement_end_node].first.emplace_back(Node(0, 0), complement_index(end_index));
        }

        std::vector<bool> complement_visited(COMPLEMENTS ? N : 0);

<<<<<<< HEAD
        std::vector<bool> visited(N);

        std::cerr << "C" << std::endl;

        for (size_t i = 0; i < N; ++i){
            std::cerr << i << ' '; print_kmer(kMers[i], K, std::cerr, K);
            std::cerr << std::endl;
        }
        std::cerr << start_index << ' ' << end_index << std::endl;

=======
>>>>>>> parent of be22c1e (Fix ILP complements)
        // Find walk from start to end
        Node node = start_node;
        size_t last_index = start_index;
<<<<<<< Updated upstream
        while (last_index != end_index){
<<<<<<< HEAD
            std::cerr << last_index << ' '; print_kmer(kMers[last_index], K, std::cerr, K);
            std::cerr << ' ' << (last_index + N / 2) % N;
            std::cerr << ' ' << edges[node].second << ' ' << edges[node].first.size() << ' ' << node.second << std::endl;
=======
        while (last_index != end_index && last_index != complement_index(end_index)){
            // print_kmer(node.first, K, std::cerr, K);
            // std::cerr << ' ' << node.second << std::endl;
>>>>>>> Stashed changes
=======
            // print_kmer(node.first, K, std::cerr, K);
            // std::cerr << ' ' << node.second << std::endl;
>>>>>>> parent of be22c1e (Fix ILP complements)

            if (edges[node].second == 0) reached.push_back(node);

            auto p = edges[node].first[edges[node].second++];
            size_t index = p.second;
<<<<<<< HEAD
            
<<<<<<< Updated upstream
            visited[index] = true;
=======
            if (COMPLEMENTS){
                while (complement_visited[complement_index(index)]){
                    p = edges[node].first[edges[node].second++];
                    index = p.second;
                }
                complement_visited[index] = true;
            }
            // std::cerr << index << std::endl;
>>>>>>> Stashed changes
=======
            // std::cerr << index << std::endl;
>>>>>>> parent of be22c1e (Fix ILP complements)

            next[last_index] = index;
            previous[index] = last_index;
            
            last_index = index;
            node = p.first;
        }

        std::cerr << "E" << std::endl;

        // Extend the walk wherever possible
        for (size_t i = 0; i < reached.size(); ++i){
            Node first_node = reached[i];
            if (edges[first_node].second == edges[first_node].first.size()) continue;
            --i;
            // std::cerr << "F: ";
            // print_kmer(first_node.first, K, std::cerr, K);
            // std::cerr << ' ' << first_node.second << std::endl;

<<<<<<< HEAD
            size_t edge_index = edges[first_node].second - 1;
<<<<<<< Updated upstream

=======
            if (COMPLEMENTS){
                while (complement_visited[complement_index(edges[first_node].first[edge_index].second)]){
                    --edge_index;
                }
            }
>>>>>>> Stashed changes
            size_t next_index = edges[first_node].first[edge_index].second;
=======
            size_t next_index = edges[first_node].first[edges[first_node].second - 1].second;
>>>>>>> parent of be22c1e (Fix ILP complements)
            size_t last_index = previous[next_index];
            // std::cerr << next_index << " <- " << last_index << std::endl;

            node = first_node;
            do {
                // print_kmer(node.first, K, std::cerr, K);
                // std::cerr << ' ' << node.second << std::endl;
                if (edges[node].second == 0) reached.push_back(node);

                auto p = edges[node].first[edges[node].second++];
                size_t index = p.second;
<<<<<<< HEAD

<<<<<<< Updated upstream
=======
                if (COMPLEMENTS){
                    while (complement_visited[complement_index(index)]){
                        p = edges[node].first[edges[node].second++];
                        index = p.second;
                    }
                    complement_visited[index] = true;
                }
    
>>>>>>> Stashed changes
=======
    
>>>>>>> parent of be22c1e (Fix ILP complements)
                next[last_index] = index;
                previous[index] = last_index;
    
                last_index = index;
                node = p.first;
            } while (node != first_node);

            previous[next_index] = last_index;
            next[last_index] = next_index;
        }
        
        std::vector<size_t> indexes(N);
        size_t actual = start_index;
<<<<<<< HEAD
        for (size_t i = 0; i < (COMPLEMENTS ? N / 2 : N); ++i){
<<<<<<< Updated upstream
=======
            // print_kmer(kMers[actual], K, std::cerr, K);
            // std::cerr << '-' << actual << std::endl;
>>>>>>> Stashed changes
=======
        for (size_t i = 0; i < N; ++i){
            // print_kmer(kMers[actual], K, std::cerr, K);
            // std::cerr << '-' << actual << std::endl;
>>>>>>> parent of be22c1e (Fix ILP complements)
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

#ifndef CORE_GADGET_GADGET_H_
#define CORE_GADGET_GADGET_H_

#include <utility>
#include <vector>
#include <time.h>
#include <assert.h>
#include <algorithm>
#include "../defs.h"
#include "../ours-csr-new/def.h"

/*the end and node_idx can be changed by other workers*/
#if 0
// #define GRAPH_EDGE(v)  for (node_t idx = graph.begin[v]; idx < graph.end[v]; idx++){ \
//         node_t u = graph.node_idx[idx];
//read edges from node v as u.  
#define GRAPH_EDGE(v, u) for (node_t idx = (graph.begin[v]); \
        idx < atomic_read(&graph.end[v]); idx++){ \
        node_t u = atomic_read(&graph.node_idx[idx]);
#define GRAPH_EDGE_END }

#else 

#define GRAPH_EDGE(v, u) for (node_t idx = (graph.begin[v]); \
        idx < (graph.end[v]); idx++){ \
        node_t u = (graph.node_idx[idx]);
#define GRAPH_EDGE_END }

#endif



#define clear_cache(x, y) __builtin___clear_cache(x, y)

#if 0
/*local defined CSR format that is used in our experiments. For edges, 
* it has enough sapce for insertion and also support removal*/
class GRAPH {
    typedef int node_t;
private: 
    int max_edge_size;
public:
    /* graph create by csr format 
    * for (node_t idx = G.begin[v];  idx < G.end(v);  idx++)
    *    u = node_idx[idx]; 
    */
    const int n; // number of vertices  
    const int m; // maximum number of edges
    std::vector<node_t> node_idx;
    std::vector<node_t> begin;
    std::vector<node_t> end; // end = next begin for full inserted graph. 
    //std::vector<node_t> mid; // begin to  < mid is all u.post; mid to < end is all u.pre 

    /*all edges is used for allocate memory*/
    GRAPH(int n, int m, std::vector<std::pair<node_t, node_t>> &edges):
            n{n}, m{m}{
        node_idx = std::vector<node_t>(m*2+1); // undirected graph
        begin=std::vector<node_t>(n+1);
        end=std::vector<node_t>(n+1);

        //allocate the memory for undirected graphs
        for (auto const e: edges) { //cnt the edge number for each node
            assert(e.first < n && e.first >=0); 
            assert(e.second < n && e.second >=0);
            end[e.first]++; end[e.second]++; 
        }
        begin[0] = 0;
        for (node_t v = 1; v < n+1; v++) {
            begin[v] = begin[v-1] + end[v-1]; //accumulate the size
        }
        // initially the graph is empty
        for (node_t v = 0; v < n+1; v++) {
            end[v]  = begin[v];
            //mid[v] = begin[v];
        }
        // for(node_t v=0; v<n; v++) {
        //     ASSERT(begin[v] <= begin[v+1]);
        //     ASSERT(begin[v] <= 2*m);
        // }

        max_edge_size = 0;
        for(node_t v = 0; v < n+1; v++) {
            int edge_size = begin[v+1]- begin[v] + 1;
            if (edge_size > max_edge_size) {
                max_edge_size = edge_size; 
            }
        }
    }

    int get_max_edge_size() {return max_edge_size; }
    
    // for undirected graph, both u and v are locked. 
    void insert(node_t u, node_t v) { 
#if 0
        node_idx[end[u]] = v;  clear_cache(&node_idx[end[u]], &node_idx[end[u]] + 1);

        end[u]++; clear_cache(&end[u], &end[u] + 1);

        //assert(end[u] <= begin[u+1]); 
        node_idx[end[v]] = u;  clear_cache(&node_idx[end[v]], &node_idx[end[v]] + 1);

        end[v]++; clear_cache(&end[v], &end[v] + 1);
        //assert(end[v] <= begin[v+1]);

#else
    node_t endu = atomic_read(&end[u]);
    node_idx[endu] = v;
    atomic_write(&end[u], endu+1);

    assert(begin[u] <= end[u]);
    assert(end[u] <= begin[u+1]); 

    node_t endv = atomic_read(&end[v]);
    node_idx[endv] = u;  
    atomic_write(&end[v], endv+1); 

    assert(begin[v] <= end[v]);
    assert(end[v] <= begin[v+1]);
#endif 

    }
    
    // for undirected graph
    void remove(node_t u, node_t v) {

        // for (int idx = begin[u]; idx < end[u]; idx++) {
        //     if (v == node_idx[idx]) { // edges are not empty
        //         node_t pos = end[u] - 1;
        //         std::swap(node_idx[idx], node_idx[pos]);

        //         clear_cache(&node_idx[idx], &node_idx[idx] + 1);
        //         clear_cache(&node_idx[pos], &node_idx[pos] + 1);

        //         --end[u]; clear_cache(&end[u], &end[u] + 1);
        //         break;
        //     }
        // }
        
        // for (int idx = begin[v]; idx < end[v]; idx++) {
        //     if (u == node_idx[idx]) { // edges are not empty
        //         node_t pos = end[v] - 1;
        //         std::swap(node_idx[idx], node_idx[pos]);

        //         clear_cache(&node_idx[idx], &node_idx[idx] + 1);
        //         clear_cache(&node_idx[pos], &node_idx[pos] + 1);

        //         --end[v]; clear_cache(&end[v], &end[v] + 1);
        //         break;
        //     }
        // }
#if 0
    node_t endu = atomic_read(&end[u]);
    for (int idx = begin[u]; idx < endu; idx++) {
        if (v == atomic_read(&node_idx[idx])) { // edges are not empty
            node_t pos = endu - 1;

            node_t w = atomic_read(&node_idx[idx]);
            atomic_write(&node_idx[idx], atomic_read(&node_idx[pos]));
            atomic_write(&node_idx[pos], w);
            //std::swap(node_idx[idx], node_idx[pos]);

            atomic_write(&end[u], endu - 1);
            break;
        }
    }
    
    node_t endv = atomic_read(&end[v]);
    for (int idx = begin[v]; idx < end[v]; idx++) {
        if (u == atomic_read(&node_idx[idx])) { // edges are not empty
            node_t pos = endv - 1;

            node_t w = atomic_read(&node_idx[idx]);
            atomic_write(&node_idx[idx], atomic_read(&node_idx[pos]));
            atomic_write(&node_idx[pos], w);
            //std::swap(node_idx[idx], node_idx[pos]);

            atomic_write(&end[v], endv - 1);
            break;
        }
    }

#else 
    
    // use std::find for speedup
    node_t endu = atomic_read(&end[u]);
    int *idx = std::find(&node_idx[begin[u]], &node_idx[endu], v); 
    if (idx != &node_idx[endu]) { // find it
        node_t pos = endu - 1;
        atomic_write(idx, atomic_read(&node_idx[pos]));
        atomic_write(&end[u], pos);
    } 
    
    node_t endv = atomic_read(&end[v]);
    idx = std::find(&node_idx[begin[v]], &node_idx[endv], u);
    if (idx != &node_idx[endv]) { //find it 
        node_t pos = endv - 1;
        atomic_write(idx, atomic_read(&node_idx[pos]));
        atomic_write(&end[v], pos);    
    }

#endif 
    }
    // the size of edge for node u; 
    int size(node_t u) {return atomic_read(&end[u]) - atomic_read(&begin[u]);}

    
};



#else  // local GRAPH 

/*local defined CSR format that is used in our experiments. For edges, 
* it has enough sapce for insertion and also support removal*/
class GRAPH {
    typedef int node_t;
public:
    /* graph create by csr format 
    * for (node_t idx = G.begin[v];  idx < G.end(v);  idx++)
    *    u = node_idx[idx]; 
    */
    const int n; // number of vertices  
    const int m; // maximum number of edges
    std::vector<node_t> node_idx;
    std::vector<node_t> begin;
    std::vector<node_t> end; // end = next begin for full inserted graph. 
    //std::vector<node_t> mid; // begin to  < mid is all u.post; mid to < end is all u.pre 
    
    
    /*all edges is used for allocate memory*/
    GRAPH(int n, int m, std::vector<std::pair<node_t, node_t>> &edges):
            n{n}, m{m}{
        node_idx = std::vector<node_t>(m*2+1); // undirected graph
        begin=std::vector<node_t>(n+1);
        end=std::vector<node_t>(n+1);

        //allocate the memory for undirected graphs
        for (auto const e: edges) { //cnt the edge number for each node
            assert(e.first < n && e.first >=0); 
            assert(e.second < n && e.second >=0);
            end[e.first]++; end[e.second]++; 
        }
        begin[0] = 0;
        for (node_t v = 1; v < n+1; v++) {
            begin[v] = begin[v-1] + end[v-1]; //accumulate the size
        }
        // initially the graph is empty
        for (node_t v = 0; v < n+1; v++) {
            end[v]  = begin[v];
            //mid[v] = begin[v];
        }
        // for(node_t v=0; v<n; v++) {
        //     ASSERT(begin[v] <= begin[v+1]);
        //     ASSERT(begin[v] <= 2*m);
        // }
    }

    // for undirected graph
    void insert(node_t u, node_t v) { 
        node_idx[end[u]] = v; end[u]++;
        assert(end[u] <= begin[u+1]); 
        node_idx[end[v]] = u; end[v]++;
        assert(end[v] <= begin[v+1]);
    }
    
    // for undirected graph
    void remove(node_t u, node_t v) {
        //node_t endu = end[u];
        int *idx = std::find(&node_idx[begin[u]], &node_idx[end[u]], v);
        if (idx != &node_idx[end[u]]) { // find it
            node_t pos = end[u] - 1;
            *idx = (node_idx[pos]);
            end[u] = pos;
        } 

        //node_t endv = end[v];
        idx = std::find(&node_idx[begin[v]], &node_idx[end[v]], u);
        if (idx != &node_idx[end[v]]) { //find it 
            node_t pos = end[v] - 1;
            *idx = (node_idx[pos]);
            end[v] = pos;    
        }


        // for (int idx = begin[u]; idx < end[u]; idx++) {
        //     if (v == node_idx[idx]) { // edges are not empty
        //         std::swap(node_idx[idx], node_idx[end[u] - 1]);
        //         --end[u];
        //         break;
        //     }
        // }
        // for (int idx = begin[v]; idx < end[v]; idx++) {
        //     if (u == node_idx[idx]) { // edges are not empty
        //         std::swap(node_idx[idx], node_idx[end[v] - 1]);
        //         --end[v];
        //         break;
        //     }
        // }
    }
    // the size of edge for node u; 
    int size(node_t u) {return end[u] - begin[u];}

    inline node_t get(int id) {return node_idx[id];}
    
};

#endif // is or not parallel graph 


namespace gadget {

void RepeatWith(const char symbol, const int repeat);
std::vector<std::vector<int>> ReadGraph(const char* const path,
                                        int* const n, int* const m);
std::vector<std::pair<int, int>> ReadTempEdgesS(const char* const path,
                                                int* const n, int* const m);
std::vector<std::pair<int, int>> ReadEdgesS(const char* const path,
                                            int* const n, int* const m);

std::vector<std::pair<int, int>> ReadEdges(const char* const path,
                                            int* const n, int* const m);
                                            
std::vector<std::pair<int, int>> CSRReadEdges(char* const path,
                    int* const n, int* const m, int shuffle); 

std::vector<std::pair<int, int>> ReadTemporalEdges(char* const path,
                    int* const n, int* const m, int shuffle = 0); 


/*read m edges from the beggining and stored as our new GRAPH*/
GRAPH CreateOurCSRGRAPH(std::vector<std::pair<int, int>> &edges, int n, int size); 
//PARGRAPH CreateOurCSRParGRAPH(std::vector<std::pair<int, int>> &edges, int n, int size);

void OutputCoreNumber(const char* const path, std::vector<int> &core, int n);


void OutputSampleEdgeCoreNumber(const char* const path, std::vector<int> &core, 
    std::vector<std::pair<int, int>> &edges, const int n);

std::vector<std::pair<int, int>> sampleEdges(const char* const path, std::vector<std::pair<int, int>> &edges, const float percent);

void CutEdges(std::vector<std::pair<int, int>> &edges, const int num);

std::vector<int> RepeateRandomRead(const char* path);

}  // namespace gadget

#endif

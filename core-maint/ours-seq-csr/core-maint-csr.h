#pragma once
#include <vector>
#include <set>
#include <queue>
#include <math.h>       /* pow */
#include <climits>
#include <map>
#include "gadget/gadget.h"
using namespace std;


//for graph
const size_t MAX_INCREASE_CORE = 1024*16; //max increased core number, should be large enough
const size_t VECTOR_CAPACITY = 1024*1024;
typedef int node_t; // the vertex of graph
const node_t NONE = -1; 
typedef int edge_t; // the edge of graph
typedef int core_t; // the core number of graph 
typedef int deg_t;  // degree
typedef int color_t;
typedef int label_t;
typedef int worker_t;
typedef unsigned int counter_t;
//typedef vector<vector<edge_t>> graph_t; 

typedef unsigned long long int tag_t; // 64bit label tag;
typedef unsigned int subtag_t; // 32bit label subtag;
const size_t MAX_TAG = 0xffffffffffffffff;
const size_t INIT_TAG_GAP = 0xffffffff; // 32 bit unsigned integer
const size_t MAX_SUBTAG = 0xffffffff; // 32 bit unsigned integer


/**priority queue for core maint begin*************/
struct DATA {
    node_t key;
    tag_t value;
    DATA(node_t key, tag_t value) : key(key), value(value){ } 
};

struct CompareOrder { 
    bool operator()(DATA const& a, DATA const& b) 
    { 
        // return "true" if "a" is ordered before "b", for example:
        return a.value > b.value;
    } 
}; 


class PRIORITY_Q {
private:
    typedef std::priority_queue<DATA, std::vector<DATA>, CompareOrder> PQ;
    PQ pq_;

public:
    PRIORITY_Q(size_t size) {
        std::vector<DATA> container;
        container.reserve(size);
        CompareOrder compare;
        pq_ = PQ(compare, std::move(container));
    } 
    PRIORITY_Q(){}
    inline void push(DATA d) { pq_.push(d);}
    inline DATA top() {return pq_.top();}
    inline void pop() {pq_.pop();}
    inline bool empty() {return pq_.empty(); }
    inline void clear() { while (!pq_.empty()){pq_.pop();} }
};
/**priority queue for core maint end *************/

/**queue for core maint begin*********************/
// this queue can be vector
class QUEUE {
private:
    std::vector<node_t> container;
public:
    QUEUE(size_t size) {
        container.reserve(size);
    }
    QUEUE() {}
    inline void push(node_t v) {container.push_back(v);}
    inline node_t top() {return container.back();}
    inline void pop() {container.pop_back();}
    inline bool empty() {return container.empty(); }
    inline void clear() { container.clear();}
};
/***queue for core maint end*********************/


enum{
    WHITE   = 0, /*initial*/
    GRAY    = 1, /*visited, not in V* but in V+*/
    BLACK   = 2, /*visited, in V* and in V+*/
};

namespace SeqCM{
    class Node {
    public:
        //Computer Core number
        deg_t degout; // deg+
        deg_t degin;  // deg*
        deg_t mcd;    //for edge removing.
       
        //Order Maintenance
        node_t pre;  //double linked list
        node_t next; //double linked list
        tag_t tag;   // label tag for each item range is 0 - n^2
        color_t color; // black, white, gray.
        label_t inQ;   // node is in priority queue PQ
        label_t inR;    // node is in R can be defined with color
        //label_t outR;   // node is already backward.
        Node(){
            degout = degin = mcd = tag = 0; // init to 0 
            pre = next = NONE;
            color = WHITE;
            inQ = false;
            inR = false; 
        }
    };

    /*first version, with one level tags (label), only tag is used.
    * take O(logn) time for insert*/
    class CoreMaint {
        private:
        size_t n; // the number of vertices in graph.
        core_t max_core; // the max core number
        GRAPH &graph;
        vector<core_t>& core;
        vector<Node> V; // node information include head and tail
        // head-> x -> y -> tail (easy for insert and remove)
        
        /*maintaining the order   head(avoid empty list check)->node->node(tail)*/
        vector<core_t> H;  //each core list has a head to avoid empty list  
        vector<core_t> T;  //each core list has a tail;
        
        /*list operation is here*/
        inline void ListDelete(node_t x);
        inline void ListInsert(node_t x, node_t y);
        inline void MultiListInsert(node_t x, vector<node_t> &y);
        //inline void ListLink(vector<node_t> &y); //link all nodes in x.

        /*order operation*/
        inline bool Order(node_t x, node_t y); // return x < y
        inline bool SameCoreOrder(node_t x, node_t y);
        inline bool TagOrder(node_t x, node_t y) {return V[x].tag < V[y].tag;}
        int OrderInsert(node_t x, node_t y); // insert y after x
        int MultiOrderInsert(node_t x, vector<node_t> &y);
        inline int OrderDelete(node_t x);   // remove x

        QUEUE R; 
        PRIORITY_Q PQ;      // queue backward
        vector<node_t> Vblack;  // the order of black nodes (may include gray)
        vector<node_t> Vcolor; // all colored vertices.
        
        public:
        CoreMaint(size_t n, GRAPH &graph, vector<core_t> &core);
        ~CoreMaint(){}
        void Init(vector<int> &odv);
        void ComputeCore(vector<int> &core, vector<int> &order_v, bool with_init);
        
        // for edge insert
        inline void Forward(node_t w);
        inline void Backward(node_t w);
        inline void DoPre(node_t u);
        inline void DoAdj(node_t u);
        int EdgeInsert(node_t x, edge_t y); // insert to double-linked list
        
        // for edge remove
        inline void BeginRemove(node_t u, node_t v);
        inline void BeginRemoveContinue(node_t u, node_t v);
        inline void FindVs(const core_t K);
        inline void AfterRemove(const core_t K);
        int EdgeRemove(node_t x, edge_t y); //remove in double-linked list.      

        // check the correctness of algorithm
        int Check(node_t, node_t, int id, vector<core_t> &tmp_core, vector<node_t> &order_v); 
    };
}




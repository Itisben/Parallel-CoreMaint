/***********************
 * This is a sequential version of our Parallel Core Maintenance
 * 1) two level of OM data structure. 
 * 2) new removal algorithms
 * 3) After finish this, we extend to Parallel Version without difficulties. 
 * 4) Using CSR format to load large graphs. 
 **************/

#pragma once
#include <vector>
#include <set>
#include <queue>
#include <math.h>       /* pow */
#include <climits>
#include <map>
#include <omp.h>
#include "gadget/gadget.h"
#include "def.h"
using namespace std;


namespace ParCM{

    enum{
        WHITE   = 0, /*initial*/
        GRAY    = 1, /*visited, not in V* but in V+*/
        BLACK   = 2, /*visited, in V* and in V+*/
    };

    /*the node can be changed by other worker*/
    class Node {
    public:
        //Computer Core number
            deg_t degout; // deg+
            deg_t degin;  // deg*
        
       
        //Order Maintenance
            node_t pre, next;  //double linked list
            node_t group; //link to group id as parent.  
            subtag_t subtag; // sub label tag 32 bits. 
       
            color_t color;  // black, white, gray.
            //bool inQ;      
            bool inR;     

        // for edge insert
            int islock;     // >= 0 is locked by worker p. NONE not locked.

        // for edge remove 
        int s;         // status shared by all workers
        deg_t mcd;    //for edge removing.
        int s_lock;  // set core-- and s= 2 atomically

        omp_lock_t omlock_omp; 
        omp_lock_t lock_omp;
        lock_t omlock_cas;
        lock_t lock_cas;

        /*buffer queueu*/
        int inB; // 0 not, 1 at the head of buffer queue, 2 at the tail.  

        Node();
        ~Node();
        // for OM data structure
        inline void OMLock();
        inline void OMUnlock();
        inline bool OMTestLock();

        // for Core Maint
        inline void Lock();
        inline void Unlock();
        inline void LockP(int p);
        inline bool TestLockP(int p);
        inline void UnlockP();
        inline bool TestLock();
        inline bool IsLocked();
        //inline bool TestCoreLock(core_t k, core_t & current_k);
        inline bool LockWithCore(const core_t k, core_t *core);
    };

    class Group {
    public:
            size_t size;  // the number of children. Group is removed if 0.
            tag_t tag;   // label tag 64 bits for each item range is 0 - n^2
            node_t pre, next; // double linked list
            
        
        omp_lock_t omlock_omp; 
        lock_t omlock_cas;

        Group();
        ~Group();
        inline void OMLock();
        inline void OMUnlock();
        //inline bool OMTestLock(core_t k, core_t & current_k);
    };

      /**priority queue for core maint begin*************/
    struct DATA {
        node_t id; tag_t tag; subtag_t subtag; int status; id_t version;
        DATA(node_t id, tag_t tag, subtag_t subtag, int status, id_t version) : 
                id(id), tag(tag), subtag(subtag), status(status), version(version){ } 
    };

    struct CompareOrder { 
        bool operator()(DATA const& a, DATA const& b) const { 
            // return "true" if "a" is ordered before "b", for example:
            if (a.tag > b.tag) {return true; }
            else if( a.tag == b.tag) {
                if (a.subtag > b.subtag) return true;
            }
            return false;
        } 
    }; 
    /*
    * we use heap for priority queue to iterate all elements 
    * in case some of them changed. */

    class PRIORITY_Q {
    public:
        vector<DATA> v;
        id_t version;
    public:
        PRIORITY_Q(size_t size) {
            v.reserve(size);
            version = 0;
        }
        PRIORITY_Q(){}
        /*if the version is changed, all nodes needs upate until version is same
        *all the vesion number is same 
        * Return False: the version is changed. all elements needs to update to same version.*/
        inline bool push(DATA &d) { 
            v.push_back(d); push_heap(v.begin(), v.end(), CompareOrder()); return true;}
        inline DATA front() {return v.front();}
        inline void pop() {pop_heap(v.begin(), v.end(),CompareOrder()); v.pop_back();}
        inline bool empty() {return v.empty(); }
        inline void clear() { v.clear(); version = 0;}
        inline void init() {make_heap(v.begin(), v.end(), CompareOrder());}
        
        /*concurrent version of operations*/
        inline void enqueue(id_t *ver, DATA &d);
        inline node_t dequeue(id_t *ver, id_t *cnt, core_t K, vector<Node> &V, vector<Group> &G, 
                vector<int> &Qdegin, vector<bool> &Qin, vector<core_t> &core, int p);
        inline void update_version(id_t *ver, id_t *cnt, vector<Node> &V, vector<Group> &G);
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

    class Count {
        public:
        size_t tag; // top-label update cnt
        size_t subtag; // bottom-lable update cnt
        size_t order_repeat; // Order operation repeat times. 
        size_t relabel; // relabel numbers
        vector<int> Vcount; 

        Count(){tag = subtag = order_repeat = relabel = 0;
                Vcount = vector<int>(10000, 0); }

        void Clear(){tag = subtag = order_repeat = relabel = 0;}
        
        Count operator+(const Count& c) {
            Count cnt;
            cnt.tag = this->tag + c.tag;
            cnt.subtag = this->subtag + c.subtag;
            cnt.order_repeat = this->order_repeat + c.order_repeat;
            cnt.relabel = this->relabel + c.relabel;

            for (int i = 0; i < Vcount.size(); i++) {
                cnt.Vcount[i] = this->Vcount[i] + c.Vcount[i];
            }
            return cnt;
        }
    };


    /*first version, with one level tags (label), only tag is used.
    * take O(logn) time for insert*/
    class CoreMaint {
        private:
        const size_t n; // the number of vertices in graph.
        core_t max_core; // the max core number
        GRAPH &graph;
        vector<core_t> &core;
        vector<Node> V; // node information include head and tail
        
        vector< Group> G; // node in group for two level OM. 
        size_t group_size; // number of 

        // head-> x -> y -> tail (easy for insert and remove)
        
        /*maintaining the order   head(avoid empty list check)->node->node(tail)*/
        vector<core_t> H;  //each core list has a head to avoid empty list  
        vector<core_t> T;  //each core list has a tail;
        vector<core_t> Hg;  //top group head. 
        vector<core_t> Tg;  //top group tail.  

        vector<subtag_t> assigned_subtags;

        /*buffer queue*/
        vector<node_t> B; 
        vector<int> Btail;
        vector<int> Bhead;

        /*double linked list operations*/
        inline void ListDelete(node_t x);
        inline void ListInsert(node_t x, node_t y);
        inline void MultiListInsert(node_t x, vector<node_t> &y, size_t size);
        //inline void ListLink(vector<node_t> &y); //link all nodes in x.
        inline void TopListDelete(group_t x);
        inline void TopListInsert(group_t x, group_t y);
        inline void MultiTopListInsert(group_t x, vector<group_t> &y);

        /*OM tag*/
        inline tag_t getTag(node_t u) {return G[V[u].group].tag;}
        inline subtag_t getSubtag(node_t u) { return V[u].subtag; }
        inline node_t getGroup(node_t u){return V[u].group;}

        /*OM data structure, two level version*/
        inline bool Order(node_t x, node_t y); // return x < y
        inline bool SameCoreOrder(node_t x, node_t y);
        void OrderInsert(node_t x, node_t y, vector<node_t> &groups, int p); // insert y after x
        void SimpleRelabel2(node_t x, node_t xnext, group_t gy, vector<node_t> &groups, int p);
        void SimpleRelabel(node_t x); // for easy implementation
        void DeepRelabel(node_t x); // fully deep relabel
        void AssignLabel();
        void MultiOrderInsert(node_t x, vector<node_t> &y);
        bool OrderDelete(node_t x);   // remove x

        /* inserting to head or tail of list. Can do without relabel.*/
        void OrderInsertHeadBatch(core_t k, vector<node_t> &y, vector<node_t> &groups);
        void OrderInsertTail(core_t k, node_t y);
        void OrderInsertTailBatch(core_t k, vector<node_t> &y, vector<node_t> &groups); 


        /*buffer queue for k-order*/
        void OrderInsertBufferHeadBatch(core_t k, vector<node_t> &y, vector<node_t> &groups);
        void OrderInsertBufferFlush(core_t k,  vector<node_t> &B, vector<node_t> &groups);
        void OrderRemoveBufferTailBatch();
        void OrderRemoveBufferFlush();

        /* current implementation.  
        * the tag version, used by PQ. all node in PQ has the same version
        * odd: lock, even unlock. use the memory transaction. 
        * */ 
        vector<id_t> om_ver; // the version for each order list.
        vector<id_t> om_cnt; // the number of access workers for each order list. 
        PRIORITY_Q PQ;      // queue backward
        void PQ_update(const core_t k);   // update all nodes with new version.

       
        omp_lock_t omlock_omp; 
        lock_t omlock_cas;


        public:
        
        /*counter*/
        vector<Count> cnt;

        CoreMaint(size_t n, GRAPH &graph, vector< core_t> &core);
        ~CoreMaint();
        // global lock
        inline void OMLock();
        inline void OMUnlock();

        void Init(vector<int> &odv);
        void InitParallel(int num_worker);
        core_t ComputeCore(vector<int> &core, vector<int> &order_v, bool with_init); // return max_core
        inline int getMaxCore(){return max_core;}

        // for edge insert
        inline void Forward(node_t w, PRIORITY_Q &PQ, vector<node_t> &Vcolor, 
                vector<int> &Qdegin, vector<bool> &Qin, int p);
        inline void Backward(node_t w, QUEUE &R, vector<int> &Qdegin, vector<bool> &Qin, 
                vector<node_t> &groups, int p);
        inline void DoPre(node_t u, QUEUE &R, int p);
        inline void DoAdj(node_t u, QUEUE &R, vector<int> &Qdegin, vector<bool> &Qin, int p);
        int EdgeInsert(node_t x, edge_t y,PRIORITY_Q &PQ, QUEUE &R,
                vector<node_t> &Vblack, vector<node_t> &Vcolor, vector<int> &Qdegin, 
                vector<bool> &Qin, vector<node_t> &groups, int p); // insert to double-linked list
        int ParallelEdgeInsert(vector<pair<int, int>> &edges, int m, int m2, bool ispar);
        
        // for edge remove
        void DoMCD(node_t u, core_t k, QUEUE &R, vector<node_t> &Vblack);
        void CheckMCD(node_t u, core_t k, node_t from_v);
        int EdgeRemove(node_t x, edge_t y, QUEUE &R, vector<node_t> &Vblack, 
                vector<bool> &A, vector<node_t> &groups, int p); //remove in double-linked list.      
        int ParallelEdgeRemove(vector<pair<int, int>> &edges, int m, int m2, bool ispar);
        
        /* test for debug 
        *  */
        bool CheckLock(bool info); // test all nodes unlock or not.
        vector<node_t> InitTestOM(size_t insert_size);            
        int TestOM(vector<node_t> &nodes); 
        int CheckOM(bool info);
        int CheckAllMCDValue(bool info);

        /* check the correctness of algorithm */
        int CheckCore(vector<core_t> &tmp_core, bool info); //check core numbers.
        int CheckDeg(bool info);
        int CheckAll(node_t, node_t, int id, vector<core_t> &tmp_core, vector<node_t> &order_v); 
        void PrintOMVersion();


    };
}

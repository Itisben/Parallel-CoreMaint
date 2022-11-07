#pragma once
#include <vector>
#include "def.h"


namespace ParOM {
    /*test cases*/
    enum TestCase{
        NO_RELABEL,
        FEW_RELABEL, // average case: randomly pick positions
        MANY_RELABEL, // worse case: pick like 100 fixed positions.
        MAX_RELABEL, // worst case: pick single one positon like after head. 
    };
    
    const int OM_SIZE           = 10000000;  // 10m  
    const int OM_INSERT_SIZE    = 10000000;  // insert to 10 m random positions.
    const int OM_POS_NO         = 10000000;  // insert to 1 m random position. 
    const int OM_POS_FEW        = 1000000;     // insert to 1000 random position. 
    const int OM_POS_MANY      =  1000;
    const int OM_POS_MAX       = 1;
    
    
    /*the node can be changed by other worker*/
    class Node {
    public:
        volatile node_t pre, next;  //double linked list
        volatile node_t group;      //link to group id as parent.  
        volatile subtag_t subtag;   //sub label tag 32 bits. 
        volatile bool live;         //false is deleted true is not deleted

        omp_lock_t omlock_omp;
        omp_lock_t lock_omp;
        volatile lock_t omlock_cas;
        volatile lock_t lock_cas;

        Node();
        ~Node();
        inline void OMLock();
        inline void OMUnlock();
        int OMTestLock();
    };

    class Group {
    public:
        volatile size_t size;  // the number of children. Group is removed if 0.
        volatile tag_t tag;   // label tag 64 bits for each item range is 0 - n^2
        volatile node_t pre, next; // double linked list
        
        omp_lock_t omlock_omp; 
        volatile lock_t omlock_cas;

        Group();
        ~Group();
        inline void OMLock();
        inline void OMUnlock();
        int OMTestLock();

    };
    
    class Count {
        public:
        size_t tag; // top-label update cnt
        size_t subtag; // bottom-lable update cnt
        size_t order_repeat; // Order operation repeat times. 
        size_t relabel;
        Count(){tag = subtag = order_repeat = relabel = 0;}
        Clear(){tag = subtag = order_repeat = relabel = 0;}
        Count operator+(const Count& c) {
            Count cnt;
            cnt.tag = this->tag + c.tag;
            cnt.subtag = this->subtag + c.subtag;
            cnt.order_repeat = this->order_repeat + c.order_repeat;
            cnt.relabel = this->relabel + c.relabel;
            return cnt;
        }
    };

    class OM {
    private:
        size_t max_size; // the max number of items
        size_t node_size;     // number of nodes.
        vector<Node> V;
        vector<Group> G; // node in group for two level OM. 
        volatile size_t group_size; // number of groups. 
        
        node_t H; // head of bottom list
        node_t T; // tail 
        group_t Hg; // head of top group list
        group_t Tg; // tail of top group list

        vector<subtag_t> assigned_subtags;

        // lock when alocate group in relabel process
        omp_lock_t omlock_omp; 
        lock_t omlock_cas;

        Count g_cnt; // global counter

        void Lock();
        void Unlock();

        /*double linked list operations*/
        inline void ListDelete(node_t x);
        inline void ListInsert(node_t x, node_t y);
        inline void MultiListInsert(node_t x, vector<node_t> &y);

        inline void TopListDelete(group_t x);
        inline void TopListInsert(group_t x, group_t y);
        inline void MultiTopListInsert(group_t x, vector<group_t> &y);

    public:

        // n is init size, maxn is the maximum size after insert
        OM(size_t max_size);
        ~OM();
        void Init(vector<node_t> &nodes);
        void InitParallel(int num_worker);


        /*OM operations*/
        inline bool Order(node_t x, node_t y, ParOM::Count &cnt); // return x < y
        void Insert(node_t x, node_t y, vector<node_t> &relabel_nodes, 
                vector<group_t> &groups, ParOM::Count &cnt); // insert y after x
        void SimpleRelabel2(node_t x, node_t xnext, group_t gy, 
                vector<node_t> &relabel_nodes,vector<group_t> &groups, ParOM::Count &cnt);
        inline bool Delete(node_t x);   // remove x
        inline bool DeleteByFlag(node_t x);
    
        vector<int> GetRepeatRandomNum(char * const path);
        vector<int> GenerateTestCase(TestCase t, vector<int> &rand, size_t size);
        vector<int> GenerateTestCase2(TestCase t, size_t insertNum);
        vector<node_t> AlocateNodes(int num);
        void TestInsert(vector<int> &pos, vector<int> &nodes, bool ispar);
        void TestOrder(vector<int> &pos, bool ispar);
        void TestDelete(vector<int> &pos, bool ispar);
        void TestMixed(vector<int> &pos, vector<int> &nodes, bool ispar);
        void PrintCount(char* str);
    };
}
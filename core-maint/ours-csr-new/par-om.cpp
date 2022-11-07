/**************************************
 * this is to test the parallel OM data structure
 * for the paper "paralle Order Maintenance Data Structure"
 * this is one part of work for our parallel Core Maint.
 * **********************************/
#include "def.h"
#include "par-om.h"
#include "../gadget/gadget.h"
#include <string>
#include <stack>
#include <stdlib.h>
using namespace std;

extern int cnt_omorder;
extern int cnt_ominsert;
extern int cnt_omdelete;
extern int cnt_ominsert_mid;
extern int cnt_omsplit;
extern int cnt_tag;
extern int cnt_rebalance_max_tag;
extern int g_debug;
extern unsigned int g_lock_type; // CSR_LOCK by default.


#if 0
#define OMP_PARALLEL_FOR  schedule(dynamic, 10000)
#else 
#define OMP_PARALLEL_FOR
#endif

const int BUSYWAIT_I = 4; // the init value of exponential backoff. 


void power_busy_wait(int &i) { //__attribute__((optimize("O0")))
        int j = i; while (j>0) --j;
        i *= 2;
}


/********************NODE****************************/
ParOM::Node::Node(){
    subtag = 0; // init to 0 
    pre = next = NONE;
    live = true; // node exist
    if (CAS_LOCK == g_lock_type) {
        omlock_cas = lock_cas = 0; //unlock
    } else if (OMP_LOCK == g_lock_type){
        omp_init_lock(&omlock_omp);
    }
}

ParOM::Node::~Node(){
    if (OMP_LOCK == g_lock_type){
        omp_destroy_lock(&omlock_omp);
    }
}

void ParOM::Node::OMLock() {
    if (CAS_LOCK == g_lock_type) {
        int i = BUSYWAIT_I;
        while (true) {
            if (atomic_read(&omlock_cas) == UNLOCK) { // this is faster
                if (cas(&omlock_cas, UNLOCK, LOCK)) {
                    return;
                }
            }
            power_busy_wait(i);
        }
    } else if (OMP_LOCK == g_lock_type){ //omp lock
        omp_set_lock(&omlock_omp);
    }
}


void ParOM::Node::OMUnlock () {
    if (CAS_LOCK == g_lock_type) {
        atomic_write(&omlock_cas, UNLOCK);
    } else if (OMP_LOCK == g_lock_type){
        omp_unset_lock(&omlock_omp);
    }
}

int ParOM::Node::OMTestLock() {
    if (OMP_LOCK == g_lock_type) {
        return omp_test_lock(&omlock_omp);
    }
    return -1;
}
/************************NODE END*********************/

/***********************Group BEGIN*******************/
ParOM::Group::Group() {
    size = 0; pre = next = NONE; tag = 0;
    if (CAS_LOCK == g_lock_type) {
        omlock_cas = 0; //unlock
    } else if (OMP_LOCK == g_lock_type){ //omp lock
        omp_init_lock(&omlock_omp);
    }
}

ParOM::Group::~Group() {
    if (OMP_LOCK == g_lock_type){ //omp lock
        omp_destroy_lock(&omlock_omp);
    }
}

void ParOM::Group::OMLock() {
    if (CAS_LOCK == g_lock_type) {
        while (true) {
            int i = BUSYWAIT_I;
            if (atomic_read(&omlock_cas) == UNLOCK) {
                if (cas(&omlock_cas, UNLOCK, LOCK)) {
                    return;
                }
            }
            power_busy_wait(i);
        }
    } else if (OMP_LOCK == g_lock_type){ //omp lock
        omp_set_lock(&omlock_omp);
    }
}

void ParOM::Group::OMUnlock() {
    if (CAS_LOCK == g_lock_type) {
        atomic_write(&omlock_cas, UNLOCK);
    } else if (OMP_LOCK == g_lock_type){ //omp lock
        omp_unset_lock(&omlock_omp);
    }
}

int ParOM::Group::OMTestLock() {
    if (OMP_LOCK == g_lock_type) {
        return omp_test_lock(&omlock_omp);
    }
    return -1;
}


/***********************Group END*****************/

/***********************OM begin******************/
ParOM::OM::OM(size_t max_size):max_size{max_size}{
    Node node; 
    V = vector<Node>(max_size + 2, node); // Vertices

    Group group; 
    G = vector<Group>(max_size + 2, group);

    H = max_size; T = max_size + 1;
    Hg = max_size; Tg = max_size + 1;

    V[H].next = T; V[T].pre = H;
    V[H].pre = EMPTY; V[T].next = EMPTY;
    
    G[Hg].next = Tg; G[Tg].pre = Hg;
    G[Hg].pre = EMPTY; G[Tg].next = EMPTY;

    V[H].group = Hg; V[T].group = Tg;
    G[Hg].size = 1; G[Tg].size = 1;

    node_size = group_size = 0; // OM is empty

    const subtag_t subtag_gap = MAX_SUBTAG / (MIN_GROUP_SIZE + 1);
    for (int i =0; i < MIN_GROUP_SIZE; i++) {
        const subtag_t subtag = subtag_gap * (i+1);
        assigned_subtags.push_back(subtag);
    }

    if (CAS_LOCK == g_lock_type) {
        omlock_cas = 0; //unlock
    } else if (OMP_LOCK == g_lock_type){ //omp lock
        omp_init_lock(&omlock_omp);
    }
}

ParOM::OM::~OM() {
    if (OMP_LOCK == g_lock_type) {
        omp_destroy_lock(&omlock_omp);
    }
}

void ParOM::OM::Init(vector<node_t> &nodes) {
    for (node_t v: nodes) {  
        // insert to boTgom list after tail.pre p.
        {  
            node_t p = V[T].pre; 
            V[p].next = v; 
            V[v].pre = p; V[v].next = T;
            V[T].pre = v;
        }
        // insert to top list after tail.pre q;
        {
            group_t p= G[Tg].pre;
            G[p].next = v;
            G[v].pre = p; G[v].next = Tg;
            G[Tg].pre = v;
        }
    }

    /*label*/
    tag_t tag_value = 0; 
    for (group_t g= G[Hg].next; g != Tg; g=G[g].next) {
        G[g].tag = tag_value;
        tag_value += INIT_TAG_GAP;
    }
    
    G[Hg].tag = 0;    
    G[Tg].tag = MAX_TAG;  // tail.tag = max tag

    /*sublabel*/
    for (node_t u: nodes) {
        V[u].subtag = INIT_SUBTAG; // all subtag initialized as middle. 
        V[u].group = u; // initially nodes and group have same ids.
        G[u].size = 1; //initially group only has one elements.
    }
    group_size = nodes.size();
}

void ParOM::OM::InitParallel(int num_worker) {
    if (num_worker < 1) {
        printf("***sequential!\n");
        printf("set number of worker: 0\n");
        return;
    }

    printf("set pin CPU\n");
    #pragma omp parallel 
    {
        pthread_t thread;
        thread = pthread_self();
        cpu_set_t CPU;
        CPU_ZERO(&CPU);
        CPU_SET(omp_get_thread_num(), &CPU);
        pthread_setaffinity_np(thread, sizeof(CPU), &CPU);
    }
    if (num_worker > 0 && num_worker < 256) {
        gm_rt_set_num_threads(num_worker); // gm_runtime.h
        printf("***parallel!\n");
        printf("set number of worker: %d\n", gm_rt_get_num_threads());
    } else {
        printf("***parallel! set number of worker error!!!\n");
        printf("set number of worker: -1\n");
    }
}

void ParOM::OM::Lock() {
    if (CAS_LOCK == g_lock_type) {
        while (true) {
            if (atomic_read(&omlock_cas) == UNLOCK) {
                if (cas(&omlock_cas, UNLOCK, LOCK)) {
                    return;
                }
            }
        }
    } else if (OMP_LOCK == g_lock_type){ //omp lock
        omp_set_lock(&omlock_omp);
    }
}


void ParOM::OM::Unlock() {
    if (CAS_LOCK == g_lock_type) {
        atomic_write(&omlock_cas, UNLOCK);
    } else if (OMP_LOCK == g_lock_type){
        omp_unset_lock(&omlock_omp);
    }
}

inline void ParOM::OM::ListDelete(node_t x) {
    node_t pre = V[x].pre; node_t next = V[x].next;
    V[pre].next = next; V[next].pre = pre;
    V[x].pre = V[x].next = NONE;
}

/*insert y after x in list*/
inline void ParOM::OM::ListInsert(node_t x, node_t y) {
    node_t next = V[x].next;
    V[x].next = y; V[y].pre = x;
    V[y].next = next; V[next].pre = y;
}



/*insert all y after x in list*/
inline void ParOM::OM::MultiListInsert(node_t x, vector<node_t> &y) {
    // link all y , scan 0 to size-2
    size_t size = y.size();
    if (1 == size) return ListInsert(x, y[0]);

    for (size_t i = 0; i < size - 1; i++) {  // like the vector
        V[y[i]].next = y[i+1];
        V[y[i+1]].pre = y[i];
    }

    node_t next = V[x].next;
    V[x].next = y[0]; V[y[0]].pre = x;
    V[y[size-1]].next = next; V[next].pre = y[size-1];
}

inline void ParOM::OM::TopListDelete(group_t x) {
    group_t pre = G[x].pre; group_t next = G[x].next;
    G[pre].next = next; G[next].pre = pre;
    //G[x].pre = G[x].next = NONE;
    G[x].pre = G[x].next = -11;  //for debug
}

inline void ParOM::OM::TopListInsert(group_t x, group_t y) {
    group_t next = G[x].next;
    G[x].next = y; G[y].pre = x;
    G[y].next = next; G[next].pre = y;
}

/*insert all y after x in list*/
inline void ParOM::OM::MultiTopListInsert(node_t x, vector<group_t> &y) {
    // link all y , scan 0 to size-2
    size_t size = y.size();
    if (1 == size) return TopListInsert(x, y[0]);

    for (size_t i = 0; i < size - 1; i++) {  // like the vector
        G[y[i]].next = y[i+1];
        G[y[i+1]].pre = y[i];
    }

    node_t next = G[x].next;
    G[x].next = y[0]; G[y[0]].pre = x;
    G[y[size-1]].next = next; G[next].pre = y[size-1];
}

/*check if it has problem
* 1 the core is ok? 
* 2 the tag and subtag may be changed during comparation,
*   which is not allowed. 
* 3 If values change, redo Order*/
inline bool ParOM::OM::Order(node_t x, node_t y, ParOM::Count &cnt) {
    //++cnt_omorder;
AGAIN:
    bool r;
    if (EMPTY == V[x].next || EMPTY == V[y].next) return false;
    tag_t xtag = G[V[x].group].tag; tag_t ytag = G[V[y].group].tag;
    if (xtag != ytag) {
        r =  (xtag < ytag); 
        // in case the relabel happen, the tag is changed by other workers.
        if (xtag != G[V[x].group].tag || ytag != G[V[y].group].tag) {
            ++cnt.order_repeat;
            goto AGAIN;
        }
        else { goto END; }
    } else { //gx == gy
        subtag_t xsub = V[x].subtag; subtag_t ysub = V[y].subtag;  
        r = (xsub < ysub);
        if (xtag != G[V[x].group].tag || ytag != G[V[y].group].tag
                || xsub != V[x].subtag || ysub != V[y].subtag) {
            ++cnt.order_repeat;
            goto AGAIN;
        }
        else { goto END; }
        
    }
END: 
    return r;
}


void ParOM::OM::Insert(node_t x, node_t y, vector<node_t> &relabel_nodes,
            vector<group_t> &groups, ParOM::Count &cnt) {

    ++cnt_ominsert;


    //lock x and x.next with order
    // after x is locked. xnext can't be changed.
    //V[y].OMLock(); // how about other worker insert after y?
    V[x].OMLock(); 
    node_t xnext = V[x].next; 
    V[xnext].OMLock();
   

    //debug
    //assert(Order(x, xnext));

AGAIN:
    group_t g0 = V[x].group;
    const node_t z = V[x].next; // x-> z 
    const subtag_t subtag0 = V[x].subtag;
    subtag_t subw;

    /*there are two cases:
    * 1. z and x in the same group. 
    * 2. z and x in the different groups, y insert append x to x's group
    * 3. z and x in the different groups, y insert before x to z's group*/
    if (likely(V[x].group == V[z].group)) {
        subw = V[z].subtag - subtag0;
    } else {
        subw = MAX_SUBTAG - subtag0; // MAX_SUBTAG + 1, max's next
        if (subw < 2) { // insert y at head of z' group. x -> z
            subw = V[z].subtag; // z at the head of z'group.
            g0 = V[z].group;
        } 
    }
    
    /* it not has enough label space. The relabel is trigered. */ 
    if (unlikely(subw < 2)) {
        // we use global version to STM for PQ. 
        //++g_tag_version[core[x]]; // doing relabel, odd number, "lock" order list
        relabel_nodes.clear(); groups.clear();
        SimpleRelabel2(x, xnext, g0, relabel_nodes, groups, cnt);
        ++cnt.relabel;
        //SimpleRelabel(x);
        //++g_tag_version[core[x]]; // finish, even number, "unlock" order list
        goto AGAIN;
    }

    V[y].subtag = (V[x].subtag + subw / 2);
    V[y].group = g0;
    ++G[g0].size; //if group size = 0, group is removed.

    cnt.subtag++;
    ListInsert(x, y);

    //V[y].OMUnlock(); 
    V[x].OMUnlock(); V[xnext].OMUnlock(); 
}

/*split g0 to miltiple groups
* x: insert y between x and x.next.
* gy: y's group
* relabel_nodes: vector alocate outside*/
void ParOM::OM::SimpleRelabel2(node_t x, node_t xnext, group_t gy, 
            vector<node_t> &relabel_nodes, vector<group_t> &groups, ParOM::Count &cnt) {

    // x and x.next is locked before
    // lock g0
    const group_t g0 = V[x].group; 
    G[g0].OMLock();

    const tag_t tag0 = G[g0].tag;

    assert(tag0 <= (MAX_TAG - 16)); // reach the limit, all nodes require reassign the labels.

    // lock all nodes that needs to split x.next with same group gy
    // relabel_nodes.clear();
    //vector<node_t> relabel_nodes;
    for(node_t u = xnext; gy == V[u].group; u = V[u].next) {
        if (u != xnext) {
            V[u].OMLock();
        }
        relabel_nodes.push_back(u);
    }

    // insert relabel_size groups after g0,
    // find enough tag space
    const size_t relabel_size = relabel_nodes.size(); // total number of relabed node in group.
    size_t relabel_group_size = relabel_size / MIN_GROUP_SIZE;
    if (relabel_size % MIN_GROUP_SIZE > 0) relabel_group_size++;
    
    group_t g0next = G[g0].next;
    G[g0next].OMLock();


    /*rebalance the group tags*/
    group_t g = g0next;
    tag_t w = G[g].tag - tag0;
    size_t j = relabel_size + 1;
    
    group_t first_locked = EMPTY, last_locked = EMPTY;
    if ( unlikely(w <= relabel_group_size) ) { 
        //the insertion walks down the list until w > j^2
        //the gap after relabeling > j. To insert n groups, j has to > n (and double???)
        j = 1;
        while (w <= j * j || j <=  relabel_group_size) { // at least insert j nodes between two. 
            g = G[g].next; 

            G[g].OMLock(); 
           
            if(EMPTY == first_locked) first_locked = g;
            last_locked = g;

            w = G[g].tag - tag0;
            j++;
        }

        // in case traverse to the tail
        if (w > INIT_TAG_GAP) w = INIT_TAG_GAP;

        //relabel the j-1 record s_1 to s_j-1 with the labels 
        tag_t taggap = w/j;
        
        int tag_num = 0;
        for (size_t k = 1, g = G[g0].next; k < j; k++, g=G[g].next) {
            // the wide gap w is splited into j parts by w/j, 
            G[g].tag = taggap * k  + tag0; 
            cnt.tag++;
            tag_num++;
        }
        if(tag_num > cnt_rebalance_max_tag) cnt_rebalance_max_tag = tag_num;

        //reset w after relabel
        g = G[g0].next;
        w = G[g].tag - tag0; 
    }

    // for(group_t g: locked_groups) { 
    //     G[g].OMUnlock();
    // }
    for (group_t  g = first_locked; g != last_locked; g=G[g].next) {
        G[g].OMUnlock();
    }
    G[last_locked].OMUnlock();

    /* in case x is the tail, scope the tag. To a reasonable range.
    * INIT_TAG_GAP / 32 is test value*/ 
    //if (w > INIT_TAG_GAP) w = min(INIT_TAG_GAP, j * j * relabel_size * 2);
    if (w > INIT_TAG_GAP) w = INIT_TAG_GAP;


    /*insert relabel_size group after g0, consider later. Ignore the group size is 0*/
    //vector<group_t> groups;

    /*allocate groups should be locked*/
    Lock();
    size_t start_group = atomic_read(&group_size);
    size_t end_group = start_group + relabel_group_size;
    atomic_write(&group_size, end_group);
    Unlock();
    
    for(size_t group_id = start_group; group_id < end_group; group_id++) {
        assert(group_id < G.size());
        groups.push_back(group_id);
    }

    MultiTopListInsert(g0, groups); //lock g0 and g0next already.
    
    /*reversely nodes to group to keep order invariant*/
    tag_t taggap = w / (relabel_group_size + 1);
    for (int k = relabel_size-1; k >=0; k--) {
        const group_t g = groups[k / MIN_GROUP_SIZE];
        const node_t u = relabel_nodes[k];
        V[u].group = g;
        
        /*assign label for each group*/
        //V[u].subtag = MAX_SUBTAG / 2;
        if (0 == k % MIN_GROUP_SIZE) { //assign label
            const int begin = k; // the fisrt nodes in the group g
            const int end = std::min(relabel_size, k + MIN_GROUP_SIZE); // the last+1 nodes in group g.
            int bottom = begin; // stack bottom. 


            G[g].tag = tag0 + taggap * (k/MIN_GROUP_SIZE + 1);
            ++cnt.tag;
            
            G[g].size = end - begin;
            for (int i = begin; i < end; i++) {
                const node_t v = relabel_nodes[i];
                if (V[v].subtag >= assigned_subtags[i%MIN_GROUP_SIZE]) 
                {
                    for (int j = i; j >=bottom; j--) {
                        const node_t v2 = relabel_nodes[j];
                        V[v2].subtag = assigned_subtags[j%MIN_GROUP_SIZE];
                        ++cnt.subtag;
                    }
                    bottom = i + 1;
                }
            }

        }

    }

    for(node_t u: relabel_nodes) { 
        if (u != xnext) {
            V[u].OMUnlock();
        } 
    }

    G[g0].OMUnlock(); G[g0next].OMUnlock();

}



/*return 0, success delete x
* return 1, x is already deleted by other workers.*/
bool ParOM::OM::Delete(node_t x) {
    //++cnt_omdelete;

DOLOCK:
    //lock x.pre, x and x.next with order
    node_t xpre = V[x].pre;
    if (NONE == xpre) return false; // x is deleted.
    V[xpre].OMLock(); 
    if (xpre != V[x].pre) { // in case xpre is changed.
        V[xpre].OMUnlock();
        //printf("delete do lock \n");
        goto DOLOCK;
    } 
    V[x].OMLock(); 
    node_t xnext = V[x].next;
    V[xnext].OMLock();

    const group_t g0 = V[x].group;
    
    // update group size, delete when group is empty. 
    // not affect the correctness but left empty groups for recycling. 
    /*
    --G[g0].size; 
    atomic_sub(&G[g0].size, 1);
    assert(G[g0].size >= 0);
    if(0 == G[g0].size) {
        TopListDelete(g0); 
    }
    */

    ListDelete(x); 
    V[x].group = EMPTY; // removed no group. 
    
    //unlock
    V[xpre].OMUnlock(); V[x].OMUnlock(); V[xnext].OMUnlock();
    return true;
}

/*by setting flag to delete, next time to recycle. 
* it avoid to use lock*/
bool ParOM::OM::DeleteByFlag(node_t x) {
    return cas(&V[x].live, true, false);
    // V[x].OMLock();
    // V[x].live = false;
    // V[x].OMUnlock();
}
/**********************OM end *******************/


/***********************Test************************/
vector<int> ParOM::OM::GetRepeatRandomNum(char * const path) {
    size_t found;
    string str(path);
    found=str.find_last_of("/\\");
    str=str.substr(0,found); //folder str.substr(found+1) file
    str+="/random-1m.txt";
    std::vector<int> rand = gadget::RepeateRandomRead(str.c_str());
    return rand;
}

vector<int> ParOM::OM::GenerateTestCase(TestCase t, vector<int> &randnum, size_t insertNum) {
    vector<int> pos;
    size_t rand_size = randnum.size();
    srand (time(NULL));
    if (NO_RELABEL == t) 
    {
        while(pos.size() < insertNum) {
            for (size_t i: randnum) {
                //i = (i*i);
                pos.push_back(i % insertNum); // 1m positions.
                if (pos.size() >= insertNum) break;
            }
        }
    }
    else if( FEW_RELABEL == t ) // repeate random insert 10 m to 1m items  
    {
         /*1m rand position*/
        vector<int> pos2;
        while (pos2.size() < OM_POS_FEW) {
            for (size_t i: randnum) {
                i = (i * i);
                pos2.push_back(i % OM_SIZE); // 1m positions.
                if (pos2.size() >= OM_POS_FEW) break;
            }
        }

        /*10m inset number*/
        while(pos.size() < insertNum) {
            for (size_t i: randnum) {
                i = (i * i);
                pos.push_back(i % OM_POS_FEW); // 1m positions.
                if (pos.size() >= insertNum) break;
            }
        }
    } 
    else if (MANY_RELABEL == t) 
    {
        /*1000 rand position*/
        vector<int> pos2;
        int k = 13;
        for (size_t i: randnum) {
            i = (i * i);
            pos2.push_back(i % OM_SIZE); // 1m positions.
            if (pos2.size() >= OM_POS_MANY) break;
        }        
        
        while(pos.size() < insertNum) {
            for (size_t i: randnum) {
                i = (i * i);
                pos.push_back(pos2[i % OM_POS_MANY]); // 1m positions.
                if (pos.size() >= insertNum) break;
            }
        }
    } 
    else if (MAX_RELABEL == t) 
    {
        // midle of the list
        int position = OM_SIZE/2;
        for (int i = 0; i < insertNum; i++) { //OM_INSERT_SIZE
            pos.push_back(position); // fixed position 
        }
    }

    return pos;
}

/* test cases with random number. 
* for scaliable test, the radio is same for different scales. 
**/
vector<int> ParOM::OM::GenerateTestCase2(TestCase t, size_t insertNum) {
    vector<int> pos;
    srand (time(NULL));

    if (NO_RELABEL == t) 
    {   
        /*10m rand position for 10m data, 1 per position*/
        for (size_t i = 0; i < insertNum; i++) {
            int j = rand() % insertNum; j = abs(j);
            pos.push_back(j);
        }
    }
    else if( FEW_RELABEL == t ) // repeate random insert 10 m to 1m items  
    {
         /*1m rand position for 10m data, 10 per position*/ 
        vector<int> pos2;
        int position_size = insertNum/10;
        for (size_t i = 0; i < position_size; i++) {
            int j = rand() % position_size; j = abs(j);
            pos2.push_back(j); // 1m positions.
        }
        

        /*10m inset number*/
        for (size_t i = 0; i < insertNum; i++) { 
            int j = rand() % position_size; j = abs(j);
            pos.push_back(pos2[j]); // 1m positions.
        }
    } 
    else if (MANY_RELABEL == t) 
    {
        /*1000 rand position for 10m data, 10,000 per position*/
        vector<int> pos2;
        int position_size = insertNum/10000;
        for (size_t i = 0; i < position_size; i++) {
            int j = rand() % position_size; j = abs(j);
            pos2.push_back(j); // 1m positions.
        }
        

        /*10m inset number*/
        for (size_t i = 0; i < insertNum; i++) { 
            int j = rand() % position_size; j = abs(j);
            pos.push_back(pos2[j]); // 1m positions.
        }
    } 
    else if (MAX_RELABEL == t) 
    {
        // midle of the list
        int position = insertNum/2;
        for (int i = 0; i < insertNum; i++) { //OM_INSERT_SIZE
            pos.push_back(position); // fixed position 
        }
    }

    return pos;
}
vector<node_t> ParOM::OM::AlocateNodes(int num) {
    vector<node_t> nodes;
    for(size_t id = node_size; id < node_size + num; id++) {
        assert(id < max_size);
        nodes.push_back(id);
    }
    node_size += num;
    return nodes;
}


void ParOM::OM::TestInsert(vector<int> &pos, vector<int> &nodes, bool ispar) {
    const size_t size = pos.size();
    if (ispar) {    
    
    #pragma omp parallel
    {
        int workerid = omp_get_thread_num();
        vector<node_t> relabel_nodes;
        vector<group_t> groups;
        ParOM::Count cnt;
        #pragma omp critical 
        {
            relabel_nodes.reserve(1024*2);
            groups.reserve(1024*2);
        }
        //printf("insert workerid %d\n", workerid);
        #pragma omp for OMP_PARALLEL_FOR
        for (int i = 0; i < size; i++) 
        {   
            Insert(pos[i], nodes[i], relabel_nodes, groups, cnt);
        }

        #pragma omp critical 
        {
            g_cnt = g_cnt + cnt;
        }
    }

    } else {
        vector<node_t> relabel_nodes;
        vector<group_t> groups;
        relabel_nodes.reserve(1024*2);
        groups.reserve(1024*2);
        for (int i = 0; i < size; i++) 
        {
            Insert(pos[i], nodes[i], relabel_nodes, groups, g_cnt);
        }
    }
}

void ParOM::OM::TestDelete(vector<int> &pos, bool ispar) {
    const size_t size = pos.size();
    if (ispar) {
    
    #pragma omp parallel
    {
        int workerid = omp_get_thread_num();
        #pragma omp for OMP_PARALLEL_FOR
        for (int i = 0; i < size; i++) 
        {
            Delete(pos[i]);
        }
    }    

    } else {
        for (int i = 0; i < size; i++) 
        {
            Delete(pos[i]);
        }
    }
}

void ParOM::OM::TestOrder(vector<int> &pos, bool ispar) {
    const size_t size = pos.size();
    if (ispar) {
    #pragma omp parallel
    {
        int workerid = omp_get_thread_num();
        ParOM::Count cnt;

        #pragma omp for OMP_PARALLEL_FOR
        for (int i = 0; i < size; i++) 
        {   
            Order(pos[i], pos[(i+1)%size], cnt);
        }

        #pragma omp critical 
        {
            g_cnt = g_cnt + cnt;
        }
    }   

    } else {
        for (int i = 0; i < size; i++) 
        {   
            Order(pos[i], pos[(i+1)%size], g_cnt);
        }
    }
}

void ParOM::OM::TestMixed(vector<int> &pos, vector<int> &nodes, bool ispar) {
    const size_t size = pos.size();

    if (ispar) {
    #pragma omp parallel
    {
        int workerid = omp_get_thread_num();
        vector<node_t> relabel_nodes;
        vector<group_t> groups;

        // vector<node_t> stack_vec; // use vector as a stack.
        // int stack_pos = 0;

        ParOM::Count cnt;

        #pragma omp critical 
        {
            groups.reserve(1000000);
            relabel_nodes.reserve(1000000);
            //stack_vec.reserve(pos.size());
        }
        #pragma omp for OMP_PARALLEL_FOR
        for (int i = 0; i < size; i++) 
        {   
            /*the mix operation is to test Order. how to test order repeat #??? */
            for (int k = 1; k <= 10; k++) {
                Order(pos[i], pos[(i+k) % size], cnt);
            }
            Insert(pos[i], nodes[i], relabel_nodes, groups, cnt);
            

            // int opt = rand[i%rand_size] % (1024*3);
            // if (opt < 1024) {
            //     Order(pos[i], pos[(i+1) % size], cnt);
            // } else if (opt >= 1024 * 2){
            //     //stack_vec.push_back(nodes[i]); ++stack_pos;
            //     Insert(pos[i], nodes[i], relabel_nodes, groups, cnt);
            // } else {
            //     // while(stack_pos > 0) {
            //     //     node_t v = stack_vec[stack_pos-1]; --stack_pos;  
            //     //     Delete(v);
            //     // }
            //     // stack_vec.clear();
            // }

            // opt = rand[abs(i*3)%rand_size] % 1024;
            // if (opt < 512) {
            //     Order(pos[size - 1 - i], pos[(size - i) % size], cnt);
            // } else {
            //     Insert(pos[i], nodes[i], relabel_nodes, groups, cnt);
            //     Delete(nodes[i]);
            // } 
        }
        
        #pragma omp critical 
        {
            g_cnt = g_cnt + cnt;
        }

    }   

    } else {
        vector<node_t> relabel_nodes;
        vector<group_t> groups;
        vector<node_t> stack_vec;
        stack<node_t, vector<node_t>> stack(stack_vec);

        relabel_nodes.reserve(1000000);
        groups.reserve(1000000);
        nodes.reserve(pos.size());


        for (int i = 0; i < size; i++) 
        { 
            for (int k = 1; k <= 10; k++) {
                Order(pos[i], pos[(i+k) % size], g_cnt);
            }
            Insert(pos[i], nodes[i], relabel_nodes, groups, g_cnt);
            
            // int opt = rand[i%rand_size] % (1024*3);
            // if (opt < 1024) {
            //     Order(pos[i], pos[(i+1) % size], g_cnt);
            // } else if (opt >= 1024 * 2){
            //     stack.push(nodes[i]);
            //     Insert(pos[i], nodes[i], relabel_nodes, groups, g_cnt);
            // } else {
            //     while(!stack.empty()) {
            //         node_t v = stack.top(); stack.pop(); 
            //         Delete(v);
            //     }
            // }

            // opt = rand[abs(i*3)%rand_size] % 1024;
            // if (opt < 512) {
            //     Order(pos[size - 1 - i], pos[(size - i) % size], g_cnt);
            // } else {
            //     Insert(pos[i], nodes[i], relabel_nodes, groups, g_cnt);
            //     //Delete(nodes[i]);
            // } 

        }
    }

   
}


void ParOM::OM::PrintCount(char* str) {
    
    printf("\n");
    printf(str); printf("tag #: %u\n", g_cnt.tag);
    printf(str); printf("subtag #: %u\n", g_cnt.subtag);
    printf(str); printf("order repeat #: %u\n", g_cnt.order_repeat);
    printf(str); printf("relabel #: %u\n", g_cnt.relabel);
    printf("\n");
}

/***********************Test END************************/

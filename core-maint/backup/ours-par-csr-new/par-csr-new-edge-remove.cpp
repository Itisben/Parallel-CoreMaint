#include "def.h"
#include "par-csr-new.h"
#include <algorithm>


//for counting 

// put all of this into the class, so that avoid to atomic add or deleted.
extern int cnt_Vs; // count of size V*
extern int cnt_Vp; // count of size V+
extern int cnt_S;  // count the number of all visited edges
extern int cnt_tag; // count all update tag. 
extern int cnt_rebalance_max_tag; // count the update tag during relabel process. 
extern int cnt_edge; // for debug
extern int cnt_edge2; // for debug

extern int cnt_omorder;
extern int cnt_ominsert;
extern int cnt_omdelete;
extern int cnt_ominsert_mid;
extern int cnt_omsplit;
extern int cnt_Q_size;

extern int g_debug;
extern unsigned int g_lock_type; // CSR_LOCK by default.

extern int optTest;


volatile int s_lock = 0; // test with the global lock.

// omp_lock_t s_lock_omp;
// omp_init_lock(&s_lock_omp);

//busy waiting 
const int BUSYWAIT_I = 1; // the init value of exponential backoff. 
static void power_busy_wait(int &i) { //__attribute__((optimize("O0")))
        int j = i; while (j>0) --j;
        i *= 2;
}

static void lock(lock_t *caslock, omp_lock_t *omplock) {
    if (CAS_LOCK == g_lock_type) {
        int i = BUSYWAIT_I;
        while (true) {
            if (atomic_read(caslock) == UNLOCK) {
                if (cas(caslock, UNLOCK, LOCK)) {
                    return;
                }
            }
            power_busy_wait(i);
        }
    } else if (OMP_LOCK == g_lock_type) {
        omp_set_lock(omplock);
    }
}

extern std::vector<int> Vcount; // count the Vblack size frequency.

/********************NODE****************************/
ParCM::Node::Node(){
    degout = degin = subtag = 0; // init to 0 
    pre = next = NONE;
    color = WHITE;
    inQ = false;
    inR = false;
    mcd = EMPTY;
    islock = NONE;

    s = 0;
    s_lock = 0;

    if (CAS_LOCK == g_lock_type) {
        omlock_cas = lock_cas = 0; //unlock
    } else if (OMP_LOCK == g_lock_type){
        omp_init_lock(&omlock_omp);
        omp_init_lock(&lock_omp);
    }
}

ParCM::Node::~Node(){
    if (OMP_LOCK == g_lock_type){
        omp_destroy_lock(&omlock_omp);
        omp_destroy_lock(&lock_omp);
    }
}

void ParCM::Node::OMLock() {
    if (likely(CAS_LOCK == g_lock_type)) {
        int i = BUSYWAIT_I;
        while (true) {
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

void ParCM::Node::OMUnlock () {
    if (likely(CAS_LOCK == g_lock_type)) {
        atomic_write(&omlock_cas, UNLOCK);
    } else if (OMP_LOCK == g_lock_type){
        omp_unset_lock(&omlock_omp);
    }
}


bool ParCM::Node::OMTestLock() {
    if (likely(CAS_LOCK == g_lock_type)) {
        if (atomic_read(&omlock_cas) == UNLOCK) {
            return cas(&omlock_cas, UNLOCK, LOCK);
        } else return false;
        return false;
    } else if (OMP_LOCK == g_lock_type) { //omp lock
        return omp_test_lock(&omlock_omp);
    }
    return true;
}


void ParCM::Node::Lock() {
    if (likely((CAS_LOCK == g_lock_type))) {
        int i = BUSYWAIT_I;
        while (true) {
            if (atomic_read(&lock_cas) == UNLOCK) {
                if (cas(&lock_cas, UNLOCK, LOCK)) {
                    return;
                }
            }
            power_busy_wait(i);
        }
    } else if (OMP_LOCK == g_lock_type){ //omp lock
        omp_set_lock(&lock_omp);
    }
}

// void ParCM::Node::NestLock() {

// }

// void ParCM::Node::NestUnlock() {

// }

// bool ParCM::Node::Lockp(int p) {
//     if ( likely(CAS_LOCK == g_lock_type) ) {
//         int i = BUSYWAIT_I;
//         while (atomic_read(&lock_cas) != p) {
//             if (atomic_read(&lock_cas) == UNLOCK && cas(&lock_cas, UNLOCK, p)) {
//                     return true; // lock successfully.
//             }
//             power_busy_wait(i);
//         }
//         return false; // already locked by worker p.
//     } else if (OMP_LOCK == g_lock_type){ //omp lock
//         omp_set_lock(&lock_omp);
//     }
// }

/*this is only work for cas lock*/
// bool ParCM::Node::TestCoreLock(core_t K, core_t & current_k) {
//     if (likely(CAS_LOCK == g_lock_type)) {
//         while (K == current_k) {
//             if (lock_cas == UNLOCK) {
//                 if (cas(&lock_cas, UNLOCK, LOCK)) {
//                     return true;
//                 }
//             }
//         }
//         return false;
//     } else if (OMP_LOCK == g_lock_type){ //omp lock
//         omp_set_lock(&omlock_omp);
//     }
// }

bool ParCM::Node::LockWithCore(const core_t k, core_t *core) {
    if (likely(CAS_LOCK == g_lock_type)) {
        int i = BUSYWAIT_I;
        while (k == atomic_read(core)) {
            if (atomic_read(&lock_cas) == UNLOCK && cas(&lock_cas, UNLOCK, LOCK)) {
                if (k == atomic_read(core)) {
                    return true;
                } else {
                    atomic_write(&lock_cas, UNLOCK);
                    return false;
                }
                
                
            }
            power_busy_wait(i);
        }
        return false;
    } else if (OMP_LOCK == g_lock_type){ //omp lock
        while (k == atomic_read(core)) {
            if (omp_test_lock(&lock_omp)) { 
                if (k == atomic_read(core)) {
                    return;
                } else {
                    omp_unset_lock(&lock_omp);
                    return false;
                }
            } 
        }
        return false;
    }
    return true;
}

void ParCM::Node::Unlock () {
    if (likely(CAS_LOCK == g_lock_type)) {
        atomic_write(&lock_cas, UNLOCK);
    } else if (OMP_LOCK == g_lock_type){
        omp_unset_lock(&lock_omp);
    }
}

/*we can add power busy wait for optimization*/
bool ParCM::Node::TestLock() {
    if (likely(CAS_LOCK == g_lock_type)) {
        if (atomic_read(&lock_cas) == UNLOCK) {
            return cas(&lock_cas, UNLOCK, LOCK);
        } else return false;
    } else if (OMP_LOCK == g_lock_type) { //omp lock
        return omp_test_lock(&lock_omp);
    }
    return true;
}

bool ParCM::Node::IsLocked() {
    if (likely(CAS_LOCK == g_lock_type)) {
        return atomic_read(&lock_cas) > 0;
    } else if (OMP_LOCK == g_lock_type) { //omp lock
        return atomic_read(&lock_omp._x[0]) > 0;
    }
    return true;
}
/************************NODE END*********************/

/***********************Group BEGIN*******************/
ParCM::Group::Group() {
    size = 0; pre = next = NONE; tag = 0;
    if (CAS_LOCK == g_lock_type) {
        omlock_cas = 0; //unlock
    } else if (OMP_LOCK == g_lock_type) {
        omp_init_lock(&omlock_omp);
    }
}

ParCM::Group::~Group() {
    if (OMP_LOCK == g_lock_type){ //omp lock
        omp_destroy_lock(&omlock_omp);
    }
}


void ParCM::Group::OMLock() {
    if (likely(CAS_LOCK == g_lock_type)) {
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

void ParCM::Group::OMUnlock() {
    if (likely(CAS_LOCK == g_lock_type)) {
        atomic_write(&omlock_cas, UNLOCK);
    } else if (OMP_LOCK == g_lock_type) {
        omp_unset_lock(&omlock_omp);
    }
}


/***********************Group END*****************/



/********ParCM************************************/

void ParCM::CoreMaint::OMLock() {
    if (likely(CAS_LOCK == g_lock_type)) {
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


void ParCM::CoreMaint::OMUnlock() {
    if (likely(CAS_LOCK == g_lock_type)) {
        atomic_write(&omlock_cas, UNLOCK);
    } else if (OMP_LOCK == g_lock_type){
        omp_unset_lock(&omlock_omp);
    }
}


ParCM::CoreMaint::CoreMaint(size_t n, PARGRAPH &graph, vector<core_t> &core):
        n{n}, graph{graph}, core{core} {
    
    const size_t reserve_size = n; //2^15
    /* allocated outside
    Vblack.reserve(reserve_size);
    Vcolor.reserve(reserve_size);
    
    R = QUEUE(reserve_size); //2^14
    PQ = PRIORITY_Q(reserve_size);

    groups.reserve(1024);
    locked_groups.reserve(1024);
    */

  
    /* init the Order List by k-order odv */
    max_core = 0;
    for(core_t k: core) {  // find max core
        if (k > max_core) max_core = k;
    }

    max_core += MAX_INCREASE_CORE; // the maximum core number after increasing. 

    Node node; 
    V = vector< Node>(n + 2*max_core + 2, node); // Vertices
    H = vector<node_t>(max_core + 1); // head of order list k
    T = vector<node_t>(max_core + 1); // tail of order list k

    // at most n group can be new inserted
    Group group; 
    G = vector<Group>(2 * n + 2*max_core + 2, group);
    group_size = n; // initially n vertices have n groups.
    Hg = vector<node_t>(max_core + 1); // head of group list k  
    Tg = vector<node_t>(max_core + 1); // tail of group list k

    g_tag_version = vector<ver_t>(max_core+1, 0); 

    /*initialized group subtages*/
    const subtag_t subtag_gap = MAX_SUBTAG / (MIN_GROUP_SIZE + 1);
    for (int i =0; i < MIN_GROUP_SIZE; i++) {
        const subtag_t subtag = subtag_gap * (i+1);
        assigned_subtags.push_back(subtag);
    }

}

ParCM::CoreMaint::CoreMaint::~CoreMaint(){
    
}

void ParCM::CoreMaint::Init(vector<int> &odv) {
    /* head include a empty node to avoid empty check.
    * tail the same as head. 
    * head -> v -> tail*/
    for (size_t k = 0; k < H.size(); k++) { // init order list
        H[k] = n + 2*k; // bottom list head. 
        T[k] = n + 2*k + 1; // bottom list tail. 
        
        //head -> tail: two nodes initially.
        V[H[k]].next = T[k];    V[H[k]].pre = NONE; 
        V[T[k]].pre = H[k];     V[T[k]].next = NONE;

        Hg[k] = 2*n + 2*k;      //top group list head
        Tg[k] = 2*n + 2*k + 1;  //top group list tail
        G[Hg[k]].next = Tg[k];  G[Hg[k]].pre = NONE;
        G[Tg[k]].pre = Hg[k];   G[Tg[k]].next = NONE;
    }

    /* link all group nodes by the k-order (odv) with head and tail
    * initially group has same number of nodes*/
    for (node_t v: odv) { 
        core_t k = core[v]; 
        // insert to bottom list after tail.pre p.
        {  
            node_t p = V[T[k]].pre; 
            V[p].next = v; 
            V[v].pre = p; V[v].next = T[k];
            V[T[k]].pre = v;
        }
        // insert to top list after tail.pre q;
        {
            group_t p= G[Tg[k]].pre;
            G[p].next = v;
            G[v].pre = p; G[v].next = Tg[k];
            G[Tg[k]].pre = v;
        }
    }


    /*assign the tag for all k list, head.tag = 0, tail.tag = maxtag
    * all tags are in the middle between head and tail
    * all subtags are same as 2^16*/
    int core_size = H.size(); // H.size is the number of cores. 
    for (size_t k = 0; k < core_size; k++) { 
        group_t head = Hg[k]; group_t tail = Tg[k];

        /*tag increasing by gap starting from middle,
        * to reduce frequency of the rebalance*/
        size_t k_core_size = 0;
        for (group_t g= G[head].next; g != tail; g=G[g].next) {
            k_core_size++;
        }
        tag_t tag_value = (MAX_TAG - k_core_size * INIT_TAG_GAP)/ 2; 
        for (group_t g= G[head].next; g != tail; g=G[g].next) {
            G[g].tag = tag_value;
            tag_value += INIT_TAG_GAP;
        }
        G[head].tag = 0;    
        G[tail].tag = MAX_TAG;  // tail.tag = max tag
        V[H[k]].group = head;   // subhead.group = tophead.group.
        V[T[k]].group = tail;   //subtail.group = toptail.group
    }
    
    for (size_t u = 0; u < n; u++) {
        V[u].subtag = INIT_SUBTAG; // all subtag initialized as middle. 
        V[u].group = u; // initially nodes and group have same ids.
        G[u].size = 1; //initially group only has one elements.
    }
    
    
    /*init the mcd at the beginning. The mcd set to EMPTY when core number changed. 
    * Remove for test*/
#if 1
    if (1 != optTest) { // test remove without init mcd. by set "-t 1"
        for (size_t u = 0; u < n; u++) {
            deg_t mcd = 0;
            //for (node_t v: graph[u]) {
            for (node_t idx = graph.begin[u]; idx < graph.end[u]; idx++){
                node_t v = graph.node_idx[idx];
                if (core[u] <= core[v]) mcd++;
            }
            V[u].mcd = mcd;
        }
    }
#endif 

    /*init for parallel*/
}

void ParCM::CoreMaint::InitParallel(int num_worker) {
    if (num_worker < 1) return;
    
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
    if (num_worker > 0 && num_worker < 65) {
        gm_rt_set_num_threads(num_worker); // gm_runtime.h
        printf("set number of worker %d, %d!\n", num_worker, 
                gm_rt_get_num_threads());
    } else {
        printf("set number of worker wrong %d!\n", num_worker);
    }
}

void ParCM::CoreMaint::ComputeCore(std::vector<int>& core, std::vector<int>& order_v, bool with_init) {
  // compute the cores
  auto& deg = core;
  int max_deg = 0;
  int n_ = n;
  for (int i = 0; i < n_; ++i) {
    //deg[i] = graph[i].size();
    deg[i] = graph.size(i);

    if (deg[i] > max_deg) {
      max_deg = deg[i];
    }
  }
  std::vector<int> bin(max_deg + 1, 0);
  for (int i = 0; i < n_; ++i) {
    ++bin[deg[i]];
  }
  int start = 0;
  for (int i = 0; i <= max_deg; ++i) {
    int temp = bin[i];
    bin[i] = start;
    start += temp;
  }
  std::vector<int> vert(n_);
  std::vector<int> pos(n_);
  for (int i = 0; i < n_; ++i) {
    pos[i] = bin[deg[i]];
    vert[pos[i]] = i;
    ++bin[deg[i]];
  }
  for (int i = max_deg; i > 0; --i) {
    bin[i] = bin[i-1];
  }
  bin[0] = 0;
  int k = 0;
  auto vis = std::vector<bool>(n_, false);
  for (int i = 0; i < n_; ++i) {
    const int v = vert[i];
    if (deg[v] > k) k = deg[v];
    ASSERT(bin[deg[v]] == i);
    ++bin[deg[v]];
    core[v] = k; // get key
    
    //init k-order
    if(with_init) {order_v[i] = v;} // get k-order O

    vis[v] = true;
    int rem = 0;

    //for (const int u : graph[v]) {
    for (node_t idx = graph.begin[v]; idx < graph.end[v]; idx++){
        node_t u = graph.node_idx[idx];
    //GRAPH_EDGE(v, u)
        if (vis[u]) continue;
        ++rem;
        const int pw = bin[deg[u]];
        const int pu = pos[u];
        if (pw != pu) {
            const int w = vert[pw];
            vert[pu] = w;
            pos[w] = pu;
            vert[pw] = u;
            pos[u] = pw;
        }
        ++bin[deg[u]];
        --deg[u];
        if (pos[u] == i + 1) {
            bin[deg[u]] = pos[u];
        }
    GRAPH_EDGE_END
    /*init degout, here is more efficient that check the order of u and v*/
    if(with_init) {V[v].degout = rem;}
  }
}


inline void ParCM::CoreMaint::ListDelete(node_t x) {
    node_t pre = V[x].pre; node_t next = V[x].next;
    V[pre].next = next; V[next].pre = pre;
    V[x].pre = V[x].next = NONE;
}

/*insert y after x in list*/
inline void ParCM::CoreMaint::ListInsert(node_t x, node_t y) {
    node_t next = V[x].next;
    V[x].next = y; V[y].pre = x;
    V[y].next = next; V[next].pre = y;
}


/*insert all y after x in list*/
inline void ParCM::CoreMaint::MultiListInsert(node_t x, vector<node_t> &y) {
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

inline void ParCM::CoreMaint::TopListDelete(group_t x) {
    group_t pre = G[x].pre; group_t next = G[x].next;
    G[pre].next = next; G[next].pre = pre;
    //G[x].pre = G[x].next = NONE;
    G[x].pre = G[x].next = -11;  //for debug
}

inline void ParCM::CoreMaint::TopListInsert(group_t x, group_t y) {
    group_t next = G[x].next;
    G[x].next = y; G[y].pre = x;
    G[y].next = next; G[next].pre = y;
}

/*insert all y after x in list*/
inline void ParCM::CoreMaint::MultiTopListInsert(node_t x, vector<group_t> &y) {
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
inline bool ParCM::CoreMaint::Order(node_t x, node_t y) {
    //++cnt_omorder;
AGAIN:
    bool r;
    core_t cx = atomic_read(&core[x]); core_t cy = atomic_read(&core[y]);
    if (cx == cy) {
        tag_t xtag = G[V[x].group].tag; tag_t ytag = G[V[y].group].tag;
        if (xtag != ytag) {
            r =  (xtag < ytag); 
            // in case the relabel happen, the tag is changed by other workers.
            if (xtag != G[V[x].group].tag || ytag != G[V[y].group].tag) {
               
                goto AGAIN;
            }
            else { goto END; }
        } else { //gx == gy
            subtag_t xsub = V[x].subtag; subtag_t ysub = V[y].subtag;  
            r = (xsub < ysub);
            if (xtag != G[V[x].group].tag || ytag != G[V[y].group].tag
                    || xsub != V[x].subtag || ysub != V[y].subtag) {
             
                goto AGAIN;
            }
            else { goto END; }
            
        }
    } else {
        r = (cx < cy);
    }
    // if core number changed redo
    if (cx != atomic_read(&core[x]) || cy !=  atomic_read(&core[y])) {
        goto AGAIN;
    }
END: 
    return r;
}

/*check if it has problem due to node has two labels*/
bool ParCM::CoreMaint::SameCoreOrder(node_t x, node_t y) {
AGAIN:
   
    core_t cx = atomic_read(&core[x]); core_t cy = atomic_read(&core[y]);
    if (cx != cy) return false;
    
    bool r;
    if (EMPTY == V[x].next || EMPTY == V[y].next) return false;
    tag_t xtag = G[V[x].group].tag; tag_t ytag = G[V[y].group].tag;
    if (xtag != ytag) {
        r =  (xtag < ytag); 
        // in case the relabel happen, the tag is changed by other workers.
        if (xtag != G[V[x].group].tag || ytag != G[V[y].group].tag) {
            
            goto AGAIN;
        }
        else { goto END; }
    } else { //gx == gy
        subtag_t xsub = V[x].subtag; subtag_t ysub = V[y].subtag;  
        r = (xsub < ysub);
        if (xtag != G[V[x].group].tag || ytag != G[V[y].group].tag
                || xsub != V[x].subtag || ysub != V[y].subtag) {
            
            goto AGAIN;
        }
        else { goto END; }
        
    }
END: 
    return r;
}

/*insert y after x. x->z 
* each list has head and tail as empty nodes. 
* j = 1
* while wj <= j*j :
*   update wj
*   j++
*
* ??? do I need to lock y (inserted node) 
* lock x.group when spliting
* NOTE: 
* 1 lock the multiple nodes folowing the order to avoid the dead-lock. 
* 2 split group, I use a simple strategy. 
* 3 
*/
void ParCM::CoreMaint::OrderInsert(node_t x, node_t y) {
    ++cnt_ominsert;

    V[x].OMLock(); 
    node_t xnext = V[x].next; 
    V[xnext].OMLock();

AGAIN:
    group_t g0 = V[x].group;
    const node_t z = V[x].next; // x-> z 
    const subtag_t subtag0 = V[x].subtag;
    subtag_t subw;

    /*there are two cases:
    * 1. z and x in the same group. 
    * 2. z and x in the different group, y insert append to x within group
    * 3. z and x in different group, */
    if (unlikely(V[x].group == V[z].group)) {
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
        SimpleRelabel2(x, xnext, g0);
        //SimpleRelabel(x);
        //++g_tag_version[core[x]]; // finish, even number, "unlock" order list
        goto AGAIN;
    }

    V[y].subtag = (V[x].subtag + subw / 2);
    V[y].group = g0;
    ++G[g0].size; //if group size = 0, group is removed.

    cnt_tag++;
    ListInsert(x, y);

    V[x].OMUnlock(); V[xnext].OMUnlock(); 
}
/*split g0 to miltiple groups
* x: insert y between x and x.next.
* gy: y's group*/
void ParCM::CoreMaint::SimpleRelabel2(node_t x, node_t xnext, group_t gy) {
    printf(" ==== Relabel by %d\n", x);

  // x and x.next is locked before
    // lock g0
    const group_t g0 = V[x].group; 
    G[g0].OMLock();

    const tag_t tag0 = G[g0].tag;

    assert(tag0 <= (MAX_TAG - 16)); // reach the limit, all nodes require reassign the labels.

    // lock all nodes that needs to split x.next with same group gy
    // relabel_nodes.clear();

    vector<node_t> relabel_nodes;
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
            //cnt.tag++;
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
    vector<group_t> groups;

    /*allocate groups should be locked*/
    OMLock();
    size_t start_group = atomic_read(&group_size);
    size_t end_group = start_group + relabel_group_size;
    atomic_write(&group_size, end_group);
    OMUnlock();
    
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
            //++cnt.tag;
            
            G[g].size = end - begin;
            for (int i = begin; i < end; i++) {
                const node_t v = relabel_nodes[i];
                if (V[v].subtag >= assigned_subtags[i%MIN_GROUP_SIZE]) 
                {
                    for (int j = i; j >=bottom; j--) {
                        const node_t v2 = relabel_nodes[j];
                        V[v2].subtag = assigned_subtags[j%MIN_GROUP_SIZE];
                        //++cnt.subtag;
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

/*the removed group are not recycled. The are maxmum 2n groups*/
void ParCM::CoreMaint::OrderDelete(node_t x) {
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


/*insert all nodes in y at hte head of k when insertting edges.
* each group includes at most 16 (MAX_GROUP_SIZE) nodes, the relabel is not triggered*/
void ParCM::CoreMaint::OrderInsertHeadBatch(core_t k, vector<node_t> &y){
    if (y.empty()) return;
    cnt_ominsert+=y.size();

#if 1
    // lock head and head.next in order, 
    // all x in y not locked since the edgeinsert is not allowed this happen
    V[H[k]].OMLock(); 
    node_t xnext = V[H[k]].next; 
    V[xnext].OMLock();

    G[Hg[k]].OMLock(); 
    group_t gnext = G[Hg[k]].next;
    G[gnext].OMLock(); 

    const group_t g_head = G[Hg[k]].next;
    vector<group_t> groups;
    for(size_t group_id = group_size; group_id < group_size + y.size(); group_id++) {
        assert(group_id < 2 * n); // at most 2n groups.
        groups.push_back(group_id);
    }
    group_size += y.size();
    MultiTopListInsert(Hg[k], groups);

    tag_t tag0;
    if (g_head == Tg[k]) { // head = tail top list is empty
        tag0 = MAX_TAG/2;
    } else {
        tag0 = G[g_head].tag;
    }

    for (int i = 0; i < y.size(); i++){
        
        ASSERT(tag0 > (groups.size() - i) * INIT_TAG_GAP);

        G[groups[i]].tag = tag0 - (groups.size() - i) * INIT_TAG_GAP;
        V[y[i]].group = groups[i]; G[groups[i]].size++;
        V[y[i]].subtag = MAX_SUBTAG/2;
    }
    MultiListInsert(H[k], y);

    V[H[k]].OMUnlock(); V[xnext].OMUnlock(); 
    G[Hg[k]].OMUnlock(); G[gnext].OMUnlock();

#else // compact all to a group.
    const group_t g_head = G[Hg[k]].next;
    const node_t n_head = V[H[k]].next;
    const tag_t tag0 = G[g_head].tag;
    
    int group_num = y.size()/MAX_GROUP_SIZE + 1;
    
    //insert groups
    groups.clear();
    for(size_t group_id = group_size; group_id < group_size + group_num; group_id++) {
        assert(group_id < 2 * n); // at most 2n groups.
        groups.push_back(group_id);
    }
    group_size += group_num;
    MultiTopListInsert(Hg[k], groups); // insert to group head

    assert(tag0 > INIT_TAG_GAP * groups.size());
    size_t j = 0; 
    for(size_t i = 0; i < groups.size(); i++) {
        const group_t g = groups[i];
        G[g].tag = tag0 - (groups.size() - i) * INIT_TAG_GAP;
        for(size_t j = 0; j < MAX_GROUP_SIZE; j++) {
            size_t idx = i * MAX_GROUP_SIZE + j;
            if (idx >= y.size()) break;
            V[y[idx]].group = g; G[g].size++;
            V[y[idx]].subtag = j * (MAX_SUBTAG/MAX_GROUP_SIZE); ++cnt_tag;
        }
    }
    MultiListInsert(H[k], y);
#endif 
}

/*insert new nodes to tail without relabeling*/
void ParCM::CoreMaint::OrderInsertTail(core_t k, node_t y) {
    ++cnt_ominsert;
#if 1
DOLOCK:
    // lock tail.pre and tail in order.
    node_t xpre = V[T[k]].pre;
    V[xpre].OMLock(); 
    if (xpre != V[T[k]].pre) { // in case xpre is changed.
        V[xpre].OMUnlock();
        goto DOLOCK;
    } 
    V[T[k]].OMLock();

    group_t gpre = V[xpre].group;
    G[gpre].OMLock();
    G[Tg[k]].OMLock();
    

    //simple add to tail, fast and easy
    const group_t group_id = group_size;
    
    //groups.push_back(group_id);
    vector<group_t> groups; groups.push_back(group_id);

    group_size++;
    assert(group_size <= 2*n);
    TopListInsert(G[Tg[k]].pre, group_id);

    const group_t g_tail = G[Tg[k]].pre;
    tag_t tag0;
    if (g_tail == Hg[k]) { // tail = head top list is empty
        tag0 = MAX_TAG/2;
    } else {
        tag0 = G[g_tail].tag;
    }
    G[group_id].tag = tag0 + INIT_TAG_GAP;
    V[y].group = group_id; G[y].size++;
    V[y].subtag = MAX_SUBTAG/2; ++cnt_tag;

    V[xpre].OMUnlock(); V[T[k]].OMUnlock(); 
    G[gpre].OMUnlock(); G[Tg[k]].OMUnlock(); 

#else // compact the tail for inserting more nodes. 
    const group_t g_tail = G[Tg[k]].pre;
    const node_t n_tail = V[T[k]].pre;
    if (n_tail == H[k]) { // list is empty
        /*insert a new group*/
        const group_t group_id = group_size;
        groups.push_back(group_id);
        group_size++;
        TopListInsert(g_tail, group_id);
        V[y].group = group_id; G[y].size++;
        V[y].subtag = MAX_SUBTAG/MAX_GROUP_SIZE; ++cnt_tag;
        return;
    }
    const subtag_t subtag0 = V[n_tail].subtag;
    if (MAX_SUBTAG - subtag0 > MAX_SUBTAG/MAX_GROUP_SIZE) {
        ListInsert(n_tail, y);
        V[y].group = g_tail; G[y].size++;
        V[y].subtag = subtag0 + MAX_SUBTAG/MAX_GROUP_SIZE;
    } else {
        /*insert a new group*/
        const group_t group_id = group_size;
        groups.push_back(group_id);
        group_size++;
        TopListInsert(g_tail, group_id);
        V[y].group = group_id; G[y].size++;
        V[y].subtag = MAX_SUBTAG/MAX_GROUP_SIZE; ++cnt_tag;
    }
#endif
}


void ParCM::CoreMaint::OrderInsertTailBatch(core_t k, vector<node_t> &y, vector<node_t> &groups) {
    if (y.empty()) return;
    cnt_ominsert+=y.size();

#if 1  
    // load tail.pre and tail in order. 
    // all y insert after the tail with y.size groups.
DOLOCK:
    // lock tail.pre and tail in order.
    node_t xpre = V[T[k]].pre;
    V[xpre].OMLock(); 
    if (xpre != V[T[k]].pre) { // in case xpre is changed.
        V[xpre].OMUnlock();
        goto DOLOCK;
    } 
    V[T[k]].OMLock();

    group_t gpre = V[xpre].group;
    G[gpre].OMLock();
    G[Tg[k]].OMLock();
    
    const group_t g_tail = G[Tg[k]].pre;
    for(size_t group_id = group_size; group_id < group_size + y.size(); group_id++) {
        ASSERT(group_id < G.size()); // at most 2n groups.
        groups.push_back(group_id);
    }
    group_size += y.size();
    MultiTopListInsert(g_tail, groups);

    tag_t tag0;
    if (g_tail == Hg[k]) { // head = tail top list is empty
        tag0 = MAX_TAG/2;
    } else {
        tag0 = G[g_tail].tag;
    }


    for (int i = 0; i < y.size(); i++){
        G[groups[i]].tag = tag0 + (i + 1) * INIT_TAG_GAP; ;
        V[y[i]].group = groups[i]; G[groups[i]].size++;
        V[y[i]].subtag = MAX_SUBTAG/2;
    }
    MultiListInsert(V[T[k]].pre, y);

    V[xpre].OMUnlock(); V[T[k]].OMUnlock(); 
    G[gpre].OMUnlock(); G[Tg[k]].OMUnlock(); 

#endif 
}



int ParCM::CoreMaint::ParallelEdgeInsert(vector<pair<int, int>> &edges, int m, int m2, bool ispar){
    if(ispar) {
    
    #pragma omp parallel 
    {
        int workerid = omp_get_thread_num();
        const size_t reserve_size = n;

        PRIORITY_Q PQ(reserve_size);
        QUEUE R(reserve_size);
        vector<node_t> Vblack; Vblack.reserve(reserve_size);
        vector<node_t> Vcolor; Vcolor.reserve(reserve_size);

        #pragma omp for // schedule(dynamic, step)
        for (int i = m - 1; i >= m2; --i) {
            node_t u = edges[i].first;
            node_t v = edges[i].second;
            EdgeInsert(u, v, PQ, R, Vblack, Vcolor, workerid);
            if(3 == g_debug){
                //printf("**************p %d: parallel insert (%d, %d)************\n", workerid, u, v);
                //CheckLock(true);
                CheckOM(true);
            }
        }
    }
    return 0;
    
    } 
}


/*
* the Vblack is all black nodes. Vcolor is for all locked nodes.  */
int ParCM::CoreMaint::EdgeInsert(node_t u, edge_t v, PRIORITY_Q &PQ, QUEUE &R,
        vector<node_t> &Vblack, vector<node_t> &Vcolor, int p) {
    
    if (11728 == u && 17444 == v) {
        printf("debug");
    }
DOLOCK:
    // lock in order by assuming order u < v
    if (!Order(u, v)) std::swap(u, v);
    V[u].Lock(); V[u].islock = p;
    if (!Order(u, v)) {
        V[u].islock = NONE; V[u].Unlock(); 
        goto DOLOCK;
    }
    V[v].Lock(); V[v].islock = p;

    const core_t K = core[u];

    /*add edges*/
    graph.insert(u, v);
    V[v].islock = NONE; V[v].Unlock();

    #pragma omp atomic update
    V[u].degout++;

    if (V[u].degout <= K ) {
        #if 0 // update mcd
        if (core[u] <= core[v]) V[u].mcd++;
        if (core[v] <= core[u]) V[v].mcd++;
        #endif
        V[u].islock = NONE; V[u].Unlock();

        #pragma omp atomic update 
        Vcount[0]++;
        
        return 1; // update mcd and return
    }


    PQ.clear(); Vblack.clear(); Vcolor.clear();

    /* u is locked, all nodes in PQ are locked.*/
    Vcolor.push_back(u); 
    PQ.push(DATA(u, getTag(u), getSubtag(u), g_tag_version[core[u]])); 
    V[u].inQ = true;

    /******* traverse with topological order***********/
    while (!PQ.empty()) {
        node_t w = PQ.front().id; PQ.pop();//w has the smallest order

        if (unlikely(false == V[w].inQ) ) { continue; } //not in PQ, this is update the value of PQ.
        else { V[w].inQ = false; }
        
        if (V[w].degin + V[w].degout > K) { // forward
            Forward(w, PQ, Vcolor, p);
        } else if (V[w].degin > 0){ // backward (w must has degin > 0, or w can't be in PQ)
            Backward(w, R, p);
        }

    }

    /**ending phase**/

    /*update core number and mcd for BLACK nodes*/
    for (const node_t w: Vcolor){
        if (likely(V[w].color == BLACK)) { // in V*
            OrderDelete(w);// remove w from the order list O_k
            
            #pragma omp atomic update
            core[w]++;

            Vblack.push_back(w);
        }
    }
#if 0
    //insert Vb2 (ordered black nodes) to the head of O_k+1
    // MultiOrderInsert(H[K+1], Vblack);
    node_t pos_x = H[K+1];
    for(const node_t v: Vblack) {
        OrderInsert(pos_x, v); pos_x = v;
    }
#else 
    OrderInsertHeadBatch(K+1, Vblack);
#endif


/*update the mcd, our new method don't 
since the mcd is generated online*/
#if 0
    if (Vblack.empty()) {
        if (core[u] <= core[v]) V[u].mcd++;
        if (core[v] <= core[u]) V[v].mcd++;
    } else {
        for (const node_t u: Vblack) {
            deg_t mcd = 0;
            //for (const node_t v: graph[u]) {
            for (int i = graph.begin[u]; i < graph.end[u]; i++){
                node_t v = graph.node_idx[i];
                if (core[u] <= core[v]) mcd++;
                //u.core +1 for v. u->v
                if (BLACK != V[v].color && core[u] == core[v]) {
                    V[v].mcd++;
                }
            }
            V[u].mcd = mcd;
        }
    }
#endif 

    //reset all color
    for(const node_t u: Vcolor) {
        V[u].color = WHITE; V[u].degin = 0;
        V[u].islock = NONE; V[u].Unlock();
    }

    #pragma omp atomic update
    Vcount[Vcolor.size()]++;
    
    // if (3 == g_debug) {    
    //     printf("Vcolor size : %u:  ", Vcolor.size());
    //     // for(const node_t u: Vcolor) {
    //     //     printf(" %d,", u);
    //     // }
    //     printf("\n");
    // }
    return 0;
}

// bool ParCM::PRIORITY_Q::push(DATA d) {
// #if 0
//     if (v.empty()) version = d.version; //init version
//     v.push_back(d); 
//     if(d.version != version) return false;
//     push_heap(v.begin(), v.end(), CompareOrder());
//     return true;
// #else
//     v.push_back(d); 
//     push_heap(v.begin(), v.end(), CompareOrder());
//     return true;
// #endif
// }

void ParCM::PRIORITY_Q::enqueue(DATA &d) {

}

void ParCM::PRIORITY_Q::dequeue() {

} 

void ParCM::PRIORITY_Q::update(vector<Node> &V) {

}

/*
*Vcolor: all locked nodes*/
void ParCM::CoreMaint::Forward(node_t w, PRIORITY_Q &PQ, vector<node_t> &Vcolor, int p){
    V[w].color = BLACK; cnt_Vs++; cnt_Vp++;// core number can be update 
    
    //for (const edge_t w2: graph[w]) { // w -> w2
    // for (int idx = graph.begin[w]; idx < graph.end[w]; idx++){
    //     node_t w2 = graph.node_idx[idx];
    GRAPH_EDGE(w, w2)
        cnt_S++;
        if (SameCoreOrder(w, w2)) { // outgoing edges
            
            // ASSERT(getTag(w) < getTag(w2) || (
            //     getTag(w) == getTag(w2) && getSubtag(w) < getSubtag(w2) ));
            
            // if (p != V[w2].islock) { // avoid repeated lock 
            //     /*the order (within core number) cannot change, avoid dead-lock*/
            //     V[w2].Lock(); V[w2].islock = p; // lock with same core and then add to PQ.
            //     if (atomic_read(&core[w]) != atomic_read(&core[w2])) {
            //         V[w2].islock = NONE; V[w2].Unlock(); 
            //         continue;
            //     }
            //     Vcolor.push_back(w2); // it has all locked nodes.
            // }



            Vcolor.push_back(w2);
           
            if (!V[w2].inQ) { // w2 is not in PQ
                /*the concurrent PQ based on OM data structure. */
                DATA &data = DATA(w2, getTag(w2), getSubtag(w2), g_tag_version[core[w2]]);
                if(!PQ.push(data)) {
                    PQ_update(core[w]);
                }
                V[w2].inQ = true; 
            }

            #pragma omp atomic update
            V[w2].degin++;  // w2 is black
        }
    GRAPH_EDGE_END
}


/*update all nodes in PQ with new tag version
* we can optimize not use v2. 
* this is also work for parallel version*/
void ParCM::CoreMaint::PQ_update(const core_t k){
    PQ.v2 = PQ.v;
AGAIN:
    PQ.v.clear();
    id_t cur_version = EMPTY;
    for(DATA d: PQ.v2){
        if(EMPTY == cur_version) cur_version = g_tag_version[k];
        if(cur_version != g_tag_version[k]) {
            PQ.v.clear();
            goto AGAIN;
        }
        node_t u = d.id;
        PQ.v.push_back(DATA(u, getTag(u), V[u].subtag, g_tag_version[k]));

    }
    if(cur_version != g_tag_version[k]) {
        PQ.v.clear();
        goto AGAIN;
    }

    PQ.init(); // re-build the heap for priority queue.
}

void ParCM::CoreMaint::Backward(node_t w, QUEUE &R, int p) { //???
    V[w].color = GRAY; cnt_Vp++; // white -> gray. 
    
    node_t prenode = w;

    DoPre(w, R, p); //first time do black pre only;
    V[w].degout += V[w].degin; V[w].degin = 0;

    V[w].color = WHITE;

    while (!R.empty()) { // Q has nodes from BLACK to GRAY 
        node_t u = R.top(); R.pop(); V[u].inR = false; 
        V[u].degout += V[u].degin; V[u].degin = 0;
        V[u].color = GRAY;
        cnt_Vs--; //V* black -> gray
        
        DoAdj(u, R, p);
        OrderDelete(u); 
        OrderInsert(prenode, u); prenode = u;
        
        V[u].color = WHITE;

        ++cnt_ominsert_mid;
    }     
}

/*We only propagate the balck nodes that is in current Vblack*/
void ParCM::CoreMaint::DoPre(node_t u, QUEUE &R, int p) {
    //for (const edge_t v: graph[u]) { 
    // for (int i = graph.begin[u]; i < graph.end[u]; i++){
    //     node_t v = graph.node_idx[i];

    GRAPH_EDGE(u, v)
        if (SameCoreOrder(v, u) && BLACK == V[v].color && p == V[v].islock) { //u<-v  pre

            #pragma omp atomic update
            V[v].degout--;
            //ASSERT(V[v].degout >= 0);

            if (BLACK == V[v].color  && V[v].degin + V[v].degout <= core[v]) {
                R.push(v); V[v].inR = true;
            }

        }
    GRAPH_EDGE_END
}


void ParCM::CoreMaint::DoAdj(node_t u, QUEUE &R, int p) {
    //for (const edge_t v: graph[u]) { 
    // for (int idx = graph.begin[u]; idx < graph.end[u]; idx++){
    //     node_t v = graph.node_idx[idx];
    GRAPH_EDGE(u, v)
        if (SameCoreOrder(v, u) && BLACK == V[v].color && p == V[v].islock) { //u<-v  pre
            
            #pragma omp atomic update
            V[v].degout--; 
            
            ASSERT(V[v].degout >= 0);
            
            if (!V[v].inR && BLACK == V[v].color && V[v].degin + V[v].degout <= core[v]) {
                /*V[v].color = GRAY; // black -> gray*/
                R.push(v); V[v].inR = true;
            }

        } else if (SameCoreOrder(u, v) && V[v].degin > 0 && p == V[v].islock) { // u->v post

            #pragma omp atomic update
            V[v].degin--;

            if (!V[v].inR && BLACK == V[v].color && V[v].degin + V[v].degout <= core[v]) {
                //V[v].color = GRAY; 
                R.push(v); V[v].inR = true;
            }
        }
    GRAPH_EDGE_END
    
}



/********************************
* edge remoal has different version.
*******************************/
int ParCM::CoreMaint::ParallelEdgeRemove(vector<pair<int, int>> &edges, 
        const int m, const int m2, bool ispar){
    if (ispar) {
    
    #pragma omp parallel 
    {
        int workerid = omp_get_thread_num();
        const size_t reserve_size = n;
        
        QUEUE R(reserve_size);
        vector<node_t> Vblack; Vblack.reserve(reserve_size);
        vector<node_t> groups; groups.reserve(reserve_size);
        
        vector<bool> A(reserve_size, false); // or  max edge size graph.get_max_edge_size()

        #pragma omp for // schedule(dynamic, step)
        for (int i = m-1; i >= m2; --i) {
        //for (int i = m2; i < m; ++i) {
            node_t u = edges[i].first; node_t v = edges[i].second;
            EdgeRemove(u, v, R, Vblack, A, groups);
 #if 0           
            if(g_debug) 
                printf("**************remove (%d, %d) by worker %d************\n", u, v, workerid);

            if (g_debug) {
                CheckAllMCDValue(true);
            }
#endif
        }
    }
    
    } 
    // else {
    //     const size_t reserve_size = 1024 * 16;
    //     QUEUE R(reserve_size);
    //     vector<node_t> Vblack; Vblack.reserve(reserve_size);

    //     for (int i = m2; i < m; ++i) {
    //         EdgeRemove(edges[i].first, edges[i].second, R, Vblack);
    //     }
    // }

    return 0;
}


int ParCM::CoreMaint::EdgeRemove(node_t u, edge_t v, QUEUE &R, 
            vector<node_t> &Vblack, vector<bool> &A, vector<node_t> &groups) {

    R.clear(); Vblack.clear();

    /*lock in order avoid dead-lock
    * here may has dead lock. The case is that when lock v, but u is in propagation to v. 
    * in this case, the deadlock happen. */

#if 0
    if (u < v) {
        V[u].Lock(); 
        V[v].Lock(); 
    } else {
        V[v].Lock(); 
        V[u].Lock(); 
    }
#else
RELOCK:
    if (u < v) {
        if ( V[u].TestLock() ) {
            if (!V[v].TestLock()) {
                V[u].Unlock();
                goto RELOCK;
            } 
        } else {
            goto RELOCK;
        }
        
    } else {
        if (V[v].TestLock() ) {
            if (!V[u].TestLock()) {
                V[v].Unlock();
                goto RELOCK;
            } 
        } else {
            goto RELOCK;
        }
    }
#endif


    const core_t coreu = atomic_read(&core[u]);
    const core_t corev = atomic_read(&core[v]);
    const core_t K = std::min(coreu, corev);

    CheckMCD(u, coreu, EMPTY); 
    CheckMCD(v, corev, EMPTY);
    graph.remove(u, v);

    if (corev > K) { 
        DoMCD(u, coreu, R, Vblack);  
        V[v].Unlock();
    }
    if (coreu > K) { 
        DoMCD(v, corev, R, Vblack);  
        V[u].Unlock();
    }

    if (coreu == corev) {
        DoMCD(u, coreu, R, Vblack);
        DoMCD(v, corev, R, Vblack);
    }

    while(!R.empty()) {
        node_t w = R.top(); R.pop(); V[w].inR = false;
     
        int repeat_num = 0;
REDO:
        atomic_sub(&V[w].s, 1);
        int a = 0;
        GRAPH_EDGE(w, w2)
			cnt_S++; // count all visited edges except calculate mcd.   
            /* make check and lock together to reduce the gap time.
            * w2 can be a busy checking the core number.*/
            
            ASSERT(a < A.size());

            if (0 == repeat_num) A[a] = false; // init A for the first time. 

            if (!A[a] && atomic_read(&core[w2]) == K) {

#if 0
                V[w2].Lock(); 
                if (atomic_read(&core[w2]) != K) {
                    V[w2].Unlock(); 
                    printf("lock %d has core != %d\n", w2, K);
                    continue; 
                }
#else 
                if (false == V[w2].LockWithCore(K, &core[w2])) continue;
#endif
                CheckMCD(w2, K, w); 
                DoMCD(w2, K, R, Vblack); 

                A[a] = true;
            }
            a++;
        GRAPH_EDGE_END
        atomic_sub(&V[w].s, 1);

        ++repeat_num;
        if (V[w].s > 0) {
            goto REDO;
        }
    }
    
    if (likely(0 == Vblack.size())) return 1;
    // append u to the tail of O_(K-1). 
    //OrderInsertTail(K-1, u);
    OrderInsertTailBatch(K-1, Vblack, groups);
    for(const node_t v: Vblack) {
        V[v].Unlock();
    }

    //printf("Vblack size: %d\n", Vblack.size());
    
    return 0;
}

/*u is locked, R is local*/
void ParCM::CoreMaint::DoMCD(node_t u, core_t k, QUEUE &R, vector<node_t> &Vblack) {

    if (!V[u].IsLocked()) {
        ASSERT( V[u].IsLocked() ); // u is locked   
    } 
    if (V[u].mcd  < k ) {
        ASSERT( V[u].mcd  >= k ); // u.mcd >= u.core
    }
    atomic_sub(&V[u].mcd, 1);
    
    if ( V[u].mcd  < k ) { // u.mcd
    
        while ( atomic_read(&V[u].s_lock) == false && cas(&V[u].s_lock, false, true)) {}
        atomic_write(&core[u], k-1); atomic_write(&V[u].s, 2);
        V[u].s_lock = false;


        ASSERT(!V[u].inR); // u not in R

        R.push(u); V[u].inR = true;
        atomic_write(&V[u].mcd, EMPTY);

        cnt_Vs++; // count V*

		OrderDelete(u);
        Vblack.push_back(u); 
    } else {
        V[u].Unlock();
    }
    
}

/*u is locked
* here has bugs, the mcd always like 113 (112 is correct), 1 is calculated more.
* I may omit one special case, e.g. edge removing.*/
void ParCM::CoreMaint::CheckMCD(node_t u, core_t k, node_t from_v) {
    if (!V[u].IsLocked()) {
        ASSERT( V[u].IsLocked() ); // u is locked   
    } 
    if (EMPTY != atomic_read(&V[u].mcd)) { return; }

    // if (19659 == u) {
    //     printf("debug\n");
    // }

    deg_t mcd = 0; 
    
    // for (int idx = graph.begin[u]; idx < graph.end[u]; idx++){
    //     node_t v = graph.node_idx[idx];
    
    GRAPH_EDGE(u, v)

        //lock
        while (atomic_read(&V[v].s_lock) == false && cas(&V[v].s_lock, false, true)) {}

        core_t k1 = atomic_read(&core[v]); int s1 = atomic_read(&V[v].s);
        if (atomic_read(&core[v]) >= k || (k - 1 == atomic_read(&core[v]) && atomic_read(&V[v].s) > 0)) {
    
            atomic_write(&V[u].s_lock, false); //unlock

            ++mcd; 
            if (k-1 == atomic_read(&core[v])) {
                if (atomic_read(&V[v].s) > 0) {
                    
                    if (v != from_v) {
                        //printf(" =====  the %d.s = %d\n", v, V[v].s);
                        if(atomic_read(&V[v].s) == 1) {
                            cas(&V[v].s, 1, 3);
                        }
                    }
                    if (atomic_read(&V[v].s) == 0) { 
                        --mcd;
                        printf(" ********  mcd-- with s=0");
                    }
                }
            }


        } else {
            atomic_write(&V[u].s_lock, false); //unlock
        } 

        core_t k2 = atomic_read(&core[v]); int s2 = atomic_read(&V[v].s);
        //if (k1 != k2) {printf(" ====== %d change core number %d -> %d\n", v, k1, k2); }
        //if (s1 != s2) {printf(" ====== %d change s %d -> %d\n", v, s1, s2); }
    GRAPH_EDGE_END

    if (!V[u].IsLocked()) {
        ASSERT( V[u].IsLocked() ); // u is locked   
    }
    //if (atomic_read(&V[u].mcd) == EMPTY) {
        ASSERT(atomic_read(&V[u].mcd) == EMPTY);
    //}
    atomic_write(&V[u].mcd, mcd);
    
}


/****************************for debug test**********************************/
bool ParCM::CoreMaint::CheckLock(bool info) {
    bool r = true;
    for (node_t u = 0; u < n; u++) {
        if (OMP_LOCK == g_lock_type) {
            auto &a = V[u].lock_omp._x;
            if (a[0] != 0 || a[1] != 0 || a[2] != 0 || a[3] != 0) {
                r = false;
                if(info) printf("Error: %d is still locked!!!\n", u);
            }

            if (V[u].islock != NONE) {
                r = false;
                if (info) printf("Error: %d islock %d !!!\n", u, V[u].islock);
            }
        }
    }
    return r;
}

vector<node_t> ParCM::CoreMaint::InitTestOM(size_t insert_size) {
     vector<int> nums(H.size(), 0);
    // for(size_t k = 0; k < H.size(); k++){
    //     for(node_t v = V[H[k]].next; v != T[k]; v = V[v].next) {
    //         nums[k]++;
    //     }
    // }

    // core_t test_k = -1;
    // for(size_t k = 0; k < nums.size(); k++) {
    //     if (nums[k] >= insert_size) test_k = k;
    // }

    // if(-1 == test_k) {
    //     printf("cannot find a list has %d size!\n", insert_size);
    //     return 1;
    // }
    const core_t test_k = 1;

    printf("k=%d insert %ld / %d nodes!\n", test_k, insert_size, nums[test_k]);
    node_t h = H[test_k]; node_t t = T[test_k];
    vector<node_t> nodes; int size = 0;
    node_t upos; 
    for(node_t u = V[h].next; u != t; u = V[u].next) {
        nodes.push_back(u);
        upos = V[u].next;
        size++; 
        if (size >insert_size) break;
    }
    return nodes;
}

/*by inserting a large number of nodes to triger relabeling.
* remove tail.pre and insert after head.next. */
int ParCM::CoreMaint::TestOM(vector<node_t> &nodes) {
   
    const core_t test_k = 1;
    node_t h = H[test_k];

    // test sequential OM
#if 0
    //int workerid = omp_get_thread_num();
    for (int i = 0; i < nodes.size(); i++) {
        OrderDelete(nodes[i]);
    }

    for(int i = 0; i < nodes.size(); i++) {
        //printf("w %d insert %d after %d\n", workerid, nodes[i], h);
        OrderInsert(h, nodes[i]);
    }
#else   
#pragma omp parallel
{
    int workerid = omp_get_thread_num();
    #pragma omp for
    for (int i = 0; i < nodes.size(); i++) 
    {
        //OrderDelete(nodes[i]);
        Order(nodes[1], nodes[2]);
    }

    // #pragma omp for
    // for(int i = 0; i < nodes.size(); i++) 
    // {
    //     //printf("w %d insert %d after %d\n", workerid, nodes[i], h);
    //     OrderInsert(h, nodes[i]);
    // }
}
//     #pragma omp parallel 
//     {
//         int workerid = omp_get_thread_num();
//         #pragma omp for
//         for (int i = 0; i < nodes.size(); i++) {
//             OrderDelete(nodes[i]);
//         }
// #if 1
//         #pragma omp for
//         for(int i = 0; i < nodes.size(); i++) {
//             printf("w %d insert %d after %d\n", workerid, nodes[i], upos);
//             OrderInsert(upos, nodes[i]);
//         }
// #en

// #if 0 // the vector cost much time
//         #pragma omp for
//         for(int i = 1; i < nodes.size(); i++) {
//             if (i % 100 == 0) {
//                 vector<node_t> y; y.reserve(16); 
//                 for(int j = i - 90; j <= i; j++)
//                     y.push_back(nodes[i]);
//                 OrderInsertHeadBatch(test_k, y);
//                 printf("batch insert %d\n", i);
//             }
//         }
// #endif 
        // #pragma omp for
        // for(size_t i = 0; i < insert_size; i++) {
        
        // #if 0  // test insert to head
        //         vector<node_t> y;
        //         for (int i = 0; i < 17; i++) {
        //             node_t v = V[t].pre;
        //             OrderDelete(v);
        //             y.push_back(v);
        //         }
        //         OrderInsertHeadBatch(test_k, y);
        // #endif
        //         node_t v = V[t].pre;
        //         OrderDelete(v);
        //         printf("delete %d\n", v);
        //         //OrderInsertTail(test_k, v);
        //         //OrderInsert(h, v);
        //         // insert at head
                
            
        //         // insert after tail
        // }
    
#endif    
    // for(node_t u = V[h].next; u != t; u = V[u].next) {
    //     nodes.push_back(u);
    // }
    
    return 0;
}

/*after test, check the k-order maintenance*/
int ParCM::CoreMaint::CheckOM(bool info) {
    //check reverse order list  
    size_t num1 = 0, num2 = 0; 
    for(size_t k = 0; k < H.size(); k++){
        vector<node_t> order, rorder;
        for(node_t v = V[H[k]].next; v != T[k]; v = V[v].next) {
            if (v == -1) {
                assert(v != -1);
            }
            order.push_back(v); num1++;
            if (!Order(v, V[v].next)) {
                printf("first: %d precequeu %d wrong!", v, V[v].next);
            }
        }

        for(node_t v = V[T[k]].pre; v != H[k]; v = V[v].pre) {
            assert(v != -1);
            rorder.push_back(v); num2++;
            if (!Order(v, V[v].next)) {
                printf("second: %d precequeu %d wrong!", v, V[v].next);
            }
        }
        std::reverse(rorder.begin(), rorder.end());

        ASSERT(order == rorder);
    }
    if (num1 != num2 || num1 != n) {
        if(info) printf("wrong! *** total number in sublist: %ld but correct %ld\n", num1, n);
        else ASSERT(false);

    }

    num1 = num2 = 0;
    for(int k = 0; k < H.size(); k++){
        vector<node_t> order, rorder;
        for(group_t g = G[Hg[k]].next; g != Tg[k]; g = G[g].next) {
            assert(g != -1);
            order.push_back(g); num1++;
        }

        for(node_t g = G[Tg[k]].pre; g != Hg[k]; g = G[g].pre) {
            assert(g != -1);
            rorder.push_back(g); num2++;
        }
        std::reverse(rorder.begin(), rorder.end());

        ASSERT(order == rorder);
    }
    
    // group_size includes non-recycled nodes.
    if (num1 != num2) {
        if(info) printf("wrong! *** total number in topgrouplist: %ld but correct %ld\n", num1, n);
        else ASSERT(false);

    }


    //check order list core number 
    for(int k = 0; k < (int)H.size(); k++){
        for(node_t v = V[H[k]].next; v != T[k]; v = V[v].next) {
            if (k != core[v]) {

            if(info) printf("wrong! *** list check %d has core number %d but should be %d\n", v, k, core[v]);
            else ASSERT(false);
            
            }
        }
    }


    //check tag
    for(core_t k = 0; k < H.size(); k++){
        for(node_t v = V[H[k]].next; v != T[k]; v = V[v].next) {
            node_t next = V[v].next;
            if (V[v].group == V[next].group) {
                if (! (getSubtag(v) < getSubtag(next)) ) {
                    if(info)  printf("k=%d, %d.subtag: %u, next %d.subtag: %u\n",
                        k, v, V[v].subtag, next, V[next].subtag);
                    else ASSERT(false);
                }
            } else {
                if (! (getTag(v) < getTag(next)) ) {
                    if(info)  printf("k=%d, %d.subtag: %llu, next %d.subtag: %llu\n",
                        k, v, getTag(v), next, getTag(next));
                    else ASSERT(false);
                }
            }

        }

        // for(node_t g = G[Hg[k]].next; g != Tg[k]; g = G[g].next) {
        //     node_t next = G[g].next;
        //     if (G[g].tag >= G[next].tag) {
        //         if(info)  printf("k=%d, %d.tag: %llu, next %d.tag: %llu\n",
        //             k, g, G[g].tag, next, G[next].tag);
        //         else ASSERT(false);
        //     }

        // }

        // if (V[H[k]].tag != 0) {
        //     node_t next = V[H[k]].next;
        //     if(g_debug)  printf("%d head.tag: %llu, next  %d.tag: %llu\n",
        //         k, V[H[k]].
        //         , next, V[next].tag);
        //     else ASSERT(false);

        // }
        // if(V[T[k]].tag != MAX_TAG) {
        //     ASSERT(false);
        // }
        
    }
}

int ParCM::CoreMaint::CheckAllMCDValue(bool info) {
    //printf("begin check mcd value: with n = %d\n", n);
    for (node_t i = 0; i < n; i++){
		if (EMPTY != atomic_read(&V[i].mcd)) {

            deg_t mcd = 0;
            
            for (node_t idx = (graph.begin[i]); \
                idx < atomic_read(&graph.end[i]); idx++){ \
                node_t v = atomic_read(&graph.node_idx[idx]);

                ASSERT(v < n);
                if (atomic_read(&core[i]) <= atomic_read(&core[v])) {
                    mcd++;
                }
            }

            if (mcd != atomic_read(&V[i].mcd)) { //check mcd
                if(info) {
                    printf("%d.mcd = %d, should be %d\n", i, atomic_read(&V[i].mcd), mcd);
                } else {
                    ASSERT(false);
                }
            }
        }
    }
    return 0;
}

/*debug: check for result*/
int ParCM::CoreMaint::CheckCore(vector<core_t> &tmp_core, bool info){
    int r = 0;
    for(node_t i = 0; i < n; i++) {
        if(tmp_core[i] != core[i]){
            r = 1;
            if(info) printf("wrong! *** %d: core is %d, but correct is %d \n", i, core[i], tmp_core[i]);
        }
    }
    return r;
}

int ParCM::CoreMaint::CheckDeg(bool info) {
    //check degout degin core mcd
    for (node_t u = 0; u < n; u++){
        // degout <= core number
        if (V[u].degout > core[u]) {
            if(info)  printf("degout < core:  %d.degout = %d, core[%d]=%d\n",
                 u, V[u].degout, u, core[u]);
            else ASSERT(false);
        }

        // check the degout and degin
        if (0 != V[u].degin){
            if(info)  printf("degin:  %d.degin = %d, should be %d\n",
                 u, V[u].degin, 0);
            else ASSERT(false);
        }
        deg_t degout = 0, mcd = 0;

        //for (const node_t v: graph[u]){
        GRAPH_EDGE(u, v)
            if (Order(u, v)) degout++;
            //if (core[u] <= core[v]) mcd++;
        GRAPH_EDGE_END

        if (degout != V[u].degout) {
            if(info) printf("degout:  %d.degout = %d, should be %d\n",
                u, V[u].degout, degout);
            else ASSERT(false);
        }
        // if (mcd != V[u].mcd) { //check mcd
        //     if(g_debug) printf("EDGE [%d (%d, %d)]:  %d.mcd = %d, should be %d\n",
        //         id, x, y, u, V[u].mcd, mcd);
        //     else ASSERT(false);
        // }
    }
    return 0;
}


void ParCM::CoreMaint::PrintOMVersion() {
    int num = 0;
    for(core_t k = 0; k < H.size(); k++){
        if(g_tag_version[k] > 0) {
            num++;
            printf("k%d:%d\t", k, g_tag_version[k]);
            if (10 == num) {
                printf("\n"); 
                num = 0;
            }
        }
    }
    printf("\n");
}




#if 0
/*x, y: inserted edge
*id: the id of inserted edge
* check each time to find bug*/
int ParCM::CoreMaint::CheckAll(node_t x, node_t y, int id, vector<core_t> &tmp_core, vector<node_t> &order_v) {

    if(g_debug) printf("*********** Our Check **************\n");

        //check core number
    if(g_debug) {
        for(int i = 0; i < n; i++) {
            if(tmp_core[i] != core[i]){
                printf("wrong! *** %d: core is %d but %d after I/R %d edges \n", i, core[i], tmp_core[i], id);
            }
        }
    } else ASSERT_INFO(tmp_core == core, "wrong result after insert");

    if(0 == g_debug)    return 0; // only check the core number for release version. 

    int error = 0;
    //check reverse order list  
    size_t total_ver_num = 0;
    for(size_t k = 0; k < H.size(); k++){
        vector<node_t> order, rorder;
        for(node_t v = V[H[k]].next; v != T[k]; v = V[v].next) {
            order.push_back(v);
            total_ver_num++;
        }

        for(node_t v = V[T[k]].pre; v != H[k]; v = V[v].pre) {
            rorder.push_back(v);
        }
        std::reverse(rorder.begin(), rorder.end());

        ASSERT(order == rorder);
    }
    if (total_ver_num != n) {
        if(g_debug) printf("wrong! *** total number in list: %ld but %ld", total_ver_num, n);
        else ASSERT(false);

    }


    //check order list core number 
    for(int k = 0; k < (int)H.size(); k++){
        for(node_t v = V[H[k]].next; v != T[k]; v = V[v].next) {
            if (k != core[v]) {

            if(g_debug) printf("wrong! *** list check %d has core number %d but should be %d\n", v, k, core[v]);
            else ASSERT(false);
            
            }
        }
    }


    //check tag
    // vector<int> order_v_index; int order = 0;
    // for (const node_t v: order_v) { // order_v of order
    //     order_v_index.push_back(order++);
    // } 

    for(size_t k = 0; k < H.size(); k++){
        for(node_t v = V[H[k]].next; v != T[k]; v = V[v].next) {
            node_t next = V[v].next;

            // if (next != T[k] && order_v_index[v] >= order_v_index[next]) {
            //     ASSERT(false);
            // }

            if (V[v].tag >= V[next].tag) {

                if(g_debug)  printf("k=%d, %d.tag: %llu, next %d.tag: %llu\n",
                    k, v, V[v].tag, next, V[next].tag);
                else ASSERT(false);

            }

        }

        if (V[H[k]].tag != 0) {
            node_t next = V[H[k]].next;
            if(g_debug)  printf("%d head.tag: %llu, next  %d.tag: %llu\n",
                k, V[H[k]].
                , next, V[next].tag);
            else ASSERT(false);

        }
        if(V[T[k]].tag != MAX_TAG) {
            ASSERT(false);
        }
        
    }

    //check degout degin core mcd
    for (size_t u = 0; u < n; u++){
        // degout <= core number
        if (V[u].degout > core[u]) {
            if(g_debug)  printf("degout [%d (%d, %d)]:  %d.degout = %d, core[%d]=%d\n",
                id, x, y, u, V[u].degout, u, core[u]);
            else ASSERT(false);
        }

        // check the degout and degin
        if (0 != V[u].degin){
            if(g_debug)  printf("degin [%d (%d, %d)]:  %d.degin = %d, should be %d\n",
                id, x, y, u, V[u].degin, 0);
            else ASSERT(false);
        }
        deg_t degout = 0, mcd = 0;

        //for (const node_t v: graph[u]){
        for (int idx = `[u]; idx < graph.end[u]; idx++){
            node_t v = graph.node_idx[idx];

            if (Order(u, v)) degout++;
            if (core[u] <= core[v]) mcd++;
        }
        if (degout != V[u].degout) {
            if(g_debug) printf("EDGE [%d (%d, %d)]:  %d.degout = %d, should be %d\n",
                id, x, y, u, V[u].degout, degout);
            else ASSERT(false);
        }
        if (mcd != V[u].mcd) { //check mcd
            if(g_debug) printf("EDGE [%d (%d, %d)]:  %d.mcd = %d, should be %d\n",
                id, x, y, u, V[u].mcd, mcd);
            else ASSERT(false);
        }
    }

    return error;

}
#endif 
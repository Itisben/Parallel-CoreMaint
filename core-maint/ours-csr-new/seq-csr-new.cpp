#include "def.h"
#include "seq-csr-new.h"
#include <algorithm>


//for counting 
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
/********SeqCM************************************/
SeqCM::CoreMaint::CoreMaint(size_t n, GRAPH &graph, vector<core_t> &core):
        n{n}, graph{graph}, core{core} {
    
    const size_t reserve_size = n; //2^15
    Vblack.reserve(reserve_size);
    Vcolor.reserve(reserve_size);
    
    R = QUEUE(reserve_size); //2^14
    PQ = PRIORITY_Q(reserve_size);
  
    /* init the Order List by k-order odv */
    max_core = 0;
    for(core_t k: core) {  // find max core
        if (k > max_core) max_core = k;
    }

    max_core += MAX_CORE_NUMBER; // the maximum core number after increasing. 

    Node node; 
    V = vector<Node>(n + 2*max_core + 2, node); // Vertices
    H = vector<node_t>(max_core + 1); // head of order list k
    T = vector<node_t>(max_core + 1); // tail of order list k

    // at most n group can be new inserted
    Group group; 
    G = vector<Group>(2 * n + 2*max_core + 2, group);
    group_size = n; // initially n vertices have n groups.
    Hg = vector<node_t>(max_core + 1); // head of group list k  
    Tg = vector<node_t>(max_core + 1); // tail of group list k

    g_tag_version = vector<ver_t>(max_core+1, 0); 
}



void SeqCM::CoreMaint::Init(vector<int> &odv) {
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
        // insert to boTgom list after tail.pre p.
        {  
            node_t &p = V[T[k]].pre; 
            V[p].next = v; 
            V[v].pre = p; V[v].next = T[k];
            V[T[k]].pre = v;
        }
        // insert to top list after tail.pre q;
        {
            group_t &p= G[Tg[k]].pre;
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
#if 0
    for (size_t u = 0; u < n; u++) {
        deg_t mcd = 0;
        //for (node_t v: graph[u]) {
        for (node_t idx = graph.begin[u]; idx < graph.end[u]; idx++){
            node_t v = graph.node_idx[idx];
            if (core[u] <= core[v]) mcd++;
        }
        V[u].mcd = mcd;
    }
#endif 

}

void SeqCM::CoreMaint::ComputeCore(std::vector<int>& core, std::vector<int>& order_v, bool with_init) {
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
    }
    /*init degout, here is more efficient that check the order of u and v*/
    if(with_init) {V[v].degout = rem;}
  }
}


inline void SeqCM::CoreMaint::ListDelete(node_t x) {
    node_t pre = V[x].pre; node_t next = V[x].next;
    V[pre].next = next; V[next].pre = pre;
    V[x].pre = V[x].next = NONE;
}

/*insert y after x in list*/
inline void SeqCM::CoreMaint::ListInsert(node_t x, node_t y) {
    node_t next = V[x].next;
    V[x].next = y; V[y].pre = x;
    V[y].next = next; V[next].pre = y;
}


/*insert all y after x in list*/
inline void SeqCM::CoreMaint::MultiListInsert(node_t x, vector<node_t> &y) {
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

inline void SeqCM::CoreMaint::TopListDelete(group_t x) {
    group_t pre = G[x].pre; group_t next = G[x].next;
    G[pre].next = next; G[next].pre = pre;
    //G[x].pre = G[x].next = NONE;
    G[x].pre = G[x].next = -11;  //for debug
}

inline void SeqCM::CoreMaint::TopListInsert(group_t x, group_t y) {
    group_t next = G[x].next;
    G[x].next = y; G[y].pre = x;
    G[y].next = next; G[next].pre = y;
}

/*insert all y after x in list*/
inline void SeqCM::CoreMaint::MultiTopListInsert(node_t x, vector<group_t> &y) {
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

inline bool SeqCM::CoreMaint::Order(node_t x, node_t y) {
    ++cnt_omorder;
    if (unlikely(core[x] == core[y]) ) {
        if(unlikely(V[x].group == V[y].group)){
            return (V[x].subtag < V[y].subtag);
        } else {
            return (G[V[x].group].tag < G[V[y].group].tag);
        }
    } else return core[x] < core[y];

}

bool SeqCM::CoreMaint::SameCoreOrder(node_t x, node_t y) {
    ++cnt_omorder;
    if (unlikely(core[x] == core[y]) ) {
        if(unlikely(V[x].group == V[y].group)) {
            return (V[x].subtag < V[y].subtag);
        } else {
            return (G[V[x].group].tag < G[V[y].group].tag);
        }
    } else return false;
}

/*insert y after x. x->z 
* each list has head and tail as empty nodes. 
* j = 1
* while wj <= j*j :
*   update wj
*   j++
*/
void SeqCM::CoreMaint::OrderInsert(node_t x, node_t y) {
    ++cnt_ominsert;
AGAIN:
    const group_t g0 = V[x].group;
    const node_t z = V[x].next; // x-> z 
    const subtag_t subtag0 = V[x].subtag;
    subtag_t subw;

    /*there are two cases:
    * 1. z and x in the same group. 
    * 2. z and x in the different group, y insert append to x within group*/
    if (V[x].group == V[z].group) {
        subw = V[z].subtag - subtag0;
    } else {
        subw = MAX_SUBTAG - subtag0;
    }
    
    /* it not has enough label space. The relabel is trigered. */ 
    if (subw < 2) {
        // we use global version to STM. 
        ++g_tag_version[core[x]]; // doing relabel, odd number, "lock" order list
        SimpleRelabel(x);
        ++g_tag_version[core[x]]; // finish, even number, "unlock" order list
        goto AGAIN;
    }

    V[y].subtag = (V[x].subtag + subw / 2);
    V[y].group = g0;
    ++G[g0].size; //if group size = 0, group is removed.

    cnt_tag++;
    ListInsert(x, y);
}

/* this is a simple implementation. does not affect the performance
* starting from x, all nodes (including x) that 
* after x in group are inserted as new groups
* the parallel version need AssignLabel algorithm*/
void SeqCM::CoreMaint::SimpleRelabel(node_t x) {
    const group_t g0 = V[x].group;
    const tag_t tag0 = G[g0].tag; 

    assert(tag0 <= (MAX_TAG - 10)); // reach the limit, all nodes require reassign the labels.

    // get the number of inserted groups.
    relabel_nodes.clear();
    for(node_t u = x; g0 == V[u].group; u = V[u].next) {
        relabel_nodes.push_back(u);
    }
    
    const size_t relabel_size = relabel_nodes.size(); // total number of relabed node in group.
    group_t g = G[g0].next;
    //assert(g != NONE);
    tag_t w = G[g].tag - tag0;
    size_t j = relabel_size + 1;

    /*relabel the group*/
    if ( unlikely(w <= relabel_size) ) { 
        //the insertion walks down the list until w > j^2
        //the gap after relabeling > j. To insert n groups, j has to > n (and double???)
        j = 1;
        while (w <= j * j || j <=  relabel_size) { // at least insert j nodes between two. 
        //while(w <= j * j * relabel_size) {
            g = G[g].next; 
            w = G[g].tag - tag0;
            j++;
        }

        // in case traverse to the tail, scope the tag to j^2 + 1 (and double ? we have enough sapce).
        if (w > INIT_TAG_GAP) w = min(INIT_TAG_GAP, j * j * 1024);

        //relabel the j-1 record s_1 to s_j-1 with the labels 
        tag_t taggap = w/j;
        
        int tag_num = 0;
        for (size_t k = 1, g = G[g0].next; k < j; k++, g=G[g].next) {
            // the wide gap w is splited into j parts by w/j, 
            G[g].tag = taggap * k  + tag0; 
            cnt_tag++;
            /************************ update PQ when tag changed********************/
            //if (V[z].inQ) {PQ.push(DATA(z, V[z].tag));}
            tag_num++;
        }
        if(tag_num > cnt_rebalance_max_tag) cnt_rebalance_max_tag = tag_num;

        //reset w after relabel
        g = G[g0].next;
        w = G[g].tag - tag0; 
    }

    /* in case x is the tail, scope the tag. To a reasonable range.
    * INIT_TAG_GAP / 32 is test value*/ 
    if (w > INIT_TAG_GAP) w = min(INIT_TAG_GAP, j * j * 1024);


    /*insert n group after g0, consider later. Ignore the group size is 0*/
    groups.clear();
    for(size_t group_id = group_size; group_id < group_size + relabel_size; group_id++) {
        assert(group_id < 2 * n); // at most 2n groups.
        groups.push_back(group_id);
    }
    group_size += relabel_size;

    MultiTopListInsert(g0, groups);   
    tag_t taggap = w / (relabel_size + 1);
    for (size_t k = 0; k < relabel_size; k++) {
        group_t g = groups[k];
        node_t u = relabel_nodes[k];
        G[g].tag = tag0 + taggap * (k + 1);
        G[g].size = 1;
        V[u].group = g;
        V[u].subtag = MAX_SUBTAG / 2;
        ++cnt_tag;
    }
    ++cnt_omsplit;
    
}

#if 0
 /* if the group tag is relabed, the group member is impossible to be notified. 
 *  so that the PQ is notified for update
 * one solution: 
 * 1 each relabel process has an version number v;
 * 2 all vertices in PQ must has same v;
 * 3 if one vertex has higher v, so that all the other vertices have to update.
 * 4 for parallel version, we have to update again if failed, until all vertices has same version number. 
 * 5 this is efficient, since relabel cost O(1) amotized time.
 * 6 As the relabel may affected many nodes, the PQ can be maintained 
 *   with a version number (add by 1 when each relabel happens). All nodes in PQ must has same version number,
 *   If failed, all nodes must update to the latest version. 
 * 6.1: ???? in parallel version, to reduce the frequent version update, one strategy is 
 *      to split all node into multiple parts (logn = 32)
 *  ** for this I don't have better idea.  
 * 
 * For implementation, we make a simple version for easy test.
 * lock all nodes after x within group, and split into 
 **/
/*this is simple implementation, for testing. We are the first one, without competitors.*/
void SeqCM::CoreMaint::DeepRelabel(node_t x) {
    //node_t z = V[x].next; // x-> z 
    const group_t g0 = V[x].group;
    const tag_t tag0 = G[g0].tag; 

    assert(tag0 == MAX_TAG-1); // reach the limit, all nodes require reassign the labels.

    //split g0 to n groups 
    size_t n = G[g0].size / MAX_GROUP_SIZE;
    
   
    group_t g = G[g0].next;
    tag_t w = G[g].tag - tag0;
    size_t j = 1;

    /*relabel the group*/
    if (unlikely(w <= n)) { 
        //the insertion walks down the list until w > j^2
        //the gap after relabeling > j. To insert n groups, j has to > n (and double???)
        while (w <= j * j && j <= 2 * n) { // could be geometric increasing. ???
            g = G[g].next; 
            w = G[g].tag - tag0;
            j++;
        }

        // in case traverse to the tail, scope the tag to j^2 + 1 (and double ? we have enough sapce).
        if (w > INIT_TAG_GAP) w = 2 * (j * j + 1); 

        //relabel the j-1 record s_1 to s_j-1 with the labels 
        tag_t taggap = w/j;
        
        for (size_t k = 1, g = G[g0].next; k < j; k++, g=G[g].next) {
            // the wide gap w is splited into j parts by w/j, 
            G[g].tag = taggap * k  + tag0; 
            cnt_tag++;
            
            /************************ update PQ when tag changed********************/
            //if (V[z].inQ) {PQ.push(DATA(z, V[z].tag));}
        }

        //reset w after relabel
        g = G[x].next;
        w = G[g].tag - tag0; 
    }

    // in case x is the tail, scope the tag.
    if (w > INIT_TAG_GAP) w = 2 * (j * j + 1);


    /*insert n group after g0, consider later. Ignore the group size is 0*/
    vector<group_t> y; size_t G_size = G.size();
    for(size_t group_id = G_size; group_id < G_size + n; group_id++) {
        Group group;  G.push_back(group);
        y.push_back(group_id);
    }
    MultiTopListInsert(g0, y);   
    tag_t taggap = w / (n + 1);
    for (size_t k = 1; k <= n; k++) {
        G[y[k]].tag = tag0 + taggap * k; 
        cnt_tag++;
    }
    
}
#endif // DeepRelabel


#if 0
void SeqCM::CoreMaint::MultiOrderInsert(node_t x, vector<node_t> &y) {
    const size_t size = y.size();
    if (0 == size) return 0 ;
    else if (1 == size) {
        OrderInsert(x, y[0]);
        return 1;
    }

    unsigned int j = 1;
    node_t z = V[x].next;
    tag_t tag0 = V[x].tag;
    tag_t w = V[z].tag - tag0;
    
    if (unlikely(w <= size)) {
        //find j: 1 - j-1
        while (w <= j * j * size) { // could be geometric increasing.
            z = V[z].next;
            w = V[z].tag - V[x].tag;
            j++;
        }

        /*scope the tag range*/
        if (w > INIT_TAG_GAP * size * j) w = INIT_TAG_GAP * size * j;

        //relabel
        tag_t tagbase = w/j;
        for (size_t k = 1, z = V[x].next; k < j; k++, z=V[z].next) {
            V[z].tag = tagbase * k + tag0; cnt_tag++;
            // push into PQ again when tage changed. 
            if (V[z].inQ) {PQ.push(DATA(z, V[z].tag));}
        }
        // reset label 
        z = V[x].next;
        w = V[z].tag - V[x].tag; //reset w
    }

    /*scope the tag range*/
    if (w > INIT_TAG_GAP * size) w = INIT_TAG_GAP * size;

    //update label
    w = w / (size + 1); 
    for (size_t i = 0; i < size; i++) {
        V[y[i]].tag = (V[x].tag + w * (i+1)); cnt_tag++;
        // push into PQ again when tag is changed. 
        if (V[y[i]].inQ) {PQ.push(DATA(y[i], V[y[i]].tag));}
    }
    
    MultiListInsert(x, y);
}
#endif 

/*the removed group are not recycled. The are maxmum 2n groups*/
void SeqCM::CoreMaint::OrderDelete(node_t x) {
    ++cnt_omdelete;

    const group_t g0 = V[x].group;
    --G[g0].size; // update group size
    assert(G[g0].size >= 0);
    if(0 == G[g0].size) {
        TopListDelete(g0); 
    }
    ListDelete(x); 
    V[x].group = EMPTY; // removed no group. 
}


/*insert all nodes in y at hte head of k when insertting edges.
* each group includes at most 16 (MAX_GROUP_SIZE) nodes, the relabel is not triggered*/
void SeqCM::CoreMaint::OrderInsertHeadBatch(core_t k, vector<node_t> &y){
    if (y.empty()) return;
    cnt_ominsert+=y.size();
#if 1
    const group_t g_head = G[Hg[k]].next;
    groups.clear();
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
        G[groups[i]].tag = tag0 - (groups.size() - i) * INIT_TAG_GAP; ;
        V[y[i]].group = groups[i]; G[groups[i]].size++;
        V[y[i]].subtag = MAX_SUBTAG/2;
    }
    MultiListInsert(H[k], y);

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
void SeqCM::CoreMaint::OrderInsertTail(core_t k, node_t y) {
    ++cnt_ominsert;
#if 1
    //simple add to tail, fast and easy
    const group_t group_id = group_size;
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

    
#else  // compact the tail for inserting more nodes. 
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

/*the insertion is the same
* only thing is how to maintain the PQ. PQ has tag and subtag as the key.
* I must find a good way. */
int SeqCM::CoreMaint::EdgeInsert(node_t u, edge_t v) {

    if(g_debug) printf("**************insert (%d, %d)************\n", u, v);

    if (214 == u && 1431 == v) {
        printf("debug");
    }
    
    // assuming order u < v
    if (!Order(u, v)) std::swap(u, v);
    const core_t K = core[u];

    /*add edges*/
    graph.insert(u, v);
    V[u].degout++;
    if (V[u].degout <= K ) {
        if (core[u] <= core[v]) V[u].mcd++;
        if (core[v] <= core[u]) V[v].mcd++;
        return 1; // update mcd and return
    }

    PQ.clear(); Vblack.clear(); Vcolor.clear();
    PQ.push(DATA(u, getTag(u), V[u].subtag, g_tag_version[core[u]])); 
    V[u].inQ = true;

    /******* traverse with topological order***********/
    while (!PQ.empty()) {
        node_t w = PQ.front().id; PQ.pop();//w has the smallest order
        if (unlikely(false == V[w].inQ) ) { continue; } //not in PQ, this is update the value of PQ.
        else { V[w].inQ = false; }
        if (w == 1771) {
            printf("debug 1771");
        }

        if (V[w].degin + V[w].degout > K) { // forward
            Forward(w);
        } else if (V[w].degin > 0){ // backward (w must has degin > 0, or w can't be in PQ)
            Backward(w);
        }

    }

    /**ending phase**/

    /*update core number and mcd for BLACK nodes*/
    for (const node_t w: Vcolor){
        if (likely(V[w].color == BLACK)) { // in V*
            OrderDelete(w);// remove w from the order list O_k
            core[w]++;
            Vblack.push_back(w);
        }
    }

    // insert Vb2 (ordered black nodes) to the head of O_k+1
    //MultiOrderInsert(H[K+1], Vblack);
    // node_t pos_x = H[K+1];
    // for(const node_t v: Vblack) {
    //     OrderInsert(pos_x, v); pos_x = v;
    // }
    OrderInsertHeadBatch(K+1, Vblack);
    
    /*update the mcd*/
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

    //reset all color
    for(const node_t u: Vcolor) {V[u].color = WHITE; V[u].degin = 0;}

    return 0;
}

void SeqCM::CoreMaint::Forward(node_t w){
    V[w].color = BLACK; cnt_Vs++; cnt_Vp++;// core number can be update 
    Vcolor.push_back(w); // all colored nodes
    //for (const edge_t w2: graph[w]) { // w -> w2
    for (int idx = graph.begin[w]; idx < graph.end[w]; idx++){
        node_t w2 = graph.node_idx[idx];
        
        cnt_S++;
        if (SameCoreOrder(w, w2)) { // outgoing edges
            V[w2].degin++;  // w is black
            if (!V[w2].inQ) { // w2 is not in PQ
                DATA data = DATA(w2, getTag(w2), V[w2].subtag, g_tag_version[core[w2]]);
                if(!PQ.push(data)) 
                    PQ_update(core[w]);
                V[w2].inQ = true; 
            }
        }
        
    }
}


/*update all nodes in PQ with new tag version
* we can optimize not use v2. 
* this is also work for parallel version*/
void SeqCM::CoreMaint::PQ_update(const core_t k){
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

void SeqCM::CoreMaint::Backward(node_t w) { //???
    V[w].color = GRAY; cnt_Vp++; // white -> gray. 
    Vcolor.push_back(w); // all colored nodes. 
    
    node_t p = w;

    DoPre(w); //first time do black pre only;
    V[w].degout += V[w].degin; V[w].degin = 0;

    V[w].color = WHITE;  /*without set E2*/

    while (!R.empty()) { // Q has nodes from BLACK to GRAY 
        node_t u = R.top(); R.pop(); V[u].inR = false; 
        V[u].degout += V[u].degin; V[u].degin = 0;
        V[u].color = GRAY;
        cnt_Vs--; // black -> gray
        

        DoAdj(u);

        V[u].color = WHITE;  /*without set E2*/

        OrderDelete(u); 
        OrderInsert(p, u); p = u;
        ++cnt_ominsert_mid;
    }     
}

void SeqCM::CoreMaint::DoPre(node_t u) {
    //for (const edge_t v: graph[u]) { 
    for (int i = graph.begin[u]; i < graph.end[u]; i++){
        node_t v = graph.node_idx[i];
        if (SameCoreOrder(v, u) && BLACK == V[v].color) { //u<-v  pre
            V[v].degout--;
            //ASSERT(V[v].degout >= 0);

            //if (BLACK == V[v].color  && V[v].degin + V[v].degout <= core[v]) {
            if (V[v].degin + V[v].degout <= core[v]) {
                R.push(v); V[v].inR = true;
            }

        }
    }   
}


void SeqCM::CoreMaint::DoAdj(node_t u) {
    //for (const edge_t v: graph[u]) { 
    for (int idx = graph.begin[u]; idx < graph.end[u]; idx++){
        node_t v = graph.node_idx[idx];

        if (SameCoreOrder(v, u) && BLACK == V[v].color) { //u<-v  pre
            V[v].degout--; 
            ASSERT(V[v].degout >= 0);
            
            //if (!V[v].inR && BLACK == V[v].color && V[v].degin + V[v].degout <= core[v]) {
            if (!V[v].inR && V[v].degin + V[v].degout <= core[v]) {
                /*V[v].color = GRAY; // black -> gray*/
                R.push(v); V[v].inR = true;
            }

        } else if (SameCoreOrder(u, v) && V[v].degin > 0) { // u->v post
            V[v].degin--;

            if (v == 1771) {
                printf("debug 1771 2");
            }

            if (!V[v].inR && BLACK == V[v].color && V[v].degin + V[v].degout <= core[v]) {
                //V[v].color = GRAY; 
                R.push(v); V[v].inR = true;
            }
        }
    }
    
}



/********************************
* edge remoal has different version.
*******************************/
int SeqCM::CoreMaint::EdgeRemove(node_t u, edge_t v) {

    if(g_debug) printf("**************remove (%d, %d)************\n", u, v);

	if(1685 == u && 3322 == v) {
		printf("start debug\n");
	}

    const core_t K = std::min(core[u], core[v]);
    
    const core_t coreu = core[u], corev = core[v];
    if (coreu <= corev) {
        CheckMCD(u, K);
		DoMCD(u, K, true); 
	}
    if (corev <= coreu) {
        CheckMCD(v, K);
		DoMCD(v, K, true);
	}

    graph.remove(u, v);

    while(!R.empty()) {
        node_t w = R.top(); R.pop();
        cnt_Vs++;
        //for (const node_t w2: graph[w]) {
        for (int idx = graph.begin[w]; idx < graph.end[w]; idx++){
            node_t w2 = graph.node_idx[idx];
			cnt_S++; // count all visited edges except calculate mcd.
            if (K == core[w2]) { 
                CheckMCD(w2, K);
                DoMCD(w2, K, true); 
            }
        } 
        V[w].inR = false;
    }

    return 0;
}


void SeqCM::CoreMaint::DoMCD(node_t u, const core_t K, bool off) {
    if(off) --V[u].mcd;
    if (V[u].mcd < K) {
        
        //all three has to be atomic in parallel
        --core[u]; 
        if(!V[u].inR) {
            R.push(u); 
            V[u].inR = true;
            cnt_Vs++;
        }
        
		V[u].mcd = EMPTY;
        // append u to the tail of O_(K-1). 
		OrderDelete(u);
        //OrderInsert(V[T[K-1]].pre, u); 
        OrderInsertTail(K-1, u);
    }
    
}

/*return 1: mcd is empty and then calculate
* return 0: mcd is not empty*/
int SeqCM::CoreMaint::CheckMCD(node_t u, const core_t K) {
    if (EMPTY == V[u].mcd) {
        deg_t mcd = 0;
        for (int idx = graph.begin[u]; idx < graph.end[u]; idx++){
            node_t v = graph.node_idx[idx];
            if(core[v] >= K) {
                ++mcd;
            }
			if (K-1 == core[v] && V[v].inR) {
				++mcd;
			}
        }
        V[u].mcd = mcd;
        return 1;
    }
    return 0;
}

/*by inserting a large number of nodes to triger relabeling.
* remove tail.pre and insert after head.next. */
int SeqCM::CoreMaint::TestOM(size_t insert_size) {
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
    group_t hg = Hg[test_k]; group_t tg = Tg[test_k];

    vector<node_t> nodes; int size = 0;
    for(node_t u = V[h].next; u != t; u = V[u].next) {
        nodes.push_back(u);
        size++; 
        if (size >1000) break;
    }

    for (int i = 0; i < nodes.size(); i++) {
        OrderDelete(nodes[i]);
    }

    for(int i = 0; i < nodes.size(); i++) {
        OrderInsert(h, nodes[i]);
    }

#if 0
    for(size_t i = 0; i < insert_size; i++) {
        
#if 0  // test insert to head
        vector<node_t> y;
        for (int i = 0; i < 17; i++) {
            node_t v = V[t].pre;
            OrderDelete(v);
            y.push_back(v);
        }
        OrderInsertHeadBatch(test_k, y);
#endif

#if 0
        node_t v = V[t].pre;
        OrderDelete(v);
        OrderInsertTail(test_k, v);
        //OrderInsert(h, v);
        // insert at head
#endif
       
        // insert after tail
    }
#endif 
    CheckOM(true);
    
    return 0;
}

/*after test, check the OM*/
int SeqCM::CoreMaint::CheckOM(bool info) {
    //check reverse order list  
    size_t num1 = 0, num2 = 0; 
    for(size_t k = 0; k < H.size(); k++){
        vector<node_t> order, rorder;
        for(node_t v = V[H[k]].next; v != T[k]; v = V[v].next) {
            if (v == -1) {
                assert(v != -1);
            }
            order.push_back(v); num1++;
        }

        for(node_t v = V[T[k]].pre; v != H[k]; v = V[v].pre) {
            assert(v != -1);
            rorder.push_back(v); num2++;
        }
        std::reverse(rorder.begin(), rorder.end());

        ASSERT(order == rorder);
    }
    if (num1 != num2 || num1 != n) {
        if(info) printf("wrong! *** total number in sublist: %ld but correct %ld", num1, n);
        else ASSERT(false);

    }

    num1 = num2 = 0;
    for(size_t k = 0; k < H.size(); k++){
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
    for(size_t k = 0; k < H.size(); k++){
        for(node_t v = V[H[k]].next; v != T[k]; v = V[v].next) {
            node_t next = V[v].next;
            if (V[v].group == V[next].group && V[v].subtag >= V[next].subtag) {
                if(info)  printf("k=%d, %d.subtag: %llu, next %d.subtag: %llu\n",
                    k, v, V[v].subtag, next, V[next].subtag);
                else ASSERT(false);
            }

        }

        for(node_t g = G[Hg[k]].next; g != Tg[k]; g = G[g].next) {
            node_t next = G[g].next;
            if (G[g].tag >= G[next].tag) {
                if(info)  printf("k=%d, %d.tag: %llu, next %d.tag: %llu\n",
                    k, g, G[g].tag, next, G[next].tag);
                else ASSERT(false);
            }

        }

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

int SeqCM::CoreMaint::CheckAllMCDValue(bool info) {
    for (size_t u = 0; u < n; u++){
		if (EMPTY == V[u].mcd) continue; // jump over the empty mcd.

        core_t mcd = 0;
        for (int idx = graph.begin[u]; idx < graph.end[u]; idx++){
            node_t v = graph.node_idx[idx];
            if (core[u] <= core[v]) mcd++;
        }

        if (mcd != V[u].mcd) { //check mcd
            if(info) printf("%d.mcd = %d, should be %d\n", u, V[u].mcd, mcd);
            else ASSERT(false);
        }
    }
}

/*debug: check for result*/
int SeqCM::CoreMaint::CheckCore(vector<core_t> &tmp_core, bool info){
    int r = 0;
    for(node_t i = 0; i < n; i++) {
        if(tmp_core[i] != core[i]){
            r = 1;
            if(info) printf("wrong! *** %d: core is %d, but correct is %d \n", i, core[i], tmp_core[i]);
        }
    }
    return r;
}

void SeqCM::CoreMaint::PrintOMVersion() {
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
int SeqCM::CoreMaint::CheckAll(node_t x, node_t y, int id, vector<core_t> &tmp_core, vector<node_t> &order_v) {

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
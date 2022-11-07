#include <unistd.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <algorithm>

#include "core.h"
#include "defs.h"
#include "gadget/gadget.h"
#include "glist/glist.h"
#include "traversal/traversal.h"
#include "ours-csr-new/seq-csr-new.h"
#include "ours-csr-new/par-csr-new.h"
#include "ours-csr-new/par-om.h"

/*counting*/
int cnt_edge = 0;
int cnt_edge2 = 0;
int cnt_Vs = 0; // count of size V*
int cnt_Vp = 0; // count of size V+
int cnt_S = 0;  // count of all visited
int cnt_tag = 0; 
int cnt_rebalance_max_tag = 0;
int cnt_gap_adjust = 0;

int cnt_omorder = 0;
int cnt_ominsert = 0;
int cnt_omdelete = 0;
int cnt_ominsert_mid = 0; // not append to head or tail. 
int cnt_omsplit = 0;
int cnt_Q_size = 0;

int global_check_num = 0;

int g_debug = 0; // for debuging
unsigned int g_lock_type = CAS_LOCK; // CAS_LOCK; //OMP_LOCK 

int optTest = 0;
int optBufferQ = 0;

ParCM::Count g_cnt;
int g_wait_i = 0; 

int main(int argc, char** argv) {
    char path[128];      // path of input graph
    int method = 0;
    double ratio = 0.0;  // ratio of edges to insert
    int optInsertNum = 0;
    int optRemoveNum = 0;
    int temp = 0;   // temporal graph 
    float edgePercent = 100;

    int optWorkerNum = 0;

    if(argc < 2) {
        printf("*************Test Core Maint*********-t < 20*****\n");
        printf("-p: the path of the data file\n");
        printf("-r: the proportion of edges that is sampled out for insertion/removal\n");
        printf("-I: the number of insert edges that for operation\n");
        printf("-R: the number of remove edges that for operation\n");
        printf("-w: the number of workers (1 - 64)\n");
        printf("-l: the type of lock, 0 CAS_LOCK, 1 OMP_LOCK(default), 2 NO lock\n");
        printf("-m: the algorithm to be used: 1 traversal, 2 order-based, 3 ours-seq, 4 ours-parallel\n");
        //printf("-M: the sub algorithm to be used: 1 ours batch insertion\n");
        printf("-T: 1 for temporal graphs, 0 for ordinary graphs sample edges, 2 for debug without sampled edges\n");
        printf("    3 csr sample edges, 4 csr without sample edges, 5 csr with repeated random \n");
        printf("    6 read temporal graphs ('u v' format sorted by time stamp)\n");
        printf("    7 read temporal graphs with smaple edges\n");
        printf("-c: 1 Output the statistic of core numbers to file *.bin.core for original graph\n");
        printf("    2 Output the statistic of sampled 10k edges core numbers to \n");
        printf("-s: [<=100]the percent samele the percent of edges \n");
        printf("    [>=100000] the number of sameled of edges, 100k, 200k ... \n");
        printf("-d: set debug.\n");
        printf("-t: 0 default. 1. remove without init mcd\n"); 
        printf("-b: buffer queue for order maint. 0 No, 1 Yes\n");
        printf("for example: core -p path -I 100000 -m 4 -T -w 16\n");
        printf("\n");

        printf("*************Test Order Maint (OM data structure)**********with -t >=20***\n");
        printf("-t: 20 (no relable) 21 test OM REPEAT_RANDOM (few), 22 FIXED_MULTIPLE (many), 23 FIXED_ONE (max)\n"); 
        printf("-p: the path of the random number file\n");
        printf("-I: the number nodes for operation\n");
        printf("-w: the number of workers (1 - 64)\n");
        printf("-l: the type of lock, 0 CAS_LOCK, 1 OMP_LOCK(default), 2 NO lock\n");
        printf("-d: set debug.\n");
        printf("-T: 1 repeated random positions, 2 on-the-fly random position.\n");
        printf("-s: 1 sorted postion, 0 unsorted position (default).\n ");
        printf("for example: ./core2 -p path -I 10000000 -T 1 -t 20 -w 16\n");

        printf("\n");
        printf("*************Trasfer graph to edge format**************** with -t = 30****\n");
        printf("-t: 30 to edge file: *.edge, 31 to inserted edge file: *.edge-10... \n");
        printf("-p: the path of the random number file\n");
        printf("-T: same as before\n");
        printf("-I: the inserted or removed number of edges: *-edge.txt\n");

        exit(0);

    }
    // initialize the options
    int option = -1;
    int outputCore = 0;
    while (-1 != (option = getopt(argc, argv, "p:r:I:R:w:l:m:M:T:c:s:d:t:b:"))) {
        switch (option) {
        case 'p':
            strcpy(path, optarg);
            break;
        case 'r':
            ASSERT(0.0 <= (ratio = atof(optarg)) && 1.0 >= ratio);
            break;
        case 'I':
            optInsertNum = atoi(optarg);
            break;
        case 'R':
            optRemoveNum = atoi(optarg);
            break;
        case 'w':
            optWorkerNum = atoi(optarg);
            break;
        case 'l':
            g_lock_type = atoi(optarg);
            break;
        case 'm':
            method = atoi(optarg);
            break;
        case 'T':
            temp = atoi(optarg);
            break;
        case 'c':
            outputCore = atoi(optarg);
            break;
        case 's':
            edgePercent = atof(optarg);
            break;
        case 'd':
            g_debug = atoi(optarg);
            break; 
        case 't':
            optTest = atoi(optarg);
            break;
        case 'b':
            optBufferQ = atoi(optarg);
            break;
        }
    }

    bool ispar = true; 
    if (optWorkerNum < 1) ispar = false;

    /********************test parallel OM******************/
    if(optTest >= 20 && optTest < 30) {
        
        printf("********* Test Parallel OM !\n");
        printf("init size %d\n", optInsertNum); // ParOM::OM_SIZE
        printf("insert size %d\n", optInsertNum);
        printf("Lock Type: %d\n", g_lock_type);
        printf("\t 0 ***CAS_LOCK, 1 OMP_LOCK(default), 2 NO lock\n");

        vector<int> pos;


        {
            /*the insertnum = initial size**/
        //ParOM::OM *parom = new ParOM::OM( (optInsertNum * 2) + 10000);
        ParOM::OM *parom = new ParOM::OM( (optInsertNum * 2) + 10000);
        //vector<node_t> nodes = parom->AlocateNodes(optInsertNum); // ParOM::OM_SIZE
        vector<node_t> nodes = parom->AlocateNodes(optInsertNum); // ParOM::OM_SIZE
        parom->Init(nodes); // insert 1m initially.
        parom->InitParallel(optWorkerNum);
        nodes = parom->AlocateNodes(optInsertNum);

        

        if (temp == 1) {
            printf("*** with repeated random position!\n");
            vector<int> rand = parom->GetRepeatRandomNum(path);

            if (20 == optTest) {
                pos = parom->GenerateTestCase(ParOM::NO_RELABEL, rand, optInsertNum);
                printf("position size %d\n", ParOM::OM_POS_NO);
            } else if (21 == optTest) { //case 1 REPEAT_RANDOM
                pos = parom->GenerateTestCase(ParOM::FEW_RELABEL, rand, optInsertNum);
                printf("position size %d\n", ParOM::OM_POS_FEW);
            } else if (22 == optTest) { // case 2 FIXED_MULTIPLE
                pos = parom->GenerateTestCase(ParOM::MANY_RELABEL, rand, optInsertNum);
                printf("position size %d\n", ParOM::OM_POS_MANY);
            } else if (23 == optTest) { // case 3 FIXED_ONE 
                pos = parom->GenerateTestCase(ParOM::MAX_RELABEL, rand, optInsertNum);
                printf("position size %d\n", ParOM::OM_POS_MAX);
            } 
        } else if (temp == 2){   
 
            printf("*** with random position!\n");
            if (20 == optTest) {
                // insert 10m in to 10m positions. 
                pos = parom->GenerateTestCase2(ParOM::NO_RELABEL,  optInsertNum);
                printf("position size %d\n", ParOM::OM_POS_NO);
            } else if (21 == optTest) { //case 1 REPEAT_RANDOM
                // insert 10m into 1m positions. 
                pos = parom->GenerateTestCase2(ParOM::FEW_RELABEL,  optInsertNum);
                printf("position size %d\n", ParOM::OM_POS_FEW);
            } else if (22 == optTest) { // case 2 FIXED_MULTIPLE
                //insert 10m into 1000 positions.
                pos = parom->GenerateTestCase2(ParOM::MANY_RELABEL,  optInsertNum);
                printf("position size %d\n", ParOM::OM_POS_MANY);
            } else if (23 == optTest) { // case 3 FIXED_ONE 
                //insert 10m into 1 position.
                pos = parom->GenerateTestCase2(ParOM::MAX_RELABEL,  optInsertNum);
                printf("position size %d\n", ParOM::OM_POS_MAX);
            } 
        }

        //for debug
        // if (edgePercent = 1) {
        //     std::sort(pos.begin(), pos.end());
        // }

        auto test_beg = std::chrono::steady_clock::now();
        parom->TestInsert(pos, nodes, ispar);
        auto test_end = std::chrono::steady_clock::now();
        auto test_dif = test_end - test_beg;
        printf("TestInsert costs(ms): %f\n", std::chrono::duration<double, std::milli>(test_dif).count());
        printf("rebalance max tag #: %d\n\n", cnt_rebalance_max_tag);
        
        // test OM order 
        test_beg = std::chrono::steady_clock::now();
        parom->TestOrder(pos, ispar);
        test_end = std::chrono::steady_clock::now();
        test_dif = test_end - test_beg;
        printf("TestOrder costs(ms): %f\n\n", std::chrono::duration<double, std::milli>(test_dif).count());
        

        // test OM delete 
        /*delete the nodes that added firstly. 
        in this case, the delete can always be in parallel.*/
        test_beg = std::chrono::steady_clock::now();
        parom->TestDelete(nodes, ispar);
        test_end = std::chrono::steady_clock::now();
        test_dif = test_end - test_beg;
        printf("TestDelete costs(ms): %f\n\n", std::chrono::duration<double, std::milli>(test_dif).count());
        
        parom->PrintCount("");
        // test OM mixed, this has bugs. 
        //delete parom;
        }

#if 1  // test the mixed operations
        {
        ParOM::OM *parom = new ParOM::OM( (optInsertNum * 2) + 10000 ); // ParOM::OM_SIZE
        vector<node_t> nodes = parom->AlocateNodes(optInsertNum); //ParOM::OM_SIZE
        parom->Init(nodes); // insert 1m initially.
        parom->InitParallel(optWorkerNum);
        nodes = parom->AlocateNodes(optInsertNum);

        // if (21 == optTest) { //case 1 REPEAT_RANDOM
        //     pos = parom->GenerateTestCase(ParOM::REPEAT_RANDOM, rand, 3 * optInsertNum);
        //     printf("position size %d\n", optInsertNum * 3);
        // } else if (22 == optTest) { // case 2 FIXED_MULTIPLE
        //     pos = parom->GenerateTestCase(ParOM::FIXED_MULTIPLE, rand, 3 *  optInsertNum);
        //     printf("position size %d\n", ParOM::OM_POS_SIZE2);
        // } else if (23 == optTest) { // case 3 FIXED_ONE 
        //     pos = parom->GenerateTestCase(ParOM::FIXED_ONE, rand, 3 * optInsertNum);
        //     //printf("position size %d\n", ParOM::OM_POS_SIZE3);
        // }


        cnt_tag = 0;
        cnt_rebalance_max_tag = 0;
        auto test_beg = std::chrono::steady_clock::now();
        parom->TestMixed(pos, nodes, ispar);
        auto test_end = std::chrono::steady_clock::now();
        auto test_dif = test_end - test_beg;
        printf("TestMixed costs(ms): %f\n", std::chrono::duration<double, std::milli>(test_dif).count());
        printf("mixed rebalance max tag #: %d\n\n", cnt_rebalance_max_tag);
        
        parom->PrintCount("Mixed ");
        
#endif
        //delete(parom);
        exit(0);
        }
    }
    /********************test parallel OM END!******************/



    // read the graph
    int n, m, m2;
    std::vector<std::pair<int, int>> edges;
    if(1 == temp) {
        //const auto read_func = temp ? gadget::ReadTempEdgesS : gadget::ReadEdgesS;
        edges = gadget::ReadTempEdgesS(path, &n, &m);
    } else if (0 == temp) {
        edges = gadget::ReadEdgesS(path, &n, &m); // random shuffle the edges. 
        
    }  else if (2 == temp) {
        edges = gadget::ReadEdges(path, &n, &m); // without shuffle
        
    } else if (3 == temp ) {
        edges = gadget::CSRReadEdges(path, &n, &m, 1); // csr with shuffle
    } else if (4 == temp ) {
        edges = gadget::CSRReadEdges(path, &n, &m, 0); // csr without shuffle.
    } else if (5 == temp) {
        edges = gadget::CSRReadEdges(path, &n, &m, 2); // csr with repeated random shuffle.
    } else if (6 == temp) {
        edges = gadget::ReadTemporalEdges(path, &n, &m, 0); // read temporal graph with time ordered edges     
    } else if (7 == temp) {
        edges = gadget::ReadTemporalEdges(path, &n, &m, 1); // read temporal graph with shuffle.
    } else {
        printf("wrong graphs\n");
        exit(1);
    }

    //sample the edges by percentage
    if (edgePercent < 100) {
        edges = gadget::sampleEdges(path, edges, edgePercent);
        printf("edges size: %f\%\n", edgePercent);
    } else if (edgePercent >= 100000){
        if (edgePercent < m) {
            // edgepercent is the number of sampled edges.
            gadget::CutEdges(edges, edgePercent);
            printf("edges size: %d\n", edges.size());
        } else {
            // all edges
            edgePercent = m;
            printf("edges size: %d\n", edges.size());
        }
    } 
    m = edges.size();

  
  // get the number of edges for inserting or deleting.
  if (optInsertNum > 0) { m2 = m - optInsertNum;} 
  if (optRemoveNum > 0) { m2 = m;} //load the whole graph for removing.
  if (ratio > 0.0){ m2 = m - static_cast<int>(m * ratio);}


#if 1  // output the to graph as edge files for test other algorithms. 
    if (optTest == 30) {

        std::string newpath(path);
        newpath += ".edge";
        auto file = fopen(newpath.c_str(), "w");
        for(auto e: edges) {
            std::string line;
            line += std::to_string(e.first);
            line += "\t";
            line += std::to_string(e.second);
            line += "\n";
            fputs(line.c_str(), file);
        }
        fclose(file);
        printf("store graph as %s\n", newpath.c_str());
    }

    if (optTest == 31) {
        std::string newpath(path);
        newpath += ".edge-";
        newpath += std::to_string(optInsertNum);
        auto file = fopen(newpath.c_str(), "w");
        for (int i = m - 1; i >= m2; --i) {
            auto e = edges[i];
            std::string line;
            line += std::to_string(e.first);
            line += "\t";
            line += std::to_string(e.second);
            line += "\n";
            fputs(line.c_str(), file);
        }
        fclose(file);
        printf("store inserted or removed edges as %s\n", newpath.c_str());
        exit(0);
    }

   
#endif 


    // print the configurations
    gadget::RepeatWith('*', 80);
    printf("# of vertices: %d\n", n);
    printf("# of (all) edges: %d\n", m);
    printf("ratio: %f\n", ratio);
    printf("# of edges to insert: %d\n", optInsertNum);
    printf("# of edges to delete: %d\n", optRemoveNum);
    printf("path of the input file: %s\n", path);
    printf("method: %d  (2 glist, 1 traversal, or 3 ours, 4 ours parallel)\n", method);
    printf("graph: %d  (1 for temporal graphs, 0 for ordinary graphs sample edges, 2 for debug without sampled edges)\n", temp);
    printf("            3 csr sample edges, 4 csr without sample edges, 5 csr with repeated random\n");
    printf("# worker: %d \n", optWorkerNum);
    gadget::RepeatWith('*', 80);
    
   
    


    printf("creat adjacent list\n");

    // compute the base core
    std::vector<core_t> core(n); 
    std::vector<core_t> parcore(n); // used in parallel version. 
    std::vector<core_t> order_v(n);
    std::vector<core_t> tmp_core(n);

    GRAPH graph = gadget::CreateOurCSRGRAPH(edges, n, m2);

    //in case graph is deleted.
    //PARGRAPH pargraph = gadget::CreateOurCSRParGRAPH(edges, n, m2);;

    // initialize the core component
    core::CoreMaintenance* cm = nullptr;
    SeqCM::CoreMaint *ourcm = nullptr;
    ParCM::CoreMaint *ourparcm = nullptr;
    if (method%10 == 1) { //traversal
        cm = new core::Traversal(n);
    } else if (method%10 == 2) { //glist

        cm = new core::GLIST(n);
    } else if (method%10 == 3 ) {
        ourcm = new SeqCM::CoreMaint(n, graph, core);
    } else if (method == 4) {
        ourparcm = new ParCM::CoreMaint(n, graph, parcore);
    } else {
        printf("wrong method! %d\n", method);
        exit(0);
    }

    printf("new the working object\n");

   

#if 0  // output the core numbers. 
     //statistic the core numbers
    if (1 == outputCore) {
        for (auto e: edges) {
            graph[e.first].push_back(e.second);
            graph[e.second].push_back(e.first);
        }
        ourcm = new SeqCM::CoreMaint(n, graph, core);
        ourcm->ComputeCore(graph, core, order_v, false);
        gadget::OutputCoreNumber(path, core, n);
        goto END;
    } else if (2 == outputCore) {
        for (auto e: edges) {
            graph[e.first].push_back(e.second);
            graph[e.second].push_back(e.first);
        }
        ourcm = new SeqCM::CoreMaint(n, graph, core);
        ourcm->ComputeCore(graph, core, order_v, false);
        gadget::OutputSampleEdgeCoreNumber(path, core, edges, n);
        goto END;
    }
#endif // output the core numbers. 


    /*initialization*/
    printf("********* initialization!\n");
    core_t max_core; 
    const auto init_beg = std::chrono::steady_clock::now();
    if (method%10 != 3) { //3 is ours method 
        //cm->ComputeCore(graph, true, core);
    }
    if (method == 3){ // ours 
        ourcm->ComputeCore(core, order_v, true);       // compute the k-order by bz. 
        ourcm->Init(order_v);   // init to our class.
    } else if (method == 4){
        max_core = ourparcm->ComputeCore(parcore, order_v, true);       // compute the k-order by bz. 
        ourparcm->Init(order_v);   // init to our class.
        ourparcm->InitParallel(optWorkerNum);

        printf("Max Core number: %d\n", ourparcm->getMaxCore());
    }
    const auto init_end = std::chrono::steady_clock::now();
    const auto init_dif = init_end - init_beg;
    printf("initialization costs(ms): %f\n", std::chrono::duration<double, std::milli>(init_dif).count());  
    
    printf("max core number: %d\n", max_core);

    if (optRemoveNum > 0) goto REMOVE;



#if 1    // insert edge 
    printf("*********begin insertion %d %d!\n", m2, m);
    const auto ins_beg = std::chrono::steady_clock::now();
    
    if (3 == method ) { // new sequential
        for (int i = m2; i < m; ++i) {
            ourcm->EdgeInsert(edges[i].first, edges[i].second);
            if (g_debug > 0) {
                // check each time. 
                ourcm->ComputeCore(tmp_core, order_v, false);
                //int r = ourcm->Check(edges[i].first, edges[i].second, i - m2, tmp_core, order_v);  
            }
        }
    } else if (4 == method) { // new parallel
        if (ispar) {
            ourparcm->ParallelEdgeInsert(edges, m, m2, ispar);
        } else {
            const size_t reserve_size = n;
            ParCM::PRIORITY_Q PQ(reserve_size);
            ParCM::QUEUE R(reserve_size);
            vector<node_t> Vblack; Vblack.reserve(reserve_size);
            vector<node_t> Vcolor; Vcolor.reserve(reserve_size);
            vector<int> Qdegin(reserve_size, 0);
            vector<bool> Qin(reserve_size, false);

            vector<node_t> groups; groups.reserve(reserve_size);

            for (int i = m - 1; i >= m2; --i) {
                node_t u = edges[i].first;
                node_t v = edges[i].second;
                ourparcm->EdgeInsert(u, v, PQ, R, Vblack, Vcolor, Qdegin, Qin, groups, 0);
                if(3 == g_debug) {
                    printf("************** sequential insert (%d, %d)************\n", u, v);
                    // ourparcm->ComputeCore(tmp_core, order_v, false);
                    // int result = ourparcm->CheckCore(tmp_core, true);
                    // if (0  == result) ERROR("check core number passed", false);
                    // else ERROR("check core number fail!!!", false);
                    // ourparcm->CheckLock(true);
                    //ourparcm->CheckOM(true);
                    ourparcm->CheckDeg(true);
                }
                
            }
        }
   
    }
    
    


  const auto ins_end = std::chrono::steady_clock::now();
  const auto ins_dif = ins_end - ins_beg;
  //printf("core insert costs \x1b[1;31m%f\x1b[0m ms\n",
         //std::chrono::duration<double, std::milli>(ins_dif).count());
    printf("\ncore IorR costs(ms): %f\n\n", std::chrono::duration<double, std::milli>(ins_dif).count());

    
    printf("V* size: %d\n", cnt_Vs);
    printf("V+ size: %d\n", cnt_Vp);
    printf("S size: %d\n", cnt_S);
    printf("our adjust tag: %d\n", cnt_tag);
    printf("our max relabel tag: %d\n", cnt_rebalance_max_tag);
    
  /*see count end*/


goto END;
#endif // insert edge. 

REMOVE:
#if 1   // remove edges
  printf("begin remove !\n");
  if (optRemoveNum > 0) { m2 = m - optRemoveNum;} 

    const auto rmv_beg = std::chrono::steady_clock::now();
  
    if (method == 3) {
        for (int i = m - 1; i >= m2; --i) {
            ourcm->EdgeRemove(edges[i].first, edges[i].second);
            if (g_debug) {
                ourcm->ComputeCore(tmp_core, order_v, false);
                int result = ourcm->CheckCore(tmp_core, true);
                if (0  == result) ERROR("check core number passed", false);
                else ERROR("check core number fail!!!", false);
                ourcm->CheckAllMCDValue(true);
                
            }
        }
    }

    if (method == 4) {
        
        if (ispar) {
            ourparcm->ParallelEdgeRemove(edges, m, m2, ispar);
            if (g_debug) {
                ourparcm->ComputeCore(tmp_core, order_v, false);
                int result = ourparcm->CheckCore(tmp_core, false);
                if (0  == result) ERROR("check core number passed", false);
                else ERROR("check core number fail!!!", false);
                ourparcm->CheckAllMCDValue(true);
                ourparcm->CheckLock(true);
            }
        } else { //not parallel
            
            // const auto mem_beg = std::chrono::steady_clock::now();

            const size_t reserve_size = n;
            ParCM::QUEUE R(n);
            vector<node_t> Vblack; Vblack.reserve(reserve_size);
            vector<node_t> groups; groups.reserve(reserve_size);
            vector<bool> A(reserve_size, false);
            
            for (int i = m-1; i >= m2; --i) {
                ourparcm->EdgeRemove(edges[i].first, edges[i].second, R, Vblack, A, groups, 0);
                if (g_debug > 1) {
                    ourparcm->ComputeCore(tmp_core, order_v, false);
                    int result = ourparcm->CheckCore(tmp_core, true);
                    if (0  == result) ERROR("check core number passed", false);
                    else ERROR("check core number fail!!!", false);
                    ourparcm->CheckAllMCDValue(true);
                    ourparcm->CheckLock(true);
                }
            }
        }
      
    }
  const auto rmv_end = std::chrono::steady_clock::now();
  const auto rmv_dif = rmv_end - rmv_beg;

  printf("\ncore IorR costs(ms): %f\n\n", std::chrono::duration<double, std::milli>(rmv_dif).count());  

  /*see count*/
    printf("temp edge number 1: %d\n", cnt_edge);
    printf("temp edge number 2: %d\n", cnt_edge2);
    
    printf("V* size: %d\n", cnt_Vs);
    printf("V+ size: %d\n", cnt_Vp);
    printf("total visited edge S size: %d\n", cnt_S);
    printf("our adjust tag: %d\n", cnt_tag);
    printf("our max relabel tag: %d\n", cnt_rebalance_max_tag);
  /*see count end*/

#endif // end remove

END:
    // verify the result for both insertion and remove
    {
        int result = 0;
        if (method%10 != 3) {

        } else {
            ourcm->ComputeCore(tmp_core, order_v, false);
            result = ourcm->CheckCore(tmp_core, false);
            if (0  == result) {
                ERROR("check core number passed", false); 
                printf("\t by method %d\n", method);
            }
            //ourcm->Check(-1, -1, -1, tmp_core, order_v);
        }

        if (method == 4) {
            printf("wait i: %d\n", g_wait_i);
            printf("\n");
            for (int i = 0; i < ourparcm->cnt.size(); i++) {
                g_cnt = g_cnt + ourparcm->cnt[i];
            }

            printf("# of relabel: %d\n", g_cnt.relabel);
            printf("# of tag: %d\n", g_cnt.tag);
            printf("# of subtag: %d\n", g_cnt.subtag);




            std::string newpath(path);
            if (optRemoveNum > 0) {
                newpath += ".R-Vcolor-";
                newpath += std::to_string(optRemoveNum);
            } else if (optInsertNum > 0) {
                newpath += ".I-Vcolor-";
                newpath += std::to_string(optInsertNum);
            }

            auto file = fopen(newpath.c_str(), "w");

            for(int i = 0; i < g_cnt.Vcount.size(); i++) {
                if (g_cnt.Vcount[i] > 0) {
                    //printf("Vcolor, %d, %d\n", i, g_cnt.Vcount[i]);

                    std::string line;
                    line += std::to_string(i);
                    line += ",";
                    line += std::to_string(g_cnt.Vcount[i]);
                    line += "\n";
                    fputs(line.c_str(), file);
                }
            }
            
            fclose(file);
            printf("\n");

          


            ourparcm->ComputeCore(tmp_core, order_v, false);
            printf("Max Core number: %d\n", ourparcm->getMaxCore());

            result = ourparcm->CheckCore(tmp_core, true);
            if (0  == result) {
                ERROR("check core number passed", false); 
                printf("\t by method %d", method);
            } else {
                ERROR("check core number fail", false); 
                printf("\t by method %d", method);
            }

            printf(" with %d workers \n", optWorkerNum);



            //check edge insert
            if (optInsertNum > 0) {
                ourparcm->CheckOM(true);
                ourparcm->CheckDeg(true);
                ourparcm->CheckLock(true);
            }

            //check edge remove    
            if (optRemoveNum > 0) {
                ourparcm->CheckOM(true);
                ourparcm->CheckLock(true);
                ourparcm->CheckAllMCDValue(true);
            }

            
        }
    }



ENDEXIT:
    if (cm != nullptr) delete cm;
    if (ourcm != nullptr) delete ourcm;
    if (ourparcm != nullptr) delete ourparcm;

}

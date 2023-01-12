#include <unistd.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>

#include "core.h"
#include "defs.h"
#include "gadget/gadget.h"
#include "glist/glist.h"
#include "traversal/traversal.h"


int main(int argc, char** argv) {
  char path[100];      // path of input graph
  char method[100];    // method
  double ratio = 0.0;  // ratio of edges to insert
  bool temp = false;   // temporal graph
  // initialize the options
  int option = -1;
  while (-1 != (option = getopt(argc, argv, "p:r:m:T:"))) {
    switch (option) {
      case 'p':
        strcpy(path, optarg);
        break;
      case 'r':
        //ASSERT(0.0 <= (ratio = atof(optarg)) && 1.0 >= ratio);
        //change to number
        ratio = atoi(optarg);
        break;
      case 'm':
        strcpy(method, optarg);
        break;
      case 'T':
        temp = atoi(optarg) != 0;
        break;
    }
  }
  // read the graph
  int n, m, m2;
  const auto read_func = temp ? gadget::ReadTempEdgesS : gadget::ReadEdgesS;
  const auto edges = read_func(path, &n, &m);
  
  //m2 = m - static_cast<int>(m * ratio);
    m2 = m - ratio; 

  // print the configurations
  gadget::RepeatWith('*', 80);
  printf("# of vertices: %d\n", n);
  printf("# of (all) edges: %d\n", m);
  printf("ratio: %f\n", ratio);
  printf("# of edges to insert: %d\n", m - m2);
  printf("path of the input file: %s\n", path);
  printf("method: %s\n", method);
  gadget::RepeatWith('*', 80);
  // initialize the core component
  core::CoreMaintenance* cm = nullptr;
  if (strcmp(method, "traversal") == 0) {
    cm = new core::Traversal(n);
  } else if (strcmp(method, "glist") == 0) {
    cm = new core::GLIST(n);
  } else {
    ASSERT(false);
  }
  // create the adjacent list representation
  std::vector<std::vector<int>> graph(n);
  for (int i = 0; i < m2; ++i) {
    int v1 = edges[i].first;
    int v2 = edges[i].second;
    graph[v1].push_back(v2);
    graph[v2].push_back(v1);
  }
  // compute the base core
  std::vector<int> core(n);
  const auto init_beg = std::chrono::steady_clock::now();
  cm->ComputeCore(graph, true, core);
  const auto init_end = std::chrono::steady_clock::now();
  const auto init_dif = init_end - init_beg;
  printf("initialization costs \x1b[1;31m%f\x1b[0m ms\n",
         std::chrono::duration<double, std::milli>(init_dif).count());
  // verify the result
  {
    cm->Check(graph, core);
    // ASSERT(false);
  }
  // insert edges
  const auto ins_beg = std::chrono::steady_clock::now();
  for (int i = m2; i < m; ++i) {
    cm->Insert(edges[i].first, edges[i].second, graph, core);
  }
  const auto ins_end = std::chrono::steady_clock::now();
  const auto ins_dif = ins_end - ins_beg;
  printf("Insert: %f\n",
        std::chrono::duration<double, std::milli>(ins_dif).count());

  // verify the result
  {
    ERROR("check: insert", false);
    std::vector<int> tmp_core(n);
    cm->ComputeCore(graph, false, tmp_core);
    ASSERT_INFO(tmp_core == core, "wrong result after insert");
    cm->Check(graph, core);
    ERROR("check passed", false);
  }
  // remove edges
  const auto rmv_beg = std::chrono::steady_clock::now();
  for (int i = m - 1; i >= m2; --i) {
    cm->Remove(edges[i].first, edges[i].second, graph, core);
  }
  const auto rmv_end = std::chrono::steady_clock::now();
  const auto rmv_dif = rmv_end - rmv_beg;
  printf("Delete: %f\n",
         std::chrono::duration<double, std::milli>(rmv_dif).count());
  // verify the result
  {
    ERROR("check: remove", false);
    std::vector<int> tmp_core(n);
    cm->ComputeCore(graph, false, tmp_core);
    ASSERT_INFO(tmp_core == core, "wrong result after remove");
    cm->Check(graph, core);
    ERROR("check passed", false);
  }
  delete cm;
}

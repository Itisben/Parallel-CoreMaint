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

int main(int, char** argv) {
  const int n = atoi(argv[1]);
  std::vector<std::pair<int,int>> edges;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      edges.push_back(std::make_pair(i, j));
    }
  }
  // randomly shuffle the edges
  srand(time(nullptr));
  for (size_t i = 0; i < edges.size(); ++i) {
    const size_t j = i + rand() % (edges.size() - i);
    std::swap(edges[i], edges[j]);
  }
  char method[100];
  strcpy(method, argv[2]);
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
  std::vector<int> core(n);
  cm->ComputeCore(graph, true, core);
  // insert edges
  const auto ins_beg = std::chrono::steady_clock::now();
  for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
    cm->Insert(edges[i].first, edges[i].second, graph, core);
  }
  const auto ins_end = std::chrono::steady_clock::now();
  const auto ins_dif = ins_end - ins_beg;
  printf("core insert costs \x1b[1;31m%f\x1b[0m ms\n",
         std::chrono::duration<double, std::milli>(ins_dif).count());
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
  for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
    cm->Remove(edges[i].first, edges[i].second, graph, core);
  }
  const auto rmv_end = std::chrono::steady_clock::now();
  const auto rmv_dif = rmv_end - rmv_beg;
  printf("core remove costs \x1b[1;31m%f\x1b[0m ms\n",
         std::chrono::duration<double, std::milli>(rmv_dif).count());
  {
    ERROR("check: insert", false);
    std::vector<int> tmp_core(n);
    cm->ComputeCore(graph, false, tmp_core);
    ASSERT_INFO(tmp_core == core, "wrong result after insert");
    cm->Check(graph, core);
    ERROR("check passed", false);
  }
  delete cm;
}

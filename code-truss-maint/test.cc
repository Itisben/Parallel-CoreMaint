#include <unistd.h>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "defs.h"
#include "graph.h"
#include "order.h"

// a sample program
int main(int argc, char** argv) {
  ASSERT(5 == argc);
  const std::string op = argv[1];
  const std::string old_index_file = argv[2];
  const std::string update_file = argv[3];
  const std::string ground_truth_file = argv[4];
  printf("*****************************************************************\n");
  printf("old index file: %s\n", old_index_file.c_str());
  printf("update file: %s\n", update_file.c_str());
  printf("ground truth file: %s\n", ground_truth_file.c_str());
  printf("*****************************************************************\n");
  // read the header
  std::ifstream infile(old_index_file, std::ios::binary);
  uint32_t n = -1, m = -1;  // # of vertices and edges
  infile.read(reinterpret_cast<char*>(&n), sizeof n)
        .read(reinterpret_cast<char*>(&m), sizeof m);
  infile.close();
  // read the graph and the index
  // We set the param "l" to "m * 2" here, where "l" is the maximum number of
  // edges that a graph can hold. We impose this constraint only for ease of
  // implementation. TODO: remove this constraint in the future.
  truss_maint::Order tm(n, m * 2, old_index_file);
  // read the edges to update
  std::vector<truss_maint::EdgT> inc_edges;
  {
    std::ifstream inc_file(update_file);
    uint32_t inc_m = 0; inc_file >> inc_m;
    for (uint32_t e = 0; e < inc_m; ++e) {
      uint32_t v1, v2;
      inc_file >> v1 >> v2;
      inc_edges.push_back({v1, v2});
    }
    inc_file.close();
  }
  printf("ratio: %f\n", inc_edges.size() / static_cast<double>(m));
  // apply the updates
  const auto beg = std::chrono::steady_clock::now();
  if (op == "uinsert") {
    printf("unit insert used.\n");
    for (const auto edge : inc_edges) {
      tm.Insert({edge});
    }
  } else if (op == "binsert") {
    printf("batch insert used.\n");
    tm.BatchInsert(inc_edges);
  } else if (op == "udelete") {
    printf("unit delete used.\n");
    for (const auto edge : inc_edges) {
      tm.Remove(edge.first, edge.second);
    }
  } else {
    printf("batch delete used.\n");
    tm.BatchRemove(inc_edges);
  }
  const auto end = std::chrono::steady_clock::now();
  const auto dif = end - beg;
  printf("Applying the updates costs \x1b[1;31m%f\x1b[0m ms.\n",
         std::chrono::duration<double, std::milli>(dif).count());
  // verify the result
  printf("Verifying the results: ");
  tm.Debug();
  tm.Check(ground_truth_file);
  printf("Done.\n");
  printf("*****************************************************************\n");
}

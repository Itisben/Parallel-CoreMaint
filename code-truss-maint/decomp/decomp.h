#ifndef TRUSS_MAINT_DECOMP_DECOMP_H_
#define TRUSS_MAINT_DECOMP_DECOMP_H_

#include <cstdint>
#include <string>
#include <vector>

namespace truss_maint {
namespace decomp {
// class Decomp is to truss-decompose a graph; as a byproduct,
// it produces a truss-decomposition order
class Decomp final {
 public:
  // file_name specifies a file that contains a graph of the following form:
  //     n  m
  //     u_1 v_1
  //     u_2 v_2
  //     ...
  //     u_m v_m
  // that is, the first line shows the # of vertices (n) and the # of edges (m);
  // the following m lines show the m edges, where (u_i, v_i) is the i-th edge;
  // vertex IDs are in the range [0, n - 1]; the graph is required to be simple.
  Decomp(const std::string& file_name);
  Decomp(const Decomp&) = delete;
  Decomp& operator=(const Decomp&) = delete;
  // write the results to disk
  void WriteToFile(const std::string& file_name) const;

 private:
  // adjacency array entry type
  typedef struct final {
    std::uint32_t vid;
    std::uint32_t eid;
  } ArrayEntry;
  // data members
  std::uint32_t n_;  // the # of vertices
  std::uint32_t m_;  // the # of edges
  // the adjacency array representation
  std::vector<std::vector<ArrayEntry>> adj_;
  // the truss numbers
  std::vector<std::uint32_t> k_;
  // the remaining supports
  std::vector<std::uint32_t> rem_;
  // the triangle supports
  std::vector<std::uint32_t> ts_;
  // the edge peeling order
  std::vector<std::uint32_t> ord_;
  // the set of edges
  std::vector<std::pair<std::uint32_t, std::uint32_t>> edges_;
};
}  // namespace decomp
}  // namespace truss_maint

#endif

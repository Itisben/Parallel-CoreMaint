#ifndef TRUSS_MAINT_GRAPH_H_
#define TRUSS_MAINT_GRAPH_H_

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>

#include "defs.h"

namespace truss_maint {
using std::int32_t;
using std::uint32_t;
// edge type
typedef std::pair<uint32_t, uint32_t> EdgT;
// graph class
class Graph final {
 public:
  // construct a graph with only n isolated vertices
  Graph(const uint32_t n, const uint32_t l, const std::vector<int32_t>& k)
      : l_(l), n_(n), m_(0), k_(k) {
    ASSERT_MSG(0 < n_ && n_ < (static_cast<uint32_t>(1) << 29),
               "invalid argument");
    ASSERT_MSG(0 < l_ && l_ < (static_cast<uint32_t>(1) << 29),
               "invalid argument");
    // available edges
    free_edges_.resize(l_);
    std::iota(free_edges_.rbegin(), free_edges_.rend(), 0);
    free_.resize(l_, true);
    // adjacency array
    adj_.resize(n_);
    // edge information
    edge_info_.resize(l_, {UINT32_MAX, UINT32_MAX});
  }
  ~Graph() {}
  // get the endpoints of the edge with ID eid
  inline EdgT Get(const uint32_t eid) const {
    ASSERT_MSG(UINT32_MAX != edge_info_.at(eid).first, "invalid edge ID");
    return edge_info_[eid];
  }
  // get the ID of the edge with endpoints v1 and v2
  uint32_t Get(uint32_t v1, uint32_t v2) const {
    if (adj_.at(v1).size() > adj_.at(v2).size()) std::swap(v1, v2);
    for (const auto ae : adj_[v1]) {
      if (ae.vid == v2) return ae.eid;
    }
    ASSERT(false);
    return UINT32_MAX;
  }
  // get the triangles containing the edge with ID eid
  std::vector<std::pair<uint32_t, uint32_t>>
  GetTriangles(const uint32_t eid) const {
    ASSERT_MSG(UINT32_MAX != edge_info_.at(eid).first, "invalid edge ID");
    const uint32_t v1 = edge_info_[eid].first;
    const uint32_t v2 = edge_info_[eid].second;
    ASSERT(v1 < n_ && v2 < n_);
    std::vector<std::pair<uint32_t, uint32_t>> triangles;
    // find common neighbors
    size_t p1 = 0, p2 = 0;
    while (p1 < adj_[v1].size() && p2 < adj_[v2].size()) {
      if (adj_[v1][p1].vid == adj_[v2][p2].vid) {
        triangles.push_back({adj_[v1][p1].eid, adj_[v2][p2].eid});
        ++p1; ++p2;
      } else if (adj_[v1][p1].vid < adj_[v2][p2].vid) {
        ++p1;
      } else {
        ++p2;
      }
    }
    return triangles;
  }
  // get the triangles which contain the edge with ID eid
  std::vector<std::pair<uint32_t, uint32_t>>
  GetTriangles(const uint32_t eid, const int32_t k) const {
    ASSERT_MSG(UINT32_MAX != edge_info_[eid].first, "invalid edge ID");
    const uint32_t v1 = edge_info_[eid].first;
    const uint32_t v2 = edge_info_[eid].second;
    ASSERT(v1 < n_ && v2 < n_);
    std::vector<std::pair<uint32_t, uint32_t>> triangles;
    // find common neighbors
    size_t p1 = 0, p2 = 0;
    while (p1 < adj_[v1].size() && p2 < adj_[v2].size()) {
      if (adj_[v1][p1].vid == adj_[v2][p2].vid) {
        // the truss numbers of the other two edges should be at least k
        if (k_[adj_[v1][p1].eid] >= k && k_[adj_[v2][p2].eid] >= k) {
          triangles.push_back({adj_[v1][p1].eid, adj_[v2][p2].eid});
        }
        ++p1; ++p2;
      } else if (adj_[v1][p1].vid < adj_[v2][p2].vid) {
        ++p1;
      } else {
        ++p2;
      }
    }
    return triangles;
  }
  // insert an edge (v1, v2) and return its edge ID;
  // the adjacency arrays affected should be sorted later by calling Rectify
  uint32_t LazyInsert(const uint32_t v1, const uint32_t v2) {
    ASSERT_MSG(v1 < n_ && v2 < n_, "invalid insertion");
    ASSERT_MSG(m_ + 1 <= l_, "# of edges exceeded");
    // the ID of (v1, v2)
    const uint32_t eid = free_edges_.back();
    free_edges_.pop_back();
    free_[eid] = false;
    ++m_;
    edge_info_[eid] = {v1, v2};
    // insert the edge to the adjacency arrays
    adj_[v1].push_back({v2, eid});
    adj_[v2].push_back({v1, eid});
    return eid;
  }
  // sort the adjacency arrays
  void Rectify() {
    for (uint32_t v = 0; v < n_; ++v) {
      std::sort(adj_[v].begin(), adj_[v].end(),
                [](const ArrayEntry& ae1, const ArrayEntry& ae2) {
                  return ae1.vid < ae2.vid;
                });
      // no duplicates
      for (size_t i = 1; i < adj_[v].size(); ++i) {
        ASSERT_MSG(adj_[v][i].vid > adj_[v][i - 1].vid,
                   "duplicate edges found");
      }
    }
  }
  // insert an edge (v1, v2) and return its edge ID
  uint32_t Insert(const uint32_t v1, const uint32_t v2) {
    ASSERT_MSG(v1 < n_ && v2 < n_, "invalid insertion");
    ASSERT_MSG(m_ + 1 <= l_, "# of edges exceeded");
    // the ID of (v1, v2)
    const uint32_t eid = free_edges_.back();
    // insert the edge to the adjacency arrays
    size_t p1 = 0;
    while (p1 < adj_[v1].size() && adj_[v1][p1].vid < v2) ++p1;
    ASSERT_MSG(p1 != adj_[v1].size() ? adj_[v1][p1].vid > v2 : true,
               "duplicate insertion");
    adj_[v1].insert(adj_[v1].begin() + p1, {v2, eid});
    size_t p2 = 0;
    while (p2 < adj_[v2].size() && adj_[v2][p2].vid < v1) ++p2;
    ASSERT_MSG(p2 != adj_[v2].size() ? adj_[v2][p2].vid > v1 : true,
               "duplicate insertion");
    adj_[v2].insert(adj_[v2].begin() + p2, {v1, eid});
    // update other information
    edge_info_[eid] = {v1, v2};
    free_edges_.pop_back();
    free_[eid] = false;
    ++m_;
    return eid;
  }
  // remove the edge with ID eid
  void Remove(const uint32_t eid) {
    const uint32_t v1 = edge_info_.at(eid).first;
    const uint32_t v2 = edge_info_.at(eid).second;
    ASSERT_MSG(UINT32_MAX != v1 && UINT32_MAX != v2, "invalid deletion");
    // update information
    free_[eid] = true;
    free_edges_.push_back(eid);
    edge_info_[eid] = {UINT32_MAX, UINT32_MAX};
    // remove the edge from the adjacency arrays
    size_t p1 = 0;
    while (adj_[v1][p1].vid != v2) ++p1;
    adj_[v1].erase(adj_[v1].begin() + p1);
    size_t p2 = 0;
    while (adj_[v2][p2].vid != v1) ++p2;
    adj_[v2].erase(adj_[v2].begin() + p2);
    // decrease the # of edges
    --m_;
  }
  // whether the ID eid is valid
  bool Contain(const uint32_t eid) const {
    return UINT32_MAX != edge_info_.at(eid).first;
  }
  // accessors
  uint32_t n() const { return n_; }
  uint32_t m() const { return m_; }
  uint32_t l() const { return l_; }

 private:
  typedef struct final {
    uint32_t vid;
    uint32_t eid;
  } ArrayEntry;
  // the maximum # of edges this strcture can hold;
  // the current implementation requires l_ < 2^29
  const uint32_t l_;
  // the # of vertices
  const uint32_t n_;
  // the # of edges at the moment; m_ <= l_
  uint32_t m_;
  // the truss numbers of the edges
  const std::vector<int32_t>& k_;
  // free_[i] = true if edge ID i can be allocated; free_.size() == l_
  std::vector<bool> free_;
  // the set of available edge IDs, i.e., the set of IDs i with free_[i] = true
  std::vector<uint32_t> free_edges_;
  // adjacency arrays
  std::vector<std::vector<ArrayEntry>> adj_;
  // edge_info_[i] records the endpoints of the edge with ID i
  std::vector<EdgT> edge_info_;
};
}  // namespace truss_maint

#endif

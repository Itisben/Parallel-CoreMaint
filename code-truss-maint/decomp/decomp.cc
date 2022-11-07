#include "decomp.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <utility>

#define ASSERT(truth) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else

#define ASSERT_MSG(truth, msg) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ << '\n' \
                << "\x1b[1;32mINFO\x1b[0m: " << msg \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else

namespace truss_maint {
namespace decomp {
// for convenience
using std::uint32_t;
// truss decomposition and the corresponding order
Decomp::Decomp(const std::string& file_name) {
  std::ifstream infile(file_name, std::ios::in);
  ASSERT_MSG(infile.is_open(), "cannot open the file");
  // read the size of the graph
  infile >> n_ >> m_;
  ASSERT_MSG(!infile.eof(), "invalid graph file");
  // read the edges
  while (true) {
    uint32_t v1, v2;
    infile >> v1 >> v2;
    if (infile.eof()) break;
    if (v1 > v2) std::swap(v1, v2);
    edges_.push_back({v1, v2});
  }
  infile.close();
  // check its validity
  ASSERT_MSG(edges_.size() == m_, "invalid graph file");
  for (const auto edge : edges_) {
    ASSERT_MSG(edge.first != edge.second, "loops exist in the graph");
    ASSERT_MSG(edge.first < n_ && edge.second < n_, "invalid vertex ID");
  }
  std::sort(edges_.begin(), edges_.end());
  ASSERT_MSG(std::unique(edges_.begin(), edges_.end()) == edges_.end(),
             "multi-edges exist in the graph");
  // initialize adjacency arrays
  adj_.resize(n_);
  for (uint32_t eid = 0; eid < m_; ++eid) {
    const uint32_t v1 = edges_[eid].first;
    const uint32_t v2 = edges_[eid].second;
    adj_[v1].push_back({v2, eid});
    adj_[v2].push_back({v1, eid});
  }
  for (uint32_t vid = 0; vid < n_; ++vid) {
    adj_[vid].shrink_to_fit();
    std::sort(adj_[vid].begin(), adj_[vid].end(),
              [](const ArrayEntry& ae1, const ArrayEntry& ae2) {
                return ae1.vid < ae2.vid;
              });
  }
  // truss decomposition
  // 1. compute the support of each edge by triangle listing
  // 1.1. define a total order over the vertices
  const auto pred = [this](const uint32_t v1, const uint32_t v2) {
    const size_t deg1 = adj_[v1].size();
    const size_t deg2 = adj_[v2].size();
    if (deg1 != deg2) return deg1 > deg2;
    else return v1 > v2;
  };
  // 1.2. sort the vertices in non-ascending order of degree
  std::vector<uint32_t> verts(n_);
  std::iota(verts.begin(), verts.end(), 0);
  std::sort(verts.begin(), verts.end(), pred);
  // 1.3. call the "forward" algorithm to list triangles
  k_.resize(m_, 0);
  std::vector<std::vector<ArrayEntry>> A(n_);
  for (const uint32_t v : verts) {
    for (const auto ae : adj_[v]) {
      const uint32_t u = ae.vid;
      const uint32_t e = ae.eid;
      if (!pred(v, u)) continue;
      size_t pv = 0, pu = 0;
      while (pv < A[v].size() && pu < A[u].size()) {
        if (A[v][pv].vid == A[u][pu].vid) {
          ++k_[A[v][pv].eid]; ++k_[A[u][pu].eid];
          ++k_[e];
          ++pv; ++pu;
        } else if (pred(A[v][pv].vid, A[u][pu].vid)) {
          ++pv;
        } else {
          ++pu;
        }
      }
      A[u].push_back({v, e});
    }
  }
  decltype(A)().swap(A);
  decltype(verts)().swap(verts);
  // 2. decomposition
  // 2.1. sort the edges according to their supports
  const uint32_t max_sup = *std::max_element(k_.cbegin(), k_.cend());
  std::vector<uint32_t> bin(max_sup + 1, 0);
  for (uint32_t eid = 0; eid < m_; ++eid) ++bin[k_[eid]];
  for (uint32_t i = 0, start = 0; i <= max_sup; ++i) {
    start += bin[i];
    bin[i] = start - bin[i];
  }
  std::vector<uint32_t> pos(m_);
  ord_.resize(m_);
  for (uint32_t eid = 0; eid < m_; ++eid) {
    pos[eid] = bin[k_[eid]];
    ord_[pos[eid]] = eid;
    ++bin[k_[eid]];
  }
  std::rotate(bin.rbegin(), bin.rbegin() + 1, bin.rend());
  bin[0] = 0;
  // 2.2. peeling
  rem_.resize(m_, 0); ts_.resize(m_, 0);
  std::vector<bool> removed(m_, false);
  uint32_t k = 0;
  // 2.2.1. process the edges one by one
  for (uint32_t i = 0; i < m_; ++i) {
    k = std::max(k, k_[ord_[i]]);
    ASSERT(bin[k] == i);
    const uint32_t eid = ord_[i];
    ++bin[k_[eid]];
    removed[eid] = true;
    // find triangles containing the edge with ID eid
    std::vector<std::pair<uint32_t, uint32_t>> tris; {
      const uint32_t v1 = edges_[eid].first;
      const uint32_t v2 = edges_[eid].second;
      size_t p1 = 0, p2 = 0;
      while (p1 < adj_[v1].size() && p2 < adj_[v2].size()) {
        if (adj_[v1][p1].vid == adj_[v2][p2].vid) {
          tris.push_back({adj_[v1][p1].eid, adj_[v2][p2].eid});
          ++p1; ++p2;
        } else if (adj_[v1][p1].vid < adj_[v2][p2].vid) {
          ++p1;
        } else {
          ++p2;
        }
      }
    }
    // update rem_[eid] and ts_[eid]
    for (const auto tri : tris) {
      const uint32_t e1 = tri.first;
      const uint32_t e2 = tri.second;
      if (k_[e1] >= k && k_[e2] >= k) ++ts_[eid];
      if (removed[e1] || removed[e2]) continue;
      ++rem_[eid];
      for (const uint32_t e : {e1, e2}) {
        if (k_[e] > k) {
          const uint32_t pe3 = bin[k_[e]];
          const uint32_t pe = pos[e];
          if (pe3 != pe) {
            const uint32_t e3 = ord_[pe3];
            ord_[pe] = e3;
            pos[e3] = pe;
            ord_[pe3] = e;
            pos[e] = pe3;
          }
          ++bin[k_[e]];
          --k_[e];
        }
      }
    }
  }
  // clear adj_
  decltype(adj_)().swap(adj_);
}
void Decomp::WriteToFile(const std::string& file_name) const {
  std::ofstream outfile(file_name, std::ios::binary);
  outfile.write(reinterpret_cast<const char*>(&n_), sizeof n_)
         .write(reinterpret_cast<const char*>(&m_), sizeof m_);
  for (const uint32_t e : ord_) {
    const uint32_t buf[] = {edges_[e].first, edges_[e].second,
                            k_[e], rem_[e], ts_[e]};
    outfile.write(reinterpret_cast<const char*>(buf), sizeof buf);
  }
  outfile.close();
}
}  // namespace decomp
}  // namespace truss_maint

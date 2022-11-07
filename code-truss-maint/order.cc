#include "order.h"

#include <algorithm>
#include <fstream>
#include <tuple>
#include <utility>

#include "defs.h"

namespace truss_maint {
Order::Order(const uint32_t n, const uint32_t l, const std::string& fn)
    : l_(l), n_(n), g_(n_, l_, k_) {
  ASSERT_MSG(64 <= l_ && l_ < (static_cast<uint32_t>(1) << 29),
             "it is required 64 <= l <= 2^29 for the ease of implementation");
  k_    = std::vector<int32_t>(l_, -1);
  node_ = std::vector<ListNode>(l_ + 1);
  new_  = std::vector<bool>(l_, false);
  chg_  = std::vector<bool>(l_, false);
  ts_   = std::vector<uint32_t>(l_, 0);
  s_    = std::vector<uint32_t>(l_, 0);
  // initialize the heap
  HPInit();
  // load index
  LoadIndex(fn);
}
void Order::LoadIndex(const std::string& fn) {
  // read data; no exception handling here
  std::ifstream infile(fn, std::ios::binary);
  uint32_t n = 0, m = 0;
  infile.read(reinterpret_cast<char*>(&n), sizeof n)
        .read(reinterpret_cast<char*>(&m), sizeof m);
  ASSERT(m <= l_ && n_ == n);
  // read the edges and their information
  for (uint32_t e = 0, buf[5]; e < m; ++e) {
    infile.read(reinterpret_cast<char*>(buf), sizeof buf);
    // insert the edge (buf[0], buf[1])
    ASSERT(g_.LazyInsert(buf[0], buf[1]) == e);
    // set the truss number, remaining support, and triangle support
    k_[e] = static_cast<int32_t>(buf[2]);
    node_[e].rem = buf[3]; ts_[e] = buf[4];
    // check
    ASSERT_MSG(e > 0 ? k_[e] >= k_[e - 1] : true, "not in order");
    ASSERT_MSG(buf[3] <= buf[2], "invalid remaining support or truss number");
  }
  infile.close();
  // rectify the graph
  g_.Rectify();
  ASSERT(g_.m() == m);
  // reconstruct the list; @l_ is always the head of the list
  uint32_t prev_e = l_;
  node_[l_].prev = node_[l_].next = UINT32_MAX;
  for (uint32_t e = 0; e < m; ++e) {
    const uint32_t k = k_[e];
    // the sublist of all edges with trussness k
    if (head_.size() <= k) {
      ASSERT(head_.size() == tail_.size());
      head_.resize(k + 1, UINT32_MAX);
      tail_.resize(k + 1, UINT32_MAX);
      head_[k] = e;
    }
    tail_[k] = e;
    // set the list node
    node_[prev_e].next = e;
    node_[e].prev = prev_e;
    node_[e].next = UINT32_MAX;
    node_[e].ext = 0;
    // set e as the new prev_e
    prev_e = e;
  }
  // order maintenance structure
  OMLoad();
}
void Order::Insert(const std::vector<EdgT>& nedges) {
  ASSERT(!nedges.empty());
  // initialization
  std::vector<uint32_t> N;
  for (const auto edge : nedges) {
    const uint32_t e = g_.Insert(edge.first, edge.second);
    N.push_back(e);
    new_[e] = true;
  }
  for (const uint32_t e : N) {
    chg_[e] = true;
    // assume the trussness is -1
    k_[e] = -1;
    // update @s and &ts
    const auto tris = g_.GetTriangles(e);
    s_[e] = ts_[e] = tris.size();
    // update @ext
    for (const auto tri : tris) {
      const uint32_t e1 = tri.first;
      const uint32_t e2 = tri.second;
      if (!new_[e1] && (chg_[e2] || (!new_[e2] && OMPred(e1, e2)))) {
        if (1 == ++node_[e1].ext) HPInsert(e1);
      }
      if (!new_[e2] && (chg_[e1] || (!new_[e1] && OMPred(e2, e1)))) {
        if (1 == ++node_[e2].ext) HPInsert(e2);
      }
    }
    node_[e].ext = 0;
  }
  for (const uint32_t e : N) new_[e] = false;
  // the last processed edge
  uint32_t le = l_;
  // simulate
  for (int32_t k = 0; true; ++k) {
    if (N.empty()) break;
    // expand the head and tail arrays
    ASSERT(head_.size() == tail_.size());
    while (head_.size() <= static_cast<size_t>(k)) {
      head_.push_back(UINT32_MAX);
      tail_.push_back(UINT32_MAX);
    }
    // a stack for DFS
    std::vector<uint32_t> S;
    const auto pred = [this, k](const uint32_t e){return uint32_t(k) >= s_[e];};
    // P1 stores the edges removed in the first shrink
    std::vector<uint32_t> P1;
    // first shrink
    std::copy_if(N.begin(), N.end(), std::back_inserter(S), pred);
    while (!S.empty()) {
      const uint32_t e = S.back(); S.pop_back();
      // update the information of e
      std::tie(node_[e].rem, s_[e], ts_[e]) = std::make_tuple(s_[e], 0, 0);
      chg_[e] = false;
      P1.push_back(e);
      // insert the edge to the order and the list
      OMInsert(e, le);
      ListInsert(e, le);
      le = e;
      // find more edges to shrink and upate the ts values of related edges
      const auto tris = g_.GetTriangles(e);
      for (const auto tri : tris) {
        const uint32_t e1 = tri.first;
        const uint32_t e2 = tri.second;
        // update @ts values
        const int32_t min_k = std::min(k_[e1], k_[e2]);
        const int32_t ori_k = std::min(min_k, k_[e]);
        const int32_t cur_k = std::min(min_k, k);
        if (ori_k < k_[e1] && cur_k >= k_[e1]) ++ts_[e1];
        if (ori_k < k_[e2] && cur_k >= k_[e2]) ++ts_[e2];
        if (cur_k >= k) ++ts_[e];
        // update @s and @ext values
        if (!chg_[e1] && !OMPred(e, e1)) continue;
        if (!chg_[e2] && !OMPred(e, e2)) continue;
        if (chg_[e1]) {
          if (static_cast<uint32_t>(k) == --s_[e1]) S.push_back(e1);
        } else if (chg_[e2] || OMPred(e1, e2)) {
          if (0 == --node_[e1].ext) HPDelete(e1);
        }
        if (chg_[e2]) {
          if (static_cast<uint32_t>(k) == --s_[e2]) S.push_back(e2);
        } else if (chg_[e1] || OMPred(e2, e1)) {
          if (0 == --node_[e2].ext) HPDelete(e2);
        }
      }
      // update the trussness
      k_[e] = k;
    } {
      const auto it = std::copy_if(N.begin(), N.end(), N.begin(),
          [this](const uint32_t e){return chg_[e];});
      N.resize(std::distance(N.begin(), it));
    }
    // update the head and tail
    if (!P1.empty()) {
      head_[k] = P1[0];
      if (UINT32_MAX == tail_[k]) tail_[k] = P1.back();
    }
    // there are edges with trussness k in the heap
    while (hp_tbl_.size() > 1 && k_[HPTop()] == k) {
      const uint32_t e = HPTop(); HPDelete(e);
      const uint32_t s = node_[e].ext + node_[e].rem;
      // e* of Type-2
      if (s > static_cast<uint32_t>(k)) {
        std::tie(s_[e], node_[e].ext) = std::make_tuple(s, 0);
        N.push_back(e);
        chg_[e] = true;
        // update the ext values
        const auto tris = g_.GetTriangles(e);
        for (const auto tri : tris) {
          const uint32_t e1 = tri.first;
          const uint32_t e2 = tri.second;
          if (!chg_[e1] && OMPred(e, e1) && (chg_[e2] || OMPred(e1, e2))) {
            if (1 == ++node_[e1].ext) HPInsert(e1);
          }
          if (!chg_[e2] && OMPred(e, e2) && (chg_[e1] || OMPred(e2, e1))) {
            if (1 == ++node_[e2].ext) HPInsert(e2);
          }
        }
        // remove the edge from the list and the order
        ListRemove(e, head_[k], tail_[k]);
        OMRemove(e);
      } else { // e* of Type-3
        std::tie(node_[e].rem, node_[e].ext) = std::make_tuple(s, 0);
        const auto tris = g_.GetTriangles(e);
        for (const auto tri : tris) {
          const uint32_t e1 = tri.first;
          const uint32_t e2 = tri.second;
          if (chg_[e1] && (chg_[e2] || OMPred(e, e2))) {
            if (static_cast<uint32_t>(k) >= --s_[e1]) S.push_back(e1);
          }
          if (chg_[e2] && (chg_[e1] || OMPred(e, e1))) {
            if (static_cast<uint32_t>(k) >= --s_[e2]) S.push_back(e2);
          }
        }
        // P3 stores the edges removed from the candidate set
        std::vector<uint32_t> P3;
        // remove edges from the candidate set
        while (!S.empty()) {
          const uint32_t ee = S.back(); S.pop_back();
          // update the status of the edge
          std::tie(node_[ee].rem, s_[ee]) = std::make_tuple(s_[ee], 0);
          chg_[ee] = false;  new_[ee] = true;
          P3.push_back(ee);
          if (k_[ee] != k) {
            ts_[ee] = 0;
          }
          const auto tris = g_.GetTriangles(ee);
          for (const auto tri : tris) {
            const uint32_t e1 = tri.first;
            const uint32_t e2 = tri.second;
            if (k_[ee] != k) {
              const int32_t min_k = std::min(k_[e1], k_[e2]);
              const int32_t ori_k = std::min(min_k, k_[ee]);
              const int32_t cur_k = std::min(min_k, k);
              if (ori_k < k_[e1] && cur_k >= k_[e1]) ++ts_[e1];
              if (ori_k < k_[e2] && cur_k >= k_[e2]) ++ts_[e2];
              if (cur_k >= k) ++ts_[ee];
            }
            if (new_[e1] || (!chg_[e1] && !OMPred(e, e1))) continue;
            if (new_[e2] || (!chg_[e2] && !OMPred(e, e2))) continue;
            if (chg_[e1]) {
              if (static_cast<uint32_t>(k) == --s_[e1]) S.push_back(e1);
            } else if (chg_[e2] || OMPred(e1, e2)) {
              if (0 == --node_[e1].ext) HPDelete(e1);
            }
            if (chg_[e2]) {
              if (static_cast<uint32_t>(k) == --s_[e2]) S.push_back(e2);
            } else if (chg_[e1] || OMPred(e2, e1)) {
              if (0 == --node_[e2].ext) HPDelete(e2);
            }
          }
          // update the trussness
          k_[ee] = k;
        }
        // insert the edges in P3 to the order and the list
        le = e;
        for (const uint32_t ee : P3) {
          new_[ee] = false;
          OMInsert(ee, le);
          ListInsert(ee, le);
          le = ee;
        }
        // update the head and the tail
        if (e == tail_[k] && !P3.empty()) tail_[k] = P3.back();
      }
    }
    // shrink N
    const auto it = std::copy_if(N.begin(), N.end(), N.begin(),
        [this](const uint32_t e){return chg_[e];});
    N.resize(std::distance(N.begin(), it));
    // update the last processed edge
    if (UINT32_MAX != tail_[k]) le = tail_[k];
  }
}
void Order::BatchInsert(const std::vector<EdgT>& nedges) {
  ASSERT(nedges.size() > size_t{g_.m()} / 100);
  // initilize the rank
  for (uint32_t r = 0, e = l_; UINT32_MAX != e; e = node_[e].next) {
    rank_[e] = ++r;
  }
  // the candidate set
  std::vector<uint32_t> N;
  for (const auto edg : nedges) N.push_back(g_.Insert(edg.first, edg.second));
  for (const uint32_t e : N) rank_[e] = 0;
  for (const uint32_t e : N) {
    chg_[e] = true;
    // assume the trussness is -1
    k_[e] = -1;
    // update @s and &ts
    const auto tris = g_.GetTriangles(e);
    s_[e] = ts_[e] = tris.size();
    // update @ext
    for (const auto tri : tris) {
      const uint32_t e1 = tri.first;
      const uint32_t e2 = tri.second;
      if (rank_[e1] && (chg_[e2] || rank_[e2] > rank_[e1])) ++node_[e1].ext;
      if (rank_[e2] && (chg_[e1] || rank_[e1] > rank_[e2])) ++node_[e2].ext;
    }
    node_[e].ext = 0;
  }
  // the last processed edge
  uint32_t le = l_;
  // simulate
  for (int32_t k = 0; true; ++k) {
    if (N.empty()) break;
    // expand the head and tail arrays
    ASSERT(head_.size() == tail_.size());
    while (head_.size() <= static_cast<size_t>(k)) {
      head_.push_back(UINT32_MAX);
      tail_.push_back(UINT32_MAX);
    }
    // a stack for DFS
    std::vector<uint32_t> S;
    // P1 stores the edges removed in the first shrink
    std::vector<uint32_t> P1;
    // first shrink
    std::copy_if(N.begin(), N.end(), std::back_inserter(S),
                 [this, k](const uint32_t e){return uint32_t(k) >= s_[e];});
    while (!S.empty()) {
      const uint32_t e = S.back(); S.pop_back();
      // update the information of e
      std::tie(node_[e].rem, s_[e], ts_[e]) = std::make_tuple(s_[e], 0, 0);
      chg_[e] = false;  new_[e] = true;
      P1.push_back(e);
      // insert the edge to the order and the list
      OMInsert(e, le);
      ListInsert(e, le);
      le = e;
      // find more edges to shrink and upate the ts values of related edges
      const auto tris = g_.GetTriangles(e);
      for (const auto tri : tris) {
        const uint32_t e1 = tri.first;
        const uint32_t e2 = tri.second;
        // update @ts values
        const int32_t min_k = std::min(k_[e1], k_[e2]);
        const int32_t ori_k = std::min(min_k, k_[e]);
        const int32_t cur_k = std::min(min_k, k);
        if (ori_k < k_[e1] && cur_k >= k_[e1]) ++ts_[e1];
        if (ori_k < k_[e2] && cur_k >= k_[e2]) ++ts_[e2];
        if (cur_k >= k) ++ts_[e];
        // update @s and @ext values
        if (!chg_[e1] && !(k_[e1] >= k && !new_[e1])) continue;
        if (!chg_[e2] && !(k_[e2] >= k && !new_[e2])) continue;
        if (chg_[e1]) {
          if (static_cast<uint32_t>(k) == --s_[e1]) S.push_back(e1);
        } else if (chg_[e2] || rank_[e2] > rank_[e1]) {
          --node_[e1].ext;
        }
        if (chg_[e2]) {
          if (static_cast<uint32_t>(k) == --s_[e2]) S.push_back(e2);
        } else if (chg_[e1] || rank_[e1] > rank_[e2]) {
          --node_[e2].ext;
        }
      }
      // update the trussness
      k_[e] = k;
    } {
      // shrink the candidate set
      const auto it = std::copy_if(N.begin(), N.end(), N.begin(),
          [this](const uint32_t e){return chg_[e];});
      N.resize(std::distance(N.begin(), it));
    }
    // update the head and tail
    if (!P1.empty()) {
      head_[k] = P1[0];
      if (UINT32_MAX == tail_[k]) tail_[k] = P1.back();
    }
    // reset the new_ array
    for (uint32_t e : P1) new_[e] = false;
    // scan the edges with trussness k
    uint32_t next_e = UINT32_MAX;
    for (uint32_t e = node_[le].next; true; e = next_e) {
      if (UINT32_MAX == e || k_[e] > k) break;
      next_e = node_[e].next;
      // e* of Type-1
      if (0 == node_[e].ext) continue;
      const uint32_t s = node_[e].ext + node_[e].rem;
      // e* of Type-2
      if (s > static_cast<uint32_t>(k)) {
        std::tie(s_[e], node_[e].ext) = std::make_tuple(s, 0);
        N.push_back(e);
        chg_[e] = true;
        // update the ext values
        const auto tris = g_.GetTriangles(e);
        for (const auto tri : tris) {
          const uint32_t e1 = tri.first;
          const uint32_t e2 = tri.second;
          if (rank_[e] < rank_[e1] && (chg_[e2] || rank_[e1] < rank_[e2])) {
            ++node_[e1].ext;
          }
          if (rank_[e] < rank_[e2] && (chg_[e1] || rank_[e2] < rank_[e1])) {
            ++node_[e2].ext;
          }
        }
        // remove the edge from the list and the order
        OMRemove(e);
        ListRemove(e, head_[k], tail_[k]);
      } else { // e* of Type-3
        std::tie(node_[e].rem, node_[e].ext) = std::make_tuple(s, 0);
        const auto tris = g_.GetTriangles(e);
        for (const auto tri : tris) {
          const uint32_t e1 = tri.first;
          const uint32_t e2 = tri.second;
          if (chg_[e1] && (chg_[e2] || rank_[e] < rank_[e2])) {
            if (static_cast<uint32_t>(k) >= --s_[e1]) S.push_back(e1);
          }
          if (chg_[e2] && (chg_[e1] || rank_[e] < rank_[e1])) {
            if (static_cast<uint32_t>(k) >= --s_[e2]) S.push_back(e2);
          }
        }
        // P3 stores the edges removed from the candidate set
        std::vector<uint32_t> P3;
        // remove edges from the candidate set
        le = e;
        while (!S.empty()) {
          const uint32_t ee = S.back(); S.pop_back();
          // update the status of the edge
          std::tie(node_[ee].rem, s_[ee]) = std::make_tuple(s_[ee], 0);
          chg_[ee] = false;
          P3.push_back(ee);
          // list the edge to the list
          OMInsert(ee, le);
          ListInsert(ee, le);
          le = ee;
          if (k_[ee] != k) {
            ts_[ee] = 0;
          }
          const auto tris = g_.GetTriangles(ee);
          for (const auto tri : tris) {
            const uint32_t e1 = tri.first;
            const uint32_t e2 = tri.second;
            if (k_[ee] != k) {
              const int32_t min_k = std::min(k_[e1], k_[e2]);
              const int32_t ori_k = std::min(min_k, k_[ee]);
              const int32_t cur_k = std::min(min_k, k);
              if (ori_k < k_[e1] && cur_k >= k_[e1]) ++ts_[e1];
              if (ori_k < k_[e2] && cur_k >= k_[e2]) ++ts_[e2];
              if (cur_k >= k) ++ts_[ee];
            }
            if (!chg_[e1] && rank_[e1] <= rank_[e]) continue;
            if (!chg_[e2] && rank_[e2] <= rank_[e]) continue;
            if (chg_[e1]) {
              if (static_cast<uint32_t>(k) == --s_[e1]) S.push_back(e1);
            } else if (chg_[e2] || rank_[e1] < rank_[e2]) {
              --node_[e1].ext;
            }
            if (chg_[e2]) {
              if (static_cast<uint32_t>(k) == --s_[e2]) S.push_back(e2);
            } else if (chg_[e1] || rank_[e2] < rank_[e1]) {
              --node_[e2].ext;
            }
          }
          // update the trussness
          k_[ee] = k;
        }
        // update the head and the tail
        if (e == tail_[k] && !P3.empty()) tail_[k] = P3.back();
      }
    }
    // shrink N
    const auto it = std::copy_if(N.begin(), N.end(), N.begin(),
        [this](const uint32_t e){return chg_[e];});
    N.resize(std::distance(N.begin(), it));
    // update the last processed edge
    if (UINT32_MAX != tail_[k]) le = tail_[k];
  }
  // reset the rank
  for (uint32_t e = l_; UINT32_MAX != e; e = node_[e].next) {
    rank_[e] = UINT32_MAX;
  }
}
void Order::Remove(const uint32_t v1, const uint32_t v2) {
  const uint32_t re = g_.Get(v1, v2);
  std::vector<uint32_t> S;
  // update the @ts and $rem values for other related edges
  const auto tris = g_.GetTriangles(re);
  for (const auto tri : tris) {
    const uint32_t e1 = tri.first;
    const uint32_t e2 = tri.second;
    // update @ts values
    const int32_t min_k = std::min({k_[re], k_[e1], k_[e2]});
    if (min_k >= k_[e1] && --ts_[e1] < uint32_t(k_[e1])) S.push_back(e1);
    if (min_k >= k_[e2] && --ts_[e2] < uint32_t(k_[e2])) S.push_back(e2);
    // update @rem values
    uint32_t min_e = re;
    if (OMPred(e1, min_e)) min_e = e1;
    if (OMPred(e2, min_e)) min_e = e2;
    --node_[min_e].rem;
  }
  // remove the edge
  g_.Remove(re);
  OMRemove(re);
  ListRemove(re, head_[k_[re]], tail_[k_[re]]);
  k_[re] = -1;
  ts_[re] = 0;
  // propagate
  while (!S.empty()) {
    const uint32_t e = S.back(); S.pop_back();
    std::tie(k_[e], ts_[e]) = std::make_pair(k_[e] - 1, 0);
    // the previous edge of @e in the new position of order
    const uint32_t prev_e = node_[head_[k_[e] + 1]].prev;
    // update @ts and @rem
    const auto tris = g_.GetTriangles(e);
    for (const auto tri : tris) {
      const uint32_t e1 = tri.first;
      const uint32_t e2 = tri.second;
      // the original trussness of the triangle
      const int32_t ori_k = std::min({k_[e] + 1, k_[e1], k_[e2]});
      // update the @ts value for @e
      if (ori_k >= k_[e]) ++ts_[e];
      // update the @ts values for other edges
      if (ori_k > k_[e]) {
        if (ori_k == k_[e1] && ts_[e1]-- == uint32_t(k_[e1])) S.push_back(e1);
        if (ori_k == k_[e2] && ts_[e2]-- == uint32_t(k_[e2])) S.push_back(e2);
      }
      // update the @rem values
      uint32_t ori_min_e = e;
      if (OMPred(e1, ori_min_e)) ori_min_e = e1;
      if (OMPred(e2, ori_min_e)) ori_min_e = e2;
      if (e != ori_min_e && OMPred(prev_e, ori_min_e)) {
        --node_[ori_min_e].rem;
        ++node_[e].rem;
      }
    }
    // remove @e from the order and the list
    OMRemove(e);
    ListRemove(e, head_[k_[e] + 1], tail_[k_[e] + 1]);
    // reinsert @e to the order and the list
    OMInsert(e, prev_e);
    ListInsert(e, prev_e);
    tail_[k_[e]] = e;
    if (UINT32_MAX == head_[k_[e]]) head_[k_[e]] = e;
  }
}
void Order::BatchRemove(const std::vector<EdgT>& redges) {
  // a stack
  std::vector<uint32_t> S;
  // indicators if an edge is in S
  auto& inS = new_;
  std::vector<uint32_t> reids;
  for (const auto edg : redges) {
    const uint32_t re = g_.Get(edg.first, edg.second);
    reids.push_back(re);
    inS[re] = true;
  }
  // remove the edges from the graph
  size_t i = 0;
  for (const auto edg : redges) {
    const uint32_t re = reids[i++];
    const auto tris = g_.GetTriangles(re);
    for (const auto tri : tris) {
      const uint32_t e1 = tri.first;
      const uint32_t e2 = tri.second;
      // update @ts values
      const int32_t min_k = std::min({k_[re], k_[e1], k_[e2]});
      if (min_k >= k_[e1]) --ts_[e1];
      if (ts_[e1] < uint32_t(k_[e1]) && !inS[e1]) {
        S.push_back(e1);
        inS[e1] = true;
      }
      if (min_k >= k_[e2]) --ts_[e2];
      if (ts_[e2] < uint32_t(k_[e2]) && !inS[e2]) {
        S.push_back(e2);
        inS[e2] = true;
      }
      // update @rem values
      uint32_t min_e = re;
      if (OMPred(e1, min_e)) min_e = e1;
      if (OMPred(e2, min_e)) min_e = e2;
      --node_[min_e].rem;
    }
    // remove the edge
    g_.Remove(re);
    OMRemove(re);
    ListRemove(re, head_[k_[re]], tail_[k_[re]]);
    k_[re] = -1;
    ts_[re] = 0;
  }
  // update trussnesses
  while (!S.empty()) {
    const uint32_t e = S.back(); S.pop_back(); inS[e] = false;
    // @tri_ts stores (ts, tri) pairs
    std::vector<std::pair<int32_t, std::pair<uint32_t, uint32_t>>> tri_ts;
    // @bin used for bin sort
    std::vector<uint32_t> bin(k_[e] + 1, 0);
    // enumerate triangles
    uint32_t v1, v2;
    std::tie(v1, v2) = g_.Get(e);
    const auto tris = g_.GetTriangles(e);
    for (const auto tri : tris) {
      const uint32_t e1 = tri.first;
      const uint32_t e2 = tri.second;
      if (uint32_t(k_[e1]) < ts_[e] || uint32_t(k_[e2]) < ts_[e]) continue;
      const int32_t min_k = std::min(k_[e1], k_[e2]);
      if (min_k > k_[e]) continue;
      tri_ts.push_back({min_k, {e1, e2}});
      ++bin[min_k];
    } { // bin sort: initialization
      uint32_t start = 0;
      for (int32_t i = 0; i <= k_[e]; ++i) {
        const uint32_t cnt = bin[i];
        bin[i] = start;
        start += cnt;
      }
    }
    // bin sort: finalization
    auto s_tri_ts = tri_ts;
    for (const auto p : tri_ts) {
      s_tri_ts[bin[p.first]] = p;
      ++bin[p.first];
    }
    const int32_t ori_k = k_[e];
    // compute new @k_[e] and @ts_[e]
    int32_t ptr = s_tri_ts.size() - 1;
    while (ptr >= 0 && s_tri_ts[ptr].first == k_[e]) --ptr;
    uint32_t prev_e = tail_[k_[e]];
    while (uint32_t(k_[e]) > ts_[e]) {
      --k_[e];
      while (ptr >= 0 && s_tri_ts[ptr].first >= k_[e]) {
        --ptr;
        ++ts_[e];
      }
      if (UINT32_MAX != head_[k_[e] + 1]) {
        prev_e = node_[head_[k_[e] + 1]].prev;
      }
    }
    for (int32_t i = s_tri_ts.size() - 1; i > ptr; --i) {
      const int32_t ori_tri_k = std::min(ori_k, s_tri_ts[i].first);
      const int32_t cur_tri_k = std::min(k_[e], s_tri_ts[i].first);
      const uint32_t e1 = s_tri_ts[i].second.first;
      const uint32_t e2 = s_tri_ts[i].second.second;
      if (ori_tri_k >= k_[e1] && cur_tri_k < k_[e1]) --ts_[e1];
      if (ts_[e1] < uint32_t(k_[e1]) && !inS[e1]) {
        S.push_back(e1);
        inS[e1] = true;
      }
      if (ori_tri_k >= k_[e2] && cur_tri_k < k_[e2]) --ts_[e2];
      if (ts_[e2] < uint32_t(k_[e2]) && !inS[e2]) {
        S.push_back(e2);
        inS[e2] = true;
      }
      // update the @rem values
      uint32_t min_e = e;
      if (OMPred(e1, min_e)) min_e = e1;
      if (OMPred(e2, min_e)) min_e = e2;
      if (e != min_e && OMPred(prev_e, min_e)) {
        --node_[min_e].rem;
        ++node_[e].rem;
      }
    }
    // remove @e from the order and the list
    OMRemove(e);
    ListRemove(e, head_[ori_k], tail_[ori_k]);
    // reinsert @e to the order and the list
    OMInsert(e, prev_e);
    ListInsert(e, prev_e);
    tail_[k_[e]] = e;
    if (UINT32_MAX == head_[k_[e]]) head_[k_[e]] = e;
  }
  // clean
  for (const uint32_t reid : reids) inS[reid] = false;
}
void Order::Debug() const {
  // check @chg_, @new_, @s_
  for (uint32_t e = 0; e < l_; ++e) {
    ASSERT(!chg_.at(e) && !new_.at(e));
    ASSERT(s_.at(e) == 0);
  }
  // check the heap
  ASSERT(1 == hp_tbl_.size());
  for (uint32_t e = 0; e < l_; ++e) {
    ASSERT(UINT32_MAX == hp_pos_.at(e));
  }
  // auxiliary structure
  uint32_t seen_cnt = 0;
  std::vector<bool> seen(l_, false);
  // check the list, i.e., @node_
  for (uint32_t e = node_.at(l_).next; UINT32_MAX != e; e = node_.at(e).next) {
    ASSERT(!seen.at(e));
    seen.at(e) = true;
    ++seen_cnt;
    // the edge indeed exists in the graph
    ASSERT(g_.Contain(e));
    // the trussness should be non-decreasing
    if (UINT32_MAX != node_.at(e).next) {
      ASSERT(k_.at(e) <= k_.at(node_.at(e).next));
    }
    // check the list node
    if (UINT32_MAX != node_.at(e).next) {
      ASSERT(node_.at(node_.at(e).next).prev == e);
    }
    if (UINT32_MAX != node_.at(e).prev) {
      ASSERT(node_.at(node_.at(e).prev).next == e);
    }
    ASSERT(0 == node_.at(e).ext);
    // check the remaining support
    const auto tris = g_.GetTriangles(e);
    uint32_t check_rem = 0;
    for (const auto tri : tris) {
      if (!seen.at(tri.first) && !seen.at(tri.second)) ++check_rem;
    }
    ASSERT(node_.at(e).rem == check_rem);
    ASSERT(node_.at(e).rem <= static_cast<uint32_t>(k_.at(e)));
  }
  for (uint32_t e = 0; e < l_; ++e) {
    ASSERT(seen.at(e) || !g_.Contain(e));
  }
  ASSERT(g_.m() == seen_cnt);
  ASSERT(UINT32_MAX == node_.at(l_).prev);
  // check the head and tail arrays
  ASSERT(head_.size() == tail_.size());
  uint32_t check_k_cnt = 0;
  for (int32_t k = 0; k < static_cast<int32_t>(head_.size()); ++k) {
    if (UINT32_MAX == head_.at(k)) {
      ASSERT(UINT32_MAX == tail_.at(k));
    } else {
      ASSERT(UINT32_MAX != tail_.at(k));
      for (uint32_t e = head_.at(k); tail_.at(k) != e; e = node_.at(e).next) {
        ASSERT(k_.at(e) == k && seen.at(e));
        ++check_k_cnt;
      }
      ASSERT(k_.at(tail_.at(k)) == k && seen.at(tail_.at(k)));
      ++check_k_cnt;
      if (l_ != node_.at(head_.at(k)).prev) {
        ASSERT(k_.at(node_.at(head_.at(k)).prev) < k);
      }
      if (UINT32_MAX != node_.at(tail_.at(k)).next) {
        ASSERT(k_.at(node_.at(tail_.at(k)).next) > k);
      }
    }
  }
  ASSERT(g_.m() == check_k_cnt);
  // check the order
  for (uint32_t e = l_; UINT32_MAX != e; e = node_.at(e).next) {
    if (UINT32_MAX != node_.at(e).next) {
      ASSERT(OMPred(e, node_.at(e).next));
    }
  }
  uint32_t check_grp_num = 0;
  for (uint32_t e = l_; UINT32_MAX != e; ) {
    const uint32_t grp = om_grp_.at(e);
    uint32_t check_cnt = 0;
    for (; UINT32_MAX != e; e = node_.at(e).next) {
      if (grp != om_grp_.at(e)) break;
      ++check_cnt;
    }
    ++check_grp_num;
    ASSERT(om_cnt_.at(grp) == check_cnt);
  }
  // check the free groups
  uint32_t check_free_grp_num = 0;
  for (uint32_t g = om_avail_; UINT32_MAX != g; g = om_nodes_.at(g).next) {
    ++check_free_grp_num;
    ASSERT(0 == om_cnt_.at(g));
    if (UINT32_MAX != om_nodes_.at(g).next) {
      ASSERT(om_nodes_.at(om_nodes_.at(g).next).prev == g);
    }
  }
  ASSERT(UINT32_MAX != om_avail_ && UINT32_MAX == om_nodes_.at(om_avail_).prev);
  // check the used groups
  uint32_t check_used_grp_num = 0;
  for (uint32_t g = 0; UINT32_MAX != g; g = om_nodes_.at(g).next) {
    ++check_used_grp_num;
    ASSERT(l_ + 1 != g ? 0 < om_cnt_.at(g) : 0 == om_cnt_.at(g));
    ASSERT(om_cnt_.at(g) <= om_grp_ub_);
    if (UINT32_MAX != om_nodes_.at(g).next) {
      ASSERT(om_nodes_.at(om_nodes_.at(g).next).tag > om_nodes_.at(g).tag);
      ASSERT(om_nodes_.at(om_nodes_.at(g).next).prev == g);
    }
  }
  ASSERT(UINT32_MAX == om_nodes_.at(0).prev);
  // sum of the # of free groups and the # of used groups
  ASSERT(check_used_grp_num == check_grp_num + 1);
  ASSERT(check_free_grp_num + check_used_grp_num == l_ + 2);
  // check the ts values
  for (uint32_t e = node_.at(l_).next; UINT32_MAX != e; e = node_.at(e).next) {
    uint32_t check_ts = 0;
    const auto tris = g_.GetTriangles(e);
    for (const auto tri : tris) {
      const uint32_t e1 = tri.first;
      const uint32_t e2 = tri.second;
      if (k_.at(e1) >= k_.at(e) && k_.at(e2) >= k_.at(e)) {
        ++check_ts;
      }
    }
    ASSERT(check_ts == ts_.at(e));
  }
  for (uint32_t e = 0; e < l_; ++e) {
    ASSERT(!seen.at(e) ? 0 == ts_.at(e) : true);
  }
  printf("Debug completed.\n");
}
void Order::Check(const std::string& fn) const {
  // load the ground truth
  std::ifstream ansfile(fn, std::ios::binary);
  // read the # of vertices and the # of edges
  uint32_t n = -1, m = -1;
  ansfile.read(reinterpret_cast<char*>(&n), sizeof n)
         .read(reinterpret_cast<char*>(&m), sizeof m);
  ASSERT(g_.m() == m && n_ == n);
  // read the edges and their truss numbers
  std::vector<std::pair<EdgT, int32_t>> answer;
  for (uint32_t e = 0; e < m; ++e) {
    uint32_t buf[5];
    ansfile.read(reinterpret_cast<char*>(buf), sizeof buf);
    uint32_t v1 = buf[0];
    uint32_t v2 = buf[1];
    ASSERT(v1 < n && v2 < n);
    if (v1 > v2) std::swap(v1, v2);
    answer.push_back({{v1, v2}, buf[2]});
  }
  ansfile.close();
  // get the results computed by ours
  std::vector<std::pair<EdgT, int32_t>> result;
  for (uint32_t e = node_.at(l_).next; UINT32_MAX != e; e = node_.at(e).next) {
    const auto info = g_.Get(e);
    uint32_t v1 = info.first;
    uint32_t v2 = info.second;
    if (v1 > v2) std::swap(v1, v2);
    result.push_back({{v1, v2}, k_.at(e)});
  }
  // compare
  std::sort(answer.begin(), answer.end());
  std::sort(result.begin(), result.end());
  ASSERT_MSG(answer == result, "wrong answer");
}
// list maintenance
void Order::ListInsert(const uint32_t e1, const uint32_t e2) {
  node_[e1].next = node_[e2].next;
  node_[e1].prev = e2;
  node_[e2].next = e1;
  if (UINT32_MAX != node_[e1].next) node_[node_[e1].next].prev = e1;
}
void Order::ListRemove(const uint32_t e, uint32_t& head, uint32_t& tail) {
  if (head == tail) {
    head = tail = UINT32_MAX;
  } else if (e == head) {
    head = node_[e].next;
  } else if (e == tail) {
    tail = node_[e].prev;
  }
  if (UINT32_MAX != node_[e].prev) {
    node_[node_[e].prev].next = node_[e].next;
  }
  if (UINT32_MAX != node_[e].next) {
    node_[node_[e].next].prev = node_[e].prev;
  }
  node_[e].prev = node_[e].next = UINT32_MAX;
}
// order maintenance
void Order::OMLoad() {
  ASSERT(l_ >= 64 && l_ < (static_cast<uint32_t>(1) << 29));
  om_nodes_ = std::vector<OMNode>(l_ + 2);
  om_tag_ = std::vector<uint64_t>(l_ + 1, 0);
  om_grp_ = std::vector<uint32_t>(l_ + 1, 0);
  om_cnt_ = std::vector<uint32_t>(l_ + 2, 0);
  // available group ids
  om_avail_ = 0;
  om_nodes_[0].prev = om_nodes_[l_].next = UINT32_MAX;
  for (uint32_t i = 1; i <= l_; ++i) {
    om_nodes_[i].prev = i - 1;
    om_nodes_[i - 1].next = i;
  }
  // bulkload
  const uint64_t step = (static_cast<uint64_t>(1) << 60) / (l_ + 1);
  uint32_t cnt = 0;
  uint32_t tgid = UINT32_MAX;
  for (uint32_t p = l_; UINT32_MAX != p; p = node_[p].next) {
    if (unlikely(0 == cnt)) {
      // initialize the group
      tgid = om_avail_;
      om_avail_ = om_nodes_[om_avail_].next;
      om_nodes_[om_avail_].prev = UINT32_MAX;
      om_nodes_[tgid].next = l_ + 1;
      if (0 == tgid) {
        om_nodes_[tgid].prev = UINT32_MAX;
        om_nodes_[tgid].tag = 0;
      } else {
        om_nodes_[tgid].prev = tgid - 1;
        om_nodes_[tgid - 1].next = tgid;
        om_nodes_[tgid].tag = om_nodes_[tgid - 1].tag + step;
      }
      om_cnt_[tgid] = 1;
      // initialize the edge
      om_tag_[p] = 0;
      om_grp_[p] = tgid;
    } else {
      om_tag_[p] = om_tag_[node_[p].prev] + (static_cast<uint64_t>(1) << 30);
      om_grp_[p] = tgid;
      ++om_cnt_[tgid];
    }
    if (++cnt == om_grp_ub_ / 2) cnt = 0;
  }
  // the end mark of the order
  om_cnt_[l_ + 1] = 0;
  om_nodes_[l_ + 1].tag = (static_cast<uint64_t>(1) << 60) - 1;
  om_nodes_[l_ + 1].next = UINT32_MAX;
  om_nodes_[l_ + 1].prev = tgid;
}
void Order::OMInsert(const uint32_t e1, const uint32_t e2) {
  // the group is full; that is, a new group needs to be created
  if (om_cnt_[om_grp_[e2]] == om_grp_ub_) {
    // the new group id
    const uint32_t ngid = om_avail_;
    om_avail_ = om_nodes_[om_avail_].next;
    om_nodes_[om_avail_].prev = UINT32_MAX;
    // the group to which e2 belongs and the next group
    const uint32_t gid1 = om_grp_[e2];
    const uint32_t gid2 = om_nodes_[gid1].next;
    // there is no vacancy; relabeling is needed
    if (om_nodes_[gid1].tag + 1 == om_nodes_[gid2].tag) {
      const uint64_t tag = om_nodes_[gid1].tag;
      uint32_t cnt = 1;
      uint32_t threshold = 1;
      uint64_t mask = UINT64_MAX;
      uint32_t ph = gid1, pe = gid1;
      // find the "smallest" enclosing tag range with low enough density
      while (cnt >= threshold) {
        threshold = threshold << 1;
        mask = mask << 2;
        while (UINT32_MAX != om_nodes_[ph].prev &&
               (om_nodes_[om_nodes_[ph].prev].tag & mask) == (tag & mask)) {
          ph = om_nodes_[ph].prev;
          ++cnt;
        }
        while (UINT32_MAX != om_nodes_[pe].next &&
               (om_nodes_[om_nodes_[pe].next].tag & mask) == (tag & mask)) {
          pe = om_nodes_[pe].next;
          ++cnt;
        }
      }
      const uint64_t step = ((~mask) + 1) / cnt;
      om_nodes_[ph].tag = (tag & mask);
      uint32_t p = om_nodes_[ph].next;
      while (--cnt > 0) {
        om_nodes_[p].tag = om_nodes_[om_nodes_[p].prev].tag + step;
        p = om_nodes_[p].next;
      }
    }
    ASSERT(om_nodes_[gid1].tag + 1 < om_nodes_[gid2].tag);
    // insert the new group to the list
    om_nodes_[ngid].prev = gid1;
    om_nodes_[ngid].next = gid2;
    om_nodes_[gid1].next = ngid;
    om_nodes_[gid2].prev = ngid;
    // assign the next tag instead of the middle of gid1 and gid2
    om_nodes_[ngid].tag = om_nodes_[gid1].tag + 1;
    // distribute the edges
    // step 1. find the first edge in the group
    uint32_t p = e2;
    while (UINT32_MAX != node_[p].prev && gid1 == om_grp_[node_[p].prev]) {
      p = node_[p].prev;
    }
    uint32_t cnt = 0;
    // step 2. assign new tags to the edges remaining in the original group
    om_tag_[p] = 0;
    p = node_[p].next;
    while (++cnt < om_grp_ub_ / 2) {
      om_tag_[p] = om_tag_[node_[p].prev] + (static_cast<uint64_t>(1) << 30);
      p = node_[p].next;
    }
    om_cnt_[gid1] = om_grp_ub_ / 2;
    // step 3. distribute the edges
    om_tag_[p] = 0;
    om_grp_[p] = ngid;
    p = node_[p].next;
    while (++cnt < om_grp_ub_) {
      om_tag_[p] = om_tag_[node_[p].prev] + (static_cast<uint64_t>(1) << 30);
      om_grp_[p] = ngid;
      p = node_[p].next;
    }
    om_cnt_[ngid] = om_grp_ub_ - om_grp_ub_ / 2;
    // in case @om_grp_ub_ is set as an odd number
    ASSERT(om_cnt_[ngid] == om_cnt_[gid1]);
  }
  // whether to relabel the edges in group @om_grp_[e2]
  bool relabel = false;
  const uint32_t e3 = node_[e2].next;
  if (UINT32_MAX == e3 || om_grp_[e3] != om_grp_[e2]) {
    relabel = (om_tag_[e2] >= (static_cast<uint64_t>(1) << 60));
  } else {
    relabel = (om_tag_[e2] + 1 == om_tag_[e3]);
  }
  if (relabel) {
    const uint32_t gid1 = om_grp_[e2];
    uint32_t p = e2;
    while (UINT32_MAX != node_[p].prev && gid1 == om_grp_[node_[p].prev]) {
      p = node_[p].prev;
    }
    om_tag_[p] = 0;
    p = node_[p].next;
    while (UINT32_MAX != p && gid1 == om_grp_[p]) {
      om_tag_[p] = om_tag_[node_[p].prev] + (static_cast<uint64_t>(1) << 30);
      p = node_[p].next;
    }
  }
  // set the new edge
  ++om_cnt_[om_grp_[e2]];
  om_grp_[e1] = om_grp_[e2];
  if (UINT32_MAX == e3 || om_grp_[e3] != om_grp_[e2]) {
    om_tag_[e1] = om_tag_[e2] + (static_cast<uint64_t>(1) << 30);
  } else {
    om_tag_[e1] = (om_tag_[e2] + om_tag_[e3]) / 2;
  }
}
void Order::OMRemove(const uint32_t e) {
  const uint32_t gid = om_grp_[e];
  ASSERT(om_cnt_[gid] >= 1);
  if (1 == om_cnt_[gid]) {
    om_nodes_[om_nodes_[gid].prev].next = om_nodes_[gid].next;
    om_nodes_[om_nodes_[gid].next].prev = om_nodes_[gid].prev;
    // update the available groups
    om_nodes_[gid].prev = UINT32_MAX;
    om_nodes_[gid].next = om_avail_;
    om_nodes_[om_avail_].prev = gid;
    om_avail_ = gid;
  }
  --om_cnt_[gid];
}
bool Order::OMPred(const uint32_t e1, const uint32_t e2) const {
  if (om_grp_[e1] == om_grp_[e2]) {
    return om_tag_[e1] < om_tag_[e2];
  } else {
    return om_nodes_[om_grp_[e1]].tag < om_nodes_[om_grp_[e2]].tag;
  }
}
// heap maintenance
void Order::HPInit() {
  ASSERT(0 < l_ && l_ < (static_cast<uint32_t>(1) << 29));
  hp_pos_.resize(l_ + 1, UINT32_MAX);
  hp_tbl_.push_back(UINT32_MAX);
}
void Order::HPUp(const uint32_t h, const uint32_t e) {
  uint32_t c = h;
  uint32_t p = c / 2;
  while (0 != p && OMPred(e, hp_tbl_[p])) {
    hp_tbl_[c] = hp_tbl_[p];
    hp_pos_[hp_tbl_[p]] = c;
    c = p;
    p = c / 2;
  }
  hp_tbl_[c] = e;
  hp_pos_[e] = c;
}
void Order::HPDown(const uint32_t h, const uint32_t e) {
  const uint32_t size = hp_tbl_.size();
  uint32_t p = h;
  uint32_t c = p * 2;
  while (c < size) {
    if (c + 1 < size && OMPred(hp_tbl_[c + 1], hp_tbl_[c])) ++c;
    if (OMPred(e, hp_tbl_[c])) break;
    hp_tbl_[p] = hp_tbl_[c];
    hp_pos_[hp_tbl_[c]] = p;
    p = c;
    c = p * 2;
  }
  hp_tbl_[p] = e;
  hp_pos_[e] = p;
}
void Order::HPInsert(const uint32_t e) {
  ASSERT(UINT32_MAX == hp_pos_.at(e));
  const uint32_t size = hp_tbl_.size();
  hp_tbl_.push_back(e);
  HPUp(size, e);
}
void Order::HPDelete(const uint32_t e) {
  ASSERT(UINT32_MAX != hp_pos_.at(e));
  const uint32_t size = hp_tbl_.size() - 1;
  const uint32_t e2 = hp_tbl_[size];
  hp_tbl_.pop_back();
  // shift down or up
  if (size != hp_pos_[e]) {
    const uint32_t h = hp_pos_[e];
    if (h / 2 != 0 && OMPred(e2, hp_tbl_[h / 2])) {
      HPUp(h, e2);
    } else {
      HPDown(h, e2);
    }
  }
  hp_pos_[e] = UINT32_MAX;
}
uint32_t Order::HPTop() const {
  ASSERT(hp_tbl_.size() > 1);
  return hp_tbl_[1];
}
}  // namespace truss_maint

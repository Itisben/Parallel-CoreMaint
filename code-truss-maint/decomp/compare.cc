#include <cassert>
#include <cstdint>
#include <fstream>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

using std::uint32_t;
using EdgeT = std::pair<uint32_t, uint32_t>;
using InfoT = std::tuple<uint32_t, uint32_t, uint32_t>;

int main(int /* argc */, char** argv) {
  std::map<EdgeT, InfoT> answer; {
    std::ifstream afile(argv[1], std::ios::in);
    uint32_t n, m, l;
    afile >> n >> m >> l;
    // space
    std::vector<EdgeT> edges;
    std::vector<uint32_t> k(m);
    std::vector<uint32_t> rem(m);
    std::vector<uint32_t> ts(m);
    // read edges and truss numbers
    for (uint32_t i = 0; i < m; ++i) {
      uint32_t v1, v2;
      afile >> v1 >> v2 >> k[i];
      if (v1 > v2) std::swap(v1, v2);
      edges.push_back({v1, v2});
    }
    for (uint32_t i = 0; i < m; ++i) {
      afile >> rem[i] >> ts[i];
    }
    // create answer
    for (uint32_t i = 0; i < m; ++i) {
      answer[edges[i]] = {k[i], 0, ts[i]};
    }
    afile.close();
  }
  std::map<EdgeT, InfoT> result; {
    std::ifstream rfile(argv[2], std::ios::binary);
    uint32_t n, m;
    rfile.read(reinterpret_cast<char*>(&n), sizeof n)
         .read(reinterpret_cast<char*>(&m), sizeof m);
    for (uint32_t buf[5], i = 0; i < m; ++i) {
      rfile.read(reinterpret_cast<char*>(buf), sizeof buf);
      if (buf[0] > buf[1]) std::swap(buf[0], buf[1]);
      result[std::make_pair(buf[0], buf[1])] = {buf[2], 0, buf[4]};
    }
    rfile.close();
  }
  // compare
  assert(answer == result);
}

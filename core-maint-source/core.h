#ifndef CORE_CORE_H_
#define CORE_CORE_H_

#include <vector>

namespace core {
class CoreMaintenance {
 public:
  virtual ~CoreMaintenance() {}

  virtual void ComputeCore(const std::vector<std::vector<int>>& graph,
                           const bool init_idx, // initialize the index?
                           std::vector<int>& core) = 0;
  virtual void Insert(const int v1, const int v2,
                      std::vector<std::vector<int>>& graph,
                      std::vector<int>& core) = 0;
  virtual void Remove(const int v1, const int v2,
                      std::vector<std::vector<int>>& graph,
                      std::vector<int>& core) = 0;
  virtual void Check(const std::vector<std::vector<int>>& graph,
                     const std::vector<int>& core) const = 0;
};
}  // namespace core

#endif

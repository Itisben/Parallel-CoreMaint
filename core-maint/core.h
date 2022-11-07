#ifndef CORE_CORE_H_
#define CORE_CORE_H_

#include <vector>
#include <assert.h>
namespace core {
class CoreMaintenance {
 public:
  virtual ~CoreMaintenance() {}

  virtual void ComputeCore( std::vector<std::vector<int>>& graph,
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


/*for debug*/
extern int cnt_edge;
extern int cnt_edge2;
extern int cnt_Vs; 
extern int cnt_Vp;
extern int cnt_S;
extern int cnt_gap_adjust;
extern int global_check_num;  
#endif

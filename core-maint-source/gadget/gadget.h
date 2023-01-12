#ifndef CORE_GADGET_GADGET_H_
#define CORE_GADGET_GADGET_H_

#include <utility>
#include <vector>

namespace gadget {
void RepeatWith(const char symbol, const int repeat);
std::vector<std::vector<int>> ReadGraph(const char* const path,
                                        int* const n, int* const m);
std::vector<std::pair<int, int>> ReadTempEdgesS(const char* const path,
                                                int* const n, int* const m);
std::vector<std::pair<int, int>> ReadEdgesS(const char* const path,
                                            int* const n, int* const m);
}  // namespace gadget

#endif

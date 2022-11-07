#include "decomp.h"

int main(int /* argc */, char** argv) {
  // read the graph and truss-decompose it
  truss_maint::decomp::Decomp index(argv[1]);
  index.WriteToFile(argv[2]);
}

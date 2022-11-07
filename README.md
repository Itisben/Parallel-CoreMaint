This is the paralle vesion of order-based core maintenance. The experiments test two methods, the edge-join method and our order-based method.

# The Structure of the Folder
1. ```~/core-maint/``` : include the source code.
* **main.cc**: the main program 
* **defs.h**: this file defines some useful macros
* **core.h**: this header file defines the base class for core maintenance algorithm
* **gadget**: this folder contains gadgets for reading files and implements several data structures
    * **disjoint.h, disjoint.cc**: a union-find-delete data structure with constant time deletion complexity
    * **gadget.h, gadget.cc**: for reading data files, the format of which will be explained in details later
    * **heap.h, heap.cc**: a minimum heap
    * **sb_tree.h, sb_tree.cc**: the size-balanced tree
    * **treap.h, treap.cc**: the treap (tree + heap)
* **traversal**: this folder contains the code for the traversal core maintenance algorithm
* **gm_graph**: this folder is to generate a library that can support loading the CSR format of graphs. 
* **glist**: this folder contains the code for the order-based core maintenance algorithm
* **ours-csr-new**: this folder contains the code for our sequential and parlalel core maintenance algorithm; it also includes the paralell OM source files that used for our parallel core maintenance method.

2. ```~/experiment/```: include all the experiment scripts.
* **ParCM-Test.sh**: the shell script that can run the experiment and get the test result as csv files.
* **ParCM-parse-output.py**: the python script called by benchmark.sh to parse the test result.
* **results**: this folder contains all the tested results. 
* **convert**: this folder contains all the python script that transfer the test results to the figures and table that can be included in the paper.


# Compiling and Executing

1. The ubuntu 18.04 with g++ are used for compiling.

1. Compile the gm_graph runtime library and the scc binary by running make in the directory ```~/core-maint/gm_graph```.

1. Compile the source code by running make in the directory ```~/core-maint```.

1. To execute application, in the directory ```~/core-main```, run:

* For insert edges: ```./core -p <graph_path> -I <number_edges> -m <method> -T <graph_format> -w <number of workers>```
* For remove edges: ```./core -p <graph_path> -R <number_edges> -m <method> -T <graph_format> -w <number of workers>```
* method includes (-m): 4 our parallel method
* graph format (-T): 3 csr sample edges, 4 csr without sample edges, 5 csr with repeated random
* worker number (-w): from 1 to the number of pysical cores, __0 means sequential running__.


# Executing with Script

1. The experiment can run with a script ```~/experiments/ParCM-Test.sh```. Before running, you should renew the paths of all tested graphs. We can vary the range of workers and the range of inserted/removed edges.  


# Graph Conversion Tool
1 **~/experiment/convert/convert**: this file can convert the graphs in edges lists into csv format. There is a graph conversion utility to convert the adjacency list or edge list graph files into the binary format used by the current application. 



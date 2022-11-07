
This is the source code for the parallel version of order-based core maintenance for the draft paper ``Parallel Order-Based Core Maintenance in Dynamic Graphs'' (https://arxiv.org/abs/2210.14290). The experiments test two methods, the edge-join method and our order-based method.

# The Structure of the Folder
1. ```~/core-maint/``` : include the source code.
   * **main.cc**: the main program 
   * **defs.h**: this file defines some useful macros
   * **core.h**: this header file defines the base class for core maintenance algorithm
   * **gadget**: this folder contains gadgets for reading files and implements several data structures
       * **disjoint.h, disjoint.cc**: a union-find-delete data structure with constant time deletion complexity
       * **gadget.h, gadget.cc**: for reading data files, the format of which will be explained in detail later
       * **heap.h, heap.cc**: a minimum heap
       * **sb_tree.h, sb_tree.cc**: the size-balanced tree
       * **treap.h, treap.cc**: the treap (tree + heap)
   * **traversal**: this folder contains the code for the traversal core maintenance algorithm
   * **gm_graph**: this folder is to generate a library that can support loading the CSR format of graphs. 
   * **glist**: this folder contains the code for the order-based core maintenance algorithm
   * **ours-csr-new**: this folder contains the code for our sequential and parallel core maintenance algorithm; it also includes the parallel OM source    files used for our parallel core maintenance method.


2. ```~/experiment/```: include all the experiment scripts.
   * **ParCM-Test.sh**: the shell script that can run the experiment and get the test result as CSV files.
   * **ParCM-parse-output.py**: the python script called by benchmark.sh to parse the test result.
   * **results**: this folder contains all the tested results. 
   * **convert**: this folder contains all the Python scripts that transfer the test results to the figures and table in the paper.

3. ```~/FasterCoreMaintenance```: the Edge Join method. 

4. ```~/ParaCoM-master```: the other compared method.


# Compiling and Executing

1. Ubuntu 18.04 with g++  11.3.0 is used for compiling.

1. Compile the gm_graph runtime library running make in the directory ```~/core-maint/gm_graph```.

1. Compile the source code by running make in the directory ```~/core-maint```.

1. To execute our method, in the directory ```~/core-main```, run:

   * For insert edges: ```./core -p <graph_path> -I <number_edges> -m <method> -T <graph_format> -w <number of workers>```
   * For remove edges: ```./core -p <graph_path> -R <number_edges> -m <method> -T <graph_format> -w <number of workers>```
   * method includes (-m): 4 our parallel method
   * graph format (-T): 3 CSR sample edges, 4 CSR without sample edges, 5 CSR with repeated random
   * worker number (-w): from 1 to the number of physical cores, __0 means sequential running__.

1. To execute the compared edge-join method, in the directory ```~/FasterCoreMaintenance``, run:
* run make to cimple
* run ```~/FasterCoreMaintenance/kcore <graph-file> <edge-file> <workers>``
* <edge-file> is the edges for insertion or removal. This file can be automatically generated when executing our method of inserting or removing edges. 


# Executing with Script

1. The experiment can run with a script ```~/experiments/ParCM-Test.sh```. Before running, you should renew the paths of all tested graphs. We can vary the range of workers and the range of inserted/removed edges.  


# Graph Conversion Tool
1 **~/experiment/convert/convert**: this file can convert the graphs in edge lists into CSR format. There is a graph conversion utility to convert the adjacency list or edge list graph files into the binary format used by the current application. 



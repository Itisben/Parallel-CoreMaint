# README #



## File Organization ##

* **order.h, order.cc**: implementation of truss maintenance
* **defs.h**: some useful macros
* **graph.h**: a class for graph management
* **test.cc**: a sample program
* **decomp/decomp.o, decomp/decomp.cc**: implementation of truss decomposition
* **decomp/sample.cc**: a sample program for truss decomposition

## How to Use the Code? ##

* First, use "decomp" to construct an index for the graph before update.
* Then, using this index as input, you can construct an "Order" object (order.h
  and order.cc) to handle update edges.

## Data Format ##

* The graph file is of the following format:

    n  m

    u_1 v_1

    u_2 v_2

    ...

    u_m v_m

where n is the # of vertices, m is the # of edges, (u_i, v_i) is the i-th edge.

* The file containing update edges is of the following format:

    m

    u_1 v_1

    u_2 v_2

    ...

    u_m v_m

where m is the # of update edges.

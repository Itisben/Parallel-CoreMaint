# Parallel core maintenance (ParaCoM) 

The core number is a fundamental index reflecting the cohesiveness of a graph, which are widely used in large-scale graph analytics. The core maintenance problem requires to update the core numbers of vertices after a set of edges and vertices are inserted into or deleted from the graph. We investigate the parallelism in the core update process when multiple edges and vertices are inserted or deleted. Specifically, we discover a structure called superior edge set, the insertion or deletion of edges in which can be processed in parallel. 

# How To Use

    make

    ./mykcore -p(parallel)/-c(centralized) graph_filename edge_filename

The graph file is supposed to be in the format of edge file and the vertices are number from 0, see the examplegraph.

# Example

    ./mykcore -p examplegraph examplegraph_10

run our parallel algorithm, the time cost is saved in examplegraph_time_p.txt

    ./mykcore -c examplegraph examplegraph_10

run the centralized algorithm, the time cost is saved in examplegraph_time_c.txt

The new core numbers are saved in examplegraph_core_del.txt and examplegraph_core_ins.txt for deletion and insertion case respectively.

1: extend the sequential version by csr format to test on the large graphs. 
    1.1 the code should be consistend with the paper. 

2: extend to the parallel version.
  the main idea is that 
  (1) when A come accrose a node that is locked by B, A should be asume that this node is black (with around 90% it is black), which can avoid the block chain. 
  (2) when the backward happen, the node is locked for backward.   
  (3) the work is O(|E+|log |E+|), but the depth is optimized from O(|E+|log |E+|) to O(E'+ log|E'+|), where E'+ is for unit insertion and E+ is for batch insertion in parallel.

  one node can do the backward once by own worker. Other worker assume it is black. 

3: A new parallel algorithm is deviced. 

4: the k-truss can also be applied. 



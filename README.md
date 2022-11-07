# ParCoreMaint

This is the paralle vesion of order-based core maintenance. 

## concurrent Order Maintenance
* we must use the two level tag to inplement the O(1) time for insertion, removing, and order checking. 
* we propose a efficient concurrent version of Order Data Structure. this has little related work.


## concurrent Core Maintenance
* based on the Order Maintenance, we propose a concurrent order.
* we have to count the size of locked vertices for each edge insertion. 
* reduce the **length** of relable the linked list of vertices.
* each core number can maintain independently. 
* show the maximum sequential length. 


## same idea can to parallel the k-trucc maintainence?



## More Thread May reduce the waiting time?


## Test
test with 1, 2, 4, 8, 16, 32, 64 cores.

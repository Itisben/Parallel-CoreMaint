# a script that tranfer the graph format from bin to edge and edge-100000
# that can be used to test compared parallel methods.
# import required module
import os
import re
# assign directory
directory = '/home/guob15/test-graph/test'
 

for filename in os.scandir(directory):
    mystr = str(filename.path)
    m = re.search(".bin$", mystr)
    if m:
        print(mystr)
        line = "./core -p " + mystr + " -T 5 -I 1000000 -t 30"
        os.system(line)
        print("")
        print("")



# import required module
import os
import re
# assign directory
directory = '/home/guob15/test-graph/test'
 
# iterate over files in
# that directory
for filename in os.listdir(directory):
    print(filename)
    name = ""
    name = filename
    match = re.search("*.bin", "thsi.bin")
    #print(match)
    # if match:
    #     f = os.path.join(directory, filename)
    #     # checking if it is a file
    #     if os.path.isfile(f):
    #         print(filename)
    #         line = "./core -p " + filename + " -T 5 -I 100000"
    #         os.system("./core")
            
        
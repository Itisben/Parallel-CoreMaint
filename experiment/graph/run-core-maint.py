import os


core_path = "../../core-maint/core"
#graph_file_path = "/home/guob15/Documents/test"
graph_file_path = "/home/guo/Documents/test"
folder = "."
#parameter = " -I 4 -m 3, -T 2"
parameter = " -T 3 -m 3 -c 2"

def TestLocal():
    print(graph_file_path)
    for root, dirs, files in os.walk(graph_file_path):
        print(files)
        for filename in files:
            #if filename[0] == 'g': continue
            #path = "./" + filename + parameter
            
            path = graph_file_path + "/" + filename
            print(path)
            cmd = core_path + " -p " + path + parameter
            
            print(cmd)
            os.system(cmd)
            print("\n")
            print("\n")

TestLocal()




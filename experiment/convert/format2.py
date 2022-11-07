#format the graph as the input graph for the experiments.

#filepath = "./out.dblp_coauthor"
filepath = "./a.txt"
filepath2 = "./dblp-temp.txt"

edges = []
with open(filepath) as fp:
    cnt = 0
    for line in fp:
        #print("line {} contents {}".format(cnt, line))
        for x in line:
            print(x)
        # v0 = int(a[0])
        # v1 = int(a[1])
        # v3 = int(a[3])
        
        # if v3 == 0: continue

        # if v0 == v1: continue
        
        # edges.append((v0, v1))


fp = open(filepath2, "w+")
for e in edges:
    fp.write(str(e[0]) + ' ' + str(e[1]) + '\n')
fp.close()




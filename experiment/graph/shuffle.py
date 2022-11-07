#format the graph as the input graph for the experiments.
import random
filepath = "/home/guo/Documents/TestG/pokec.txt"
filepath2 = "/home/guo/Documents/TestG/pokec-shuffle.txt"
haszero = True

edgeset = []
i = 0
n = 0; m = 0
with open(filepath) as fp:
    for line in fp:
        #print("line {} contents {}".format(cnt, line))
        a = line.split(' ')
        v0 = int(a[0])
        v1 = int(a[1])
        if i is 0: n = v0; m = v1
        else: edgeset.append((v0, v1))
        i+=1

random.shuffle(edgeset)
fp = open(filepath2, "w")
head = str(n) + ' ' + str(m) + '\n'
fp.write(head)
for e in edgeset:
    fp.write(str(e[0]) + ' ' + str(e[1]) + '\n')
fp.close()






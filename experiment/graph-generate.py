#this is the graph generator by SNAP stanford. 
import snap

n = 1000000
deg = 8
m = n * 8
path = "/home/guo/Documents/SCC-test-data/"
def ER():
    Graph = snap.GenRndGnm(snap.PNGraph, n, m)

    with open(path + 'er-1m-8m.txt', 'a') as f:
        for EI in Graph.Edges():
            #print("edge: (%d, %d)" % (EI.GetSrcNId(), EI.GetDstNId()))
            f.write( str(EI.GetSrcNId()) + ' ' + str(EI.GetDstNId()) + '\n')

def BA():
    Rnd = snap.TRnd()
    UGraph = snap.GenPrefAttach(n, deg, Rnd)
    i = 0
    with open(path + 'ba-1m-8m.txt', 'a') as f:
        for EI in UGraph.Edges():
            #print("edge: (%d, %d)" % (EI.GetSrcNId(), EI.GetDstNId())) 
            f.write("%d %d\n" % (EI.GetSrcNId(), EI.GetDstNId()) )
            i += 1
    print ("num: %d" % (i))

def SW():
    Rnd = snap.TRnd(1,0)
    UGraph1 = snap.GenSmallWorld(n, deg, 0, Rnd)
    i = 0
    with open(path + 'sw-1m-8m.txt', 'a') as f:
        for EI in UGraph1.Edges():
            f.write("%d %d\n" % (EI.GetSrcNId(), EI.GetDstNId()) )
            i += 1
    print("num: %d" % (i))

def RMAT():
    Rnd = snap.TRnd()
    Graph = snap.GenRMat(n, m, .6, .1, .15, Rnd)
    s = set()
    with open(path + 'rmat-1m-8m.txt', 'a') as f:
        for EI in Graph.Edges():
            a = EI.GetSrcNId()
            b = EI.GetDstNId()
            if a > b: 
                a,b=b,a
            if (a,b) not in s:  
                f.write("%d %d\n" % (a,b) )
                s.add((a,b))
SW()
#BA()
#ER()
#RMAT()

# def main():
#     ER()

# if __name__ == "__main__":
#     main()

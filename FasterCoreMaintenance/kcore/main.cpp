#include <ctime>
#include <iostream>
#include "kcore.h"

using namespace std;

int main(int argc,char* argv[])
{

    if(argc != 4) {
        fprintf(stderr,"Usage: %s <graph-file> <edge-file> <threads>\n",argv[0]);
        return -1;
    }
    
    int thread = atoi(argv[3]);
    if (thread <= 0) thread = 1;
    Decomposition dec(thread);

    if(!dec.Load(argv[1],argv[2])) {
        fprintf(stderr,"Load file failed\n");
        return -1;
    }

    dec.Delete();
    dec.Insert();
}

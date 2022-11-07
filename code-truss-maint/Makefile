CC=g++ -march=native -O3
CFLAGS=-c -I. -std=c++0x -Wfatal-errors

all: tm

tm: test.o order.o
	$(CC) test.o order.o -o tm
	rm *.o

test.o: test.cc
	$(CC) $(CFLAGS) test.cc -o test.o

order.o: order.cc
	$(CC) $(CFLAGS) order.cc -o order.o

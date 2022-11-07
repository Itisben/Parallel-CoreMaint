CC=g++ -march=native -O3
CFLAGS=-c -I. -std=c++11 -Wfatal-errors

all: sample

sample: sample.o decomp.o
	$(CC) sample.o decomp.o -o sample
	rm *.o

sample.o: sample.cc
	$(CC) $(CFLAGS) sample.cc -o sample.o

decomp.o: decomp.cc
	$(CC) $(CFLAGS) decomp.cc -o decomp.o

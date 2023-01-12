CC=g++ -march=native -O3
CFLAGS=-c -I. -std=c++0x -Wfatal-errors

all: core

core: main.o gadget.o heap.o treap.o traversal.o glist.o
	$(CC) main.o gadget.o heap.o treap.o traversal.o glist.o -o core
	rm *.o

test: test.o gadget.o heap.o treap.o traversal.o glist.o
	$(CC) test.o gadget.o heap.o treap.o traversal.o glist.o -o test
	rm *.o

main.o: main.cc
	$(CC) $(CFLAGS) main.cc -o main.o

test.o: test.cc
	$(CC) $(CFLAGS) test.cc -o test.o

gadget.o: gadget/gadget.cc
	$(CC) $(CFLAGS) gadget/gadget.cc -o gadget.o

heap.o: gadget/heap.cc
	$(CC) $(CFLAGS) gadget/heap.cc -o heap.o

treap.o: gadget/treap.cc
	$(CC) $(CFLAGS) gadget/treap.cc -o treap.o

traversal.o: traversal/traversal.cc
	$(CC) $(CFLAGS) traversal/traversal.cc -o traversal.o

glist.o: glist/glist.cc
	$(CC) $(CFLAGS) glist/glist.cc -o glist.o


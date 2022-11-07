#pragma once
#include <gm.h>
#include <omp.h>
#include <vector>
#include <set>
#include <queue>
#include <math.h>       /* pow */
#include <climits>
#include <atomic>
#include <assert.h>
#include <stdio.h>

using namespace std;

//atomic operations
#if 0
#define atomic_read(v)      (*(volatile typeof(*v) *)(v))
#define atomic_write(v,a)   (*(volatile typeof(*v) *)(v) = (a))


//The valid memory order variants are __ATOMIC_RELAXED, __ATOMIC_SEQ_CST, __ATOMIC_ACQUIRE, and __ATOMIC_CONSUME.
#define atomic_read(v)      __atomic_load_n(v)
#define atomic_write(v, a)     __atomic_store_n(v, a)


#else


#define atomic_read(v)      __atomic_load_n(v, __ATOMIC_RELAXED)
#define atomic_write(v, a)     __atomic_store_n(v, a, __ATOMIC_RELAXED)
#endif 


#define fetch_or(a, b)      __sync_fetch_and_or(a,b)
#define cas(a, b, c)        __sync_bool_compare_and_swap(a,b,c)
#define casval(a, b, c)     __sysc_val_compare_and_swap(a, c, c)
#define atomic_add(a, b)        __sync_fetch_and_add(a, b)
#define atomic_sub(a, b)        __sync_fetch_and_sub(a, b)
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

#define clear_cache(x, y) __builtin___clear_cache(x, y)

const size_t MAX_WORKER_NUM = 128;

// open mp atomic 
// #define omp_atomic_update(a, b) {
//     #pragma omp atomic update 
//         a

// #define omp_atomic_read(a) #pragma omp atomic read \
//                                 a

// #define omp_atomic_write(a) #pragma omp atomic write 
//                                a

//for graph
const size_t MAX_CORE_NUMBER = 1024; //max increased core number, should be large enough
const size_t CORE_BUFFER_SIZE = 1024;
const size_t VECTOR_CAPACITY = 1024*1024;

typedef int node_t; // the vertex of graph
const node_t NONE = -1; 
const int EMPTY = -1;
typedef int edge_t; // the edge of graph
typedef int core_t; // the core number of graph 
typedef int deg_t;  // degree
//typedef int id_t;  // vertex id, group id, worker id. -1 is NONE 
typedef int color_t;
typedef int label_t;
typedef size_t ver_t;
typedef unsigned int lock_t;

typedef unsigned int counter_t;

// tag is the label for OM data structure. 
typedef unsigned long long int tag_t; // 64bit label tag;
typedef unsigned int subtag_t; // 32bit label subtag;
typedef int group_t;
const size_t MAX_TAG = 0xffffffffffffffff; //64 bit
const size_t INIT_TAG_GAP = 0xffffffff; // 32 bit unsigned integer
const size_t MAX_SUBTAG = 0xffffffff; // 32 bit unsigned integer
const size_t INIT_SUBTAG = 0xefffffff; // 2^31 -1
const size_t MAX_CORE_SIZE = 0xffffffff; //32 bit unsigned integer
const size_t MIN_GROUP_SIZE = 15; //(log n)/2

#define CAS_LOCK 0
#define OMP_LOCK 1
#define UNLOCK 0
#define LOCK 1

#if 1
//assert 
// @ASSERT will print the code line and file where the assert is false, say,
// @truth is false. After printing the code infomation, it will halt the
// program.
#define ASSERT(truth) \
    if (!(truth)) { \
      printf("\x1b[1;31mASSERT\x1b[0m, LINE:%d, FILE:%s\n", \
             __LINE__, __FILE__); \
      exit(EXIT_FAILURE); \
    } else

// Similar to @ASSERT, but @ASSERT_INFO will print another line called INFO.
// The output contents are specified by @info.
#define ASSERT_INFO(truth, info) \
    if (!(truth)) { \
      printf("\x1b[1;31mASSERT\x1b[0m, LINE:%d, FILE:%s\n", \
             __LINE__, __FILE__); \
      printf("\x1b[1;32mINFO\x1b[0m: %s\n", info); \
      exit(EXIT_FAILURE); \
    } else

// @ERROR will print @msg and then according to the truth value of @to_exit,
// it will halt the program or not. To be specific, when @to_exit is true,
// it will halt the program, otherwise, not.
#define ERROR(msg, to_exit) \
    if (true) { \
      printf("\x1b[1;31mERROR\x1b[0m: %s\n", msg); \
      if (to_exit) { \
        exit(EXIT_FAILURE); \
      } \
    } else  
#else

#define ASSERT(truth) 
#define ASSERT_INFO(truth, info)
#define ERROR(msg, to_exit) 

#endif








/*this stack is more efficient than the STL stack
* push to tail and pop the head
* ??? this may out of the bound */
// class QUEUE {
// private:
//     std::vector<int> array_;
//     int head_; 
//     int tail_; // append position.
// public:
//     QUEUE(std::size_t size) { 
//         head_ = 0; tail_ = 0;
//         reserve(size);
//     }

//     QUEUE() { 
//         head_ = 0; tail_ = 0;
//         //reserve(size);
//     }

//     ~QUEUE() {}
//     void reserve(std::size_t size) {array_.reserve(size);}
//     inline bool empty() { return head_ >= tail_;}
//     inline void clear() { head_ = 0; tail_ = 0; array_.clear();}
//     inline void push(int v) { array_.push_back(v); tail_++; }
//     inline void pop() { head_++; }
//     inline void poptail(){array_.pop_back(); tail_--;}
//     inline void remove() {}
//     inline int top() { return array_[head_]; }
//     inline int size() {return tail_ - head_;}
// };



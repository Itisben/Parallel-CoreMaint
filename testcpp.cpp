// atomic::load/store example
#include <iostream>       // std::cout
#include <atomic>         // std::atomic, std::memory_order_relaxed
#include <thread>         // std::thread
#include <vector>
#include <assert.h>

int x, y;
std::atomic<bool> ready{false};

void init()
{
    // for(int i = 0; i < 100; i++){
    //     printf("init %d\n", i);
    // }
  x = 2;
  y = 3;
  //atomic_thread_fence(std::memory_order_release);
  ready.store(true, std::memory_order_relaxed);
}
void use()
{
  while(!ready.load(std::memory_order_relaxed)){}

  {
    //atomic_thread_fence(std::memory_order_acquire);
    std::cout << x + y << std::endl;;
  }
}

int main ()
{

  
//   std::thread second (set_foo,10);
//   std::thread first (print_foo);
//   first.join();
//   second.join();


    std::thread first(init);
    std::thread second(use);

    first.join();
    second.join();
  return 0;
}
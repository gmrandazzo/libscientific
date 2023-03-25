#include <stdio.h>
#include <math.h>
#include "numeric.h"
#include <time.h>

void test3()
{
  size_t i;
  double v;
  puts("Test3: Generate 10M random numbers between -2.8 and 3.8");
  srand_(time(NULL));
  for(i = 0; i < 10000000; i++){
    v = randDouble(-2.8, 3.8);
    if(v > -2.8 && v < 3.8){
      continue;
    }
    else{
      printf("ERROR: [%lu] %f\n", i,  v);
      abort();
    }
  }
}

void test2()
{
  int i;
  puts("Test2: Random double between -1.5 and 2.3");
  srand_(1245);
  for(i = 0; i < 5; i++)
    printf("%.4f\n", randDouble(-1.5, 2.3));
}

void test1()
{
  int i;
  puts("Test 1: Random int between 1 10");
  srand_(12345);
  for(i = 0; i < 5; i++)
    printf("%d\n", randInt(1, 10));
}

int main()
{
  test1();
  test2();
  test3();
  return 0;
}

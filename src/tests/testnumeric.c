#include <stdio.h>
#include <math.h>
#include "numeric.h"
#include <time.h>

void test2()
{
  int i;
  puts("Test2: Random double between -1.5 and 2.3");
  srand(1245);
  for(i = 0; i < 5; i++)
    printf("%.4f\n", randDouble(-1.5, 2.3));
}

void test1()
{
  int i;
  puts("Test 1: Random int between 1 10");
  srand(12345);
  for(i = 0; i < 5; i++)
    printf("%d\n", randInt(1, 10));
}

int main()
{
  test1();
  test2();
  return 0;
}

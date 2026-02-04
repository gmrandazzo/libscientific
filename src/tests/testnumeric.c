/* Unit tests for the numeric module.
 * Copyright (C) 2020-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <math.h>
#include "numeric.h"
#include <time.h>

/* MISSING TESTS
 * void StochasticUniversalSample
 * void RouletteWheelselection
 * void Combinations
 */
void test5()
{
  puts("Test5: factorial of a number");
  if(FLOAT_EQ(Factorial(12), 479001600, 1e-1)){
    puts("Factorial OK!");
  }
  else{
    puts("Factorial ERROR!");
    abort();
  }
}

void test4()
{
  puts("Test4: Random number generator");
  srand_(12345);
  size_t i;
  double r[] = {3239304236.000000,
                1907938550.0000000000,
                3743148142.0000000000,
                2778857055.0000000000,
                3762851757.0000000000,
                795925487.0000000000,
                340349905.0000000000,
                1721361429.0000000000,
                2420684656.0000000000,
                1428642710.0000000000};
  for(i = 0; i < 10; i++){
    double a = rand_();
    if(FLOAT_EQ(a, r[i], 1e-2))
      continue;
    else{
      printf("error random number generator rand_: %f != %f\n", a, r[i]);
      abort();
    }
  }
}

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
      printf("ERROR: [%zu] %f\n", i,  v);
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
  test4();
  test5();
  return 0;
}

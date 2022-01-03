/* numeric.c
*
* Copyright (C) <2016>  Giuseppe Marco Randazzo
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "interpolate.h"


double missing_value(){ return MISSING; }

/* Random Generator
 * No thread safe....
 */
//#ifdef WIN32
int myrand_r (unsigned int *seed)
{
  unsigned int next = *seed;
  int result;

  next *= 1103515245;
  next += 12345;
  result = (unsigned int) (next / 65536) % 2048;

  next *= 1103515245;
  next += 12345;
  result <<= 10;
  result ^= (unsigned int) (next / 65536) % 1024;

  next *= 1103515245;
  next += 12345;
  result <<= 10;
  result ^= (unsigned int) (next / 65536) % 1024;

  *seed = next;

  return result;
}
//#endif

/* Xorshift 128 is thread safe */
uint32_t xor128(void) {
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t;

  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = w ^ (w >> 19) ^ (t ^ (t >> 8));
}

size_t Factorial(size_t x)
{
  size_t f;
  size_t i;
  f = x;
  i = x-1;

  while(i > 0){
    f *= i;
    i--;
  }
  return f;
}

int randInt(int low, int high)
{
  return (int) (xor128() % ((high) - low) + low);
}

double randDouble(double low, double high)
{
  /*
   * xor128() cannot return 4294967296
   */
  double range = (high - low);
  double div = 4294967296.0 / range;
  return low + (xor128() / div);
}

double square(double x){ return x*x; }

void StochasticUniversalSample(dvector *fitness, size_t nselect, size_t init, uivector **selection)
{
  size_t i, k;
  double sum, ptr;
  k = 0;
  sum = 0.f;
  for(i = 0; i <  fitness->size; i++){
    sum += getDVectorValue(fitness, i);
  }

  /*CRITICAL BUG! IF ALL VALUES ARE < 0 then there is an infinite loop
  if(FLOAT_EQ(sum, 0.f, 1e-3) || sum < 0){
    srand(time(0));
    do{
      k = rand() % (fitness->size-1);
      if(FLOAT_EQ(getDVectorValue(fitness, k), 0.f, 1e-3) || getDVectorValue(fitness, k) < 0){
        continue;
      }
      else{
        UIVectorAppend(selection, k);
        break;
      }
    }while(1);
  }
  else{*/
    sum = 0.f;
    srand(init);
    ptr = (double)rand()/(double)RAND_MAX;
    for(i = 0; i < fitness->size; i++){
      if(k < nselect){
        sum += getDVectorValue(fitness, i);
        while(sum > ptr){
          // select individual k;
          UIVectorAppend(selection, i);
          ptr += 1;
          k++;
        }
      }
      else{
        break;
      }
    }
  /*}*/
}

void RouletteWheelselection(dvector *fitness, size_t nselect, size_t init, uivector **selection)
{
  size_t i, k;
  double sum, sumfitness, ptr;
  k = 0;
  sum = 0.f;
  sumfitness = 0.f;
  srand(init);

  for(i = 0; i < fitness->size; i++){
    sumfitness += getDVectorValue(fitness, i);
  }

  for(k = 0; k < nselect; k++){
    ptr = (double)rand()/(double)RAND_MAX;
    for(i = 0; i < fitness->size; i++){
      sum += getDVectorValue(fitness, i);
      if(sum > ptr){
        UIVectorAppend(selection, i);
        break;
      }
      else{
        continue;
      }
    }
  }
}

void Combinations(uivector *num, matrix **comb)
{
  size_t temp;
  size_t i, j;
  for(j = 1; j <= num->size; j++){
    for(i = 0; i < num->size-1; i++){
      temp = num->data[i];
      num->data[i] = num->data[i+1];
      num->data[i+1] = temp;
      MatrixAppendUIRow(comb, num);
    }
  }
}

/*
 * area of the curve via the  trapezoid rule
 */
double curve_area(matrix *xy, size_t intervals)
{
  size_t i;
  double base, height;
  double area = 0.f;
  matrix *interp_xy;
  initMatrix(&interp_xy);
  /*If intervals > 0 interpolate with natural cubic splines to have more fine area */
  if(intervals > 0){
    /* If two points have different y but share same x the algorithm will fail*/
    interpolate(xy, intervals, &interp_xy);
  }
  else{
    MatrixCopy(xy, &interp_xy);
    intervals = xy->row;
  }

  for(i = 0; i < intervals-1; i++){
      /* Trapezoidal method */
      base = interp_xy->data[i+1][0]-interp_xy->data[i][0];
      height = ((interp_xy->data[i][1]+interp_xy->data[i+1][1])/2.f);
      area += base*height;
  }
  DelMatrix(&interp_xy);
  return area;
}

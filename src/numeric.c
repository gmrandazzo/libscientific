/* Provides numerical analysis methods and solvers.
 * Copyright (C) 2016-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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

uint32_t generate_seed(uint32_t seed) {
    uint32_t a = 0x7AFB2C23ULL;
    uint32_t c = 0x894C3ULL;
    uint32_t new_seed = (a * seed + c);
    return new_seed;
}

/* The state must be initialized to non-zero */
uint32_t xorshift128(struct xorshift128_state *state)
{
	/* Algorithm "xor128" from p. 5 of Marsaglia, "Xorshift RNGs" */
	uint32_t t  = state->x[3];
  uint32_t s  = state->x[0];  /* Perform a contrived 32-bit shift. */
	state->x[3] = state->x[2];
	state->x[2] = state->x[1];
	state->x[1] = s;

	t ^= t << 11;
	t ^= t >> 8;
	return state->x[0] = t ^ s ^ (s >> 19);
}

uint32_t XOR128_SEED = 0;

void srand_(uint32_t seed)
{
  XOR128_SEED = generate_seed(seed);
}

double rand_()
{
  struct xorshift128_state state;
  if(XOR128_SEED  == 0)
    XOR128_SEED = time(NULL);
  state.x[0] = XOR128_SEED;
  state.x[1] = XOR128_SEED ^ 0x5a7b96158bd42e27ULL;
  state.x[2] = XOR128_SEED ^ 0x3a8e9f2baf7e592bULL;
  state.x[3] = XOR128_SEED ^ 0x0b243e4b4b2aa8d3ULL;
  XOR128_SEED = generate_seed(XOR128_SEED);
  return xorshift128(&state);
}

int randInt(int low, int high)
{
  struct xorshift128_state state;
  if(XOR128_SEED  == 0)
    XOR128_SEED = time(NULL);
  state.x[0] = XOR128_SEED;
  state.x[1] = XOR128_SEED ^ 0x5a7b96158bd42e27ULL;
  state.x[2] = XOR128_SEED ^ 0x3a8e9f2baf7e592bULL;
  state.x[3] = XOR128_SEED ^ 0x0b243e4b4b2aa8d3ULL;
  XOR128_SEED = generate_seed(XOR128_SEED);
  return (int) (xorshift128(&state) % ((high) - low) + low);
}

double randDouble(double low, double high)
{
  /*
   * xor128() cannot return 4294967296
   */
  struct xorshift128_state state;
  if(XOR128_SEED  == 0)
    XOR128_SEED = time(NULL);
  state.x[0] = XOR128_SEED;
  state.x[1] = XOR128_SEED ^ 0x5a7b96158bd42e27ULL;
  state.x[2] = XOR128_SEED ^ 0x3a8e9f2baf7e592bULL;
  state.x[3] = XOR128_SEED ^ 0x0b243e4b4b2aa8d3ULL;
  XOR128_SEED = generate_seed(XOR128_SEED);
  double range = (high - low);
  double div = 4294967296.0 / range;
  return low + (xorshift128(&state) / div);
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

double square(double x){ return x*x; }

void StochasticUniversalSample(dvector *fitness, size_t nselect, size_t init, uivector *selection)
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
    srand_(init);
    ptr = (double)rand_()/(double)RAND_MAX;
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

void RouletteWheelselection(dvector *fitness, size_t nselect, size_t init, uivector *selection)
{
  size_t i, k;
  double sum, ptr;
  k = 0;
  sum = 0.f;
  srand_(init);

  for(k = 0; k < nselect; k++){
    ptr = (double)rand_()/(double)RAND_MAX;
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

void Combinations(uivector *num, matrix *comb)
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
    interpolate(xy, intervals, interp_xy);
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

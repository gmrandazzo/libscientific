/* numeric.h
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

#ifndef NUMERIC_H
#define NUMERIC_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "vector.h"
#include "matrix.h"

/* FLOATING POINT NUMBER PRECISION */
#define EPSILON 1e-3  /* Define your own tolerance*/
#define FLOAT_EQ(x,v, EPSILON) (((v - EPSILON) < x) && (x <( v + EPSILON)))

#define _isnan_(a) (a != a)
#define _isinf_(a) (!_isnan_(a) && _isnan_(a - a))
#define _pi_ 3.14159265358979323846264338327950288419716939937510

#define MISSING 99999999

/* use this like:
 *    NumOne and NumTwo are the number to compare
 *         if(FLOAT_EQ(NumOne, NumTwo, EPSILON))...
 * END FLOATING POINT NUMBER PRECISION */

double missing_value();

// #ifdef WIN32
int myrand_r(unsigned int *seed);
// #endif
size_t Factorial(size_t x);
int randInt(int low, int high);
double randDouble(double low, double high);
double square(double x);
/*
 *Stocastic Universal Sample
 *
 * Input:
 * fitness = vector of weight to select...
 * nselect = number of objects to select
 * init = initializator for randomness...
 *
 * Output:
 * selection = vector of selected id
 *
 */
void StochasticUniversalSample(dvector *fitness, size_t nselect, size_t init, uivector *selection);

void RouletteWheelselection(dvector *fitness, size_t nselect, size_t init, uivector *selection);

void Combinations(uivector *num, matrix *comb);

/* Description: calculate area of a curve.
 * If intervals > 0 will interpolate the curve
 */
double curve_area(matrix *xy, size_t intervals);

#endif

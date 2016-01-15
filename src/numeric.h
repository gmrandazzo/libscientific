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
#define pi 3.1415926535897932384626433

/* use this like:
 *    NumOne and NumTwo are the number to compare
 *         if(FLOAT_EQ(NumOne, NumTwo, EPSILON))...
 * END FLOATING POINT NUMBER PRECISION */

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
void StochasticUniversalSample(dvector *fitness, size_t nselect, size_t init, uivector **selection);

void RouletteWheelselection(dvector *fitness, size_t nselect, size_t init, uivector **selection);

#endif

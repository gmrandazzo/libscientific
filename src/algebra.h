#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <stdio.h>
#include <stdlib.h>

#include "vector.h"
#include "matrix.h"

void SolveLSE(matrix *X, dvector **solution);
void OrdinaryLeastSquares(matrix *x, dvector *y, dvector *coefficients);

#endif

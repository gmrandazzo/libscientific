/* algebra.h
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

#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <stdio.h>
#include <stdlib.h>

#include "vector.h"
#include "matrix.h"

void SolveLSE(matrix *mx, dvector *solution);
void OrdinaryLeastSquares(matrix *x, dvector *y, dvector *coefficients);

#endif

/* statistic.h
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

#ifndef STATISTIC_H
#define STATISTIC_H

#include "matrix.h"
#include "vector.h"

/*
 * TP = True PositivePredictedValue
 * FN = False Negative
 * 
 * Sensitivity = TP / TP + FN
 */
void Sensitivity(dvector *tp, double thmin, double thmax, double thstep, matrix **s);

/*
 * TP = True Positive
 * FP = False Positive
 * 
 * PPV = TP / TP + FP
 */
void PositivePredictedValue(dvector *dtp, dvector *dtn, double thmin, double thmax, double thstep, matrix **p);


/*
 * This function code the matrix to a sign matrix 
 */
void MatrixCode(matrix *inmx, matrix *outmx);

void BifactorialMatrixExpansion(matrix* inmx, matrix* outmx);

/*
 * This function study the variable effect through the yates algorithm
 */
void YatesVarEffect(matrix *mx, dvector *veff);

#endif

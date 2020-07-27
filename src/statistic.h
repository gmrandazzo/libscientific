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
 * ytrue = True values
 * ypred = Predicted values
 *
 */
double R2(dvector *ytrue, dvector *ypred);

/*
 * ytrue = True values
 * ypred = Predicted values
 *
 */
double MAE(dvector *ytrue, dvector *ypred);

/*
 * ytrue = True values
 * ypred = Predicted values
 *
 */
double MSE(dvector *ytrue, dvector *ypred);

/*
 * ytrue = True values
 * ypred = Predicted values
 *
 */
double RMSE(dvector *ytrue, dvector *ypred);

/*
 * ytrue = True values
 * ypred = Predicted values
 *
 */
double BIAS(dvector *ytrue, dvector *ypred);

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

/* Description: calculate the ROC curve giving an y_true and an y_score */
void ROC(dvector *y_true, dvector *y_score, matrix **roc, double *auc);

/* Description: calculate the Precision-Recall curve giving an y_true and an y_score */
void PrecisionRecall(dvector *y_true, dvector *y_score,  matrix **pr, double *ap);

#endif

/* mlr.h
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

#include "matrix.h"
#include "scientificinfo.h"

/* mlr.h
 * Multiple Linear Regression
 * - Model Building
 * - Model Validation
 */

typedef struct{
  matrix *b;
  matrix *recalculated_y;
  matrix *predicted_y;
  matrix *recalc_residuals;
  matrix *pred_residuals;
  dvector *ymean;
  dvector *r2y_model;
  dvector *q2y;
  dvector *sdep;/* Standard Deviation over Prediction */
  dvector *sdec; /* Standard Deviation over Recalculating */
  dvector *bias; /* model bias */
  matrix *r2q2scrambling;
} MLRMODEL;

void NewMLRModel(MLRMODEL **m);
void DelMLRModel(MLRMODEL **m);

void MLR(matrix *mx, matrix *my, MLRMODEL *model, ssignal *s);

void MLRPredictY(matrix* mx, matrix *my, MLRMODEL* model, matrix** predicted_y, matrix** predicted_residuals, dvector** r2y, dvector** sdep);

/*
 * Description: Calculate the correlation coefficient (ccoeff),
 *              the standard deviation of the prediction (stdev),
 *              the bias of the prediction (bias) in a regression model.
 *              mx and my could be the training or the test datasets.
 */
void MLRRegressionStatistics(matrix *my_true, matrix *my_pred, dvector** ccoeff, dvector **stdev, dvector **bias);

/*Description:
 * Generate Random Models by putting randomly y and check if models have correlations...
 */
void MLRYScrambling(matrix *mx, matrix *my,
                        size_t block, size_t valtype, size_t rgcv_group, size_t rgcv_iterations,
                        matrix **r2q2scrambling, ssignal *s);

void MLRRandomGroupsCV(matrix *mx, matrix *my,
                        size_t group, size_t iterations,
                        dvector **q2y, dvector **sdep, dvector **bias, matrix **predicted_y, matrix** pred_residuals, ssignal *s); /* *cv_ are the cross validated coefficients used to plot predicted vs experimental. if is null is not calculated.*/

void MLRLOOCV(matrix *mx, matrix *my,
                        dvector **q2y, dvector **sdep, dvector **bias, matrix **predicted_y, matrix **pred_residuals, ssignal *s); /* *loov_ are the loo validated coefficients used to plot predicted vs experimental. if is null is not calculated.*/


void PrintMLR(MLRMODEL *m);

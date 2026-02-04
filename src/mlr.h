/* Implements Multiple Linear Regression (MLR).
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

#ifndef MLR_H
#define MLR_H

#include "matrix.h"
#include "scientificinfo.h"

/**
 * Multiple Linear Regression data structure
 * 
 * - **b** beta coefficients
 * - **recalculated_y** y recalculated
 * - **predicted_y** y predicted
 * - **recalc_residuals** residuals using recalculated values
 * - **pred_residuals** residuals using predicted values
 * - **ymean** y averages
 * - **r2y_model** model r-squared using y recalculated 
 * - **q2y** model q-squared using y predicted
 * - **sdep** standard deviation over prediction
 * - **sdec** standard deviation over recalculation
 * - **bias** model bias
 * - **r2q2scrambling** y scrambling results
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

/**
 * Initialize an MLR model
*/
void NewMLRModel(MLRMODEL **m);

/**
 * Delete a MLR model 
 */
void DelMLRModel(MLRMODEL **m);

/**
 * Calculate a MLR model using OLS algorithm
 * 
 * @param [in] mx x independent variables input
 * @param [in] my y dependent variables input 
 * @param [out] model initialized model using NewMLRModel(...). The datastructure will be populated with results
 * 
 */
void MLR(matrix *mx,
         matrix *my,
         MLRMODEL *model,
         ssignal *s);

/**
 * Predict Y targets using a MLR model
 * 
 * @param [in] mx x independent variables input
 * @param [in] my y dependent variables input (if known)
 * @param [in] model computed MLRMODEL
 * @param [out] predicted_y y predicted
 * @param [out] predicted_residuals residuals using y predicted
 * @param [out] r2y r-squared of the prediction if my known
 * @param [out] sdep standard deviation over prediction if my known
 * 
 */
void MLRPredictY(matrix* mx,
                 matrix *my,
                 MLRMODEL *model,
                 matrix *predicted_y,
                 matrix *predicted_residuals,
                 dvector *r2y,
                 dvector *sdep);

/**
 * Calculate the correlation coefficient (ccoeff),
 * the standard deviation of the prediction (stdev),
 * the bias of the prediction (bias) in a regression model.
 * mx and my could be the training or the test datasets.
 * 
 * @param [in] my_true input matrix of y true values
 * @param [in] my_pred input matrix of y predicted or recalculated values
 * @param [out] ccoeff output - correlation coefficients q-squared or r-squared if my_pred is respectivelly predicted or recalculated 
 * @param [out] rmse output - root mean square error
 * @param [out] bias output - bias defined as how distant are we from the diagonal
 */
 void MLRRegressionStatistics(matrix *my_true,
                              matrix *my_pred,
                              dvector *ccoeff,
                              dvector *rmse,
                              dvector *bias);

/**
 * @brief Print MRLMODEL to video.
 *
 * @param [in] m computed mlr model
 *
 * @par Returns
 *    Nothing.
 */
void PrintMLR(MLRMODEL *m);

#endif

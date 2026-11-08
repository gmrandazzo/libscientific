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
 * Multiple Linear Regression model (ordinary least squares).
 *
 * Create it with :c:func:`NewMLRModel`, fill it with :c:func:`MLR` and free
 * it with :c:func:`DelMLRModel`. With ``n`` objects, ``f`` features and
 * ``ny`` responses:
 *
 * Members:
 *   b: Coefficients, ``(f+1) x ny``: row 0 is the intercept, row ``j+1`` the
 *     coefficient of feature ``j``; column ``k`` refers to response ``k``.
 *   recalculated_y: Fitted responses of the training objects, ``n x ny``.
 *   predicted_y: Cross-validated predictions (filled by the validation
 *     functions of ``modelvalidation.h``), ``n x ny``.
 *   recalc_residuals: Fitted minus true responses, ``n x ny``.
 *   pred_residuals: Cross-validated residuals, ``n x ny``.
 *   ymean: Mean of each training response (size ``ny``).
 *   r2y_model: R² of the training fit, one value per response.
 *   q2y: Cross-validated R² (Q²), one value per response.
 *   sdep: Standard deviation of the prediction errors, per response.
 *   sdec: Standard deviation of the training errors, per response.
 *   bias: Bias of the predictions, per response.
 *   r2q2scrambling: Y-scrambling results.
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
 * Allocate an empty MLR model.
 *
 * Args:
 *   m: [out] Address of the pointer that receives the new model.
 *
 * Note:
 *   Release the model with :c:func:`DelMLRModel`.
 */
void NewMLRModel(MLRMODEL **m);

/**
 * Free an MLR model and all its members, and set the pointer to ``NULL``.
 *
 * Args:
 *   m: [in,out] Address of the model pointer created with :c:func:`NewMLRModel`.
 */
void DelMLRModel(MLRMODEL **m);

/**
 * Fit a multiple linear regression model by ordinary least squares.
 *
 * Every response is fitted independently with an intercept:
 * ``y_k = b[0][k] + sum_j b[j+1][k] x_j``. The training fit (fitted values,
 * residuals, R², SDEC) is stored in the model. Fitting an already fitted
 * model replaces its results.
 *
 * Args:
 *   mx: [in] Predictors, ``n x f`` (``n > f``). NaN/inf values are replaced
 *     by ``MISSING``.
 *   my: [in] Responses, ``n x ny``.
 *   model: [out] Model created with :c:func:`NewMLRModel`; filled with the
 *     results (see :c:type:`MLRMODEL`).
 *   s: [in] Unused, pass ``NULL``.
 *
 * Returns:
 *   Nothing, the results are written into ``model``.
 *
 * Example:
 *   .. code-block:: c
 *
 *      MLRMODEL *model;
 *      matrix *ypred;
 *      NewMLRModel(&model);
 *      MLR(x, y, model, NULL);
 *      // model->b->data[0][k]: intercept of response k
 *      initMatrix(&ypred);
 *      MLRPredictY(xnew, NULL, model, ypred, NULL, NULL, NULL);
 *      DelMatrix(&ypred);
 *      DelMLRModel(&model);
 */
void MLR(matrix *mx,
         matrix *my,
         MLRMODEL *model,
         scisignal *s);

/**
 * Predict the responses of new objects and, if their true values are known,
 * compute the prediction statistics.
 *
 * Args:
 *   mx: [in] Objects, ``m x f`` (same features as the training data).
 *   my: [in] True responses ``m x ny``, or ``NULL`` if unknown.
 *   model: [in] Model fitted with :c:func:`MLR`.
 *   predicted_y: [out] Matrix created with ``initMatrix``; resized to
 *     ``m x ny`` and filled with the predictions.
 *   predicted_residuals: [out] Predicted minus true responses, ``m x ny``.
 *     Only computed when ``my`` is given; ``NULL`` to skip.
 *   r2y: [out] Empty vector created with ``initDVector``: one R² value per
 *     response is **appended** (relative to the training mean). Only when
 *     ``my`` is given; ``NULL`` to skip.
 *   sdep: [out] Empty vector: one standard deviation of the errors per
 *     response is appended. Only when ``my`` is given; ``NULL`` to skip.
 *
 * Returns:
 *   Nothing, the results are written into the outputs. If the number of
 *   features differs from the model, an error is printed and nothing is
 *   computed.
 */
void MLRPredictY(matrix* mx,
                 matrix *my,
                 MLRMODEL *model,
                 matrix *predicted_y,
                 matrix *predicted_residuals,
                 dvector *r2y,
                 dvector *sdep);

/**
 * Compute R², RMSE and bias of regression predictions.
 *
 * Args:
 *   my_true: [in] True responses, ``n x ny``.
 *   my_pred: [in] Predicted (or fitted) responses, ``n x ny``.
 *   ccoeff: [out] Vector created with ``initDVector``; resized to ``ny``:
 *     R² of each response (Q² for cross-validated predictions).
 *     ``NULL`` to skip.
 *   rmse: [out] Root mean squared error of each response. ``NULL`` to skip.
 *   bias: [out] Bias of each response. ``NULL`` to skip.
 *
 * Returns:
 *   Nothing, the statistics are written into the output vectors.
 */
 void MLRRegressionStatistics(matrix *my_true,
                              matrix *my_pred,
                              dvector *ccoeff,
                              dvector *rmse,
                              dvector *bias);

/**
 * Print the content of an MLR model to standard output.
 *
 * Args:
 *   m: [in] The model to print.
 */
void PrintMLR(MLRMODEL *m);

#endif

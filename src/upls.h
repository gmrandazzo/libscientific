/* Implements UPLS algorithms.
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

#ifndef UPLS_H
#define UPLS_H
#include "tensor.h"
#include "matrix.h"
#include "vector.h"
#include "preprocessing.h"
#include "scientificinfo.h"

#define UPLSCONVERGENCE 1e-3

/**
 * Multi-way (unfolded) Partial Least Squares regression model.
 *
 * X and Y are tensors of layers with the same objects: ``X->m[k]`` is
 * ``n x f`` (``Kx`` layers), ``Y->m[k]`` is ``n x ny`` (``Ky`` layers; a single
 * response matrix is a tensor with one layer). Create the model with
 * :c:func:`NewUPLSModel`, fill it with :c:func:`UPLS` and free it with
 * :c:func:`DelUPLSModel`. With ``nlv`` latent variables:
 *
 * Members:
 *   xscores: X scores, ``n x nlv``.
 *   xloadings: One ``f x Kx`` X loading matrix per latent variable.
 *   xweights: One ``f x Kx`` X weight matrix per latent variable.
 *   xvarexp: Percentage of the X variance explained by each latent variable.
 *   yscores: Y scores, ``n x nlv``.
 *   yloadings: One ``ny x Ky`` Y loading matrix per latent variable.
 *   b: Inner regression coefficient of each latent variable (1 for a single
 *     response: the coefficient is carried by the Y loadings).
 *   xcolaverage: Column means of each X layer.
 *   xcolscaling: Column scalings of each X layer.
 *   ycolaverage: Column means of each Y layer.
 *   ycolscaling: Column scalings of each Y layer.
 *   r2x_model: R² of X for 1..nlv latent variables (:c:func:`UPLSRSquared`).
 *   r2x_validation: Cross-validated R² of X.
 *   r2y_model: R² of Y: one ``ny x Ky`` matrix per latent variable.
 *   r2y_validation: Cross-validated R² of Y, same layout.
 *   q2y: Q² of Y, same layout (validation functions).
 *   sdep: Standard deviation of the prediction errors, same layout.
 *   sdec: Standard deviation of the recalculation errors, same layout.
 *   recalculated_y: Fitted Y of the training objects: one matrix per Y
 *     layer, ``n x (ny*nlv)``.
 *   predicted_y: Cross-validated Y, same layout.
 *   recalc_residuals: Fitted minus true Y, same layout.
 *   pred_residuals: Cross-validated residuals, same layout.
 *   q2y_yscrambling: Q² of the y-scrambling models.
 *   sdep_yscrambling: SDEP of the y-scrambling models.
 */
typedef struct{
  matrix *xscores;
  tensor *xloadings;
  tensor *xweights;
  dvector *xvarexp;
  matrix *yscores;
  tensor *yloadings;
  dvector *b;
  dvectorlist *xcolaverage;
  dvectorlist *xcolscaling;
  dvectorlist *ycolaverage;
  dvectorlist *ycolscaling;
  dvector *r2x_model;
  dvector *r2x_validation;
  tensor *r2y_model; /* each order correspond to a principal component; each row correspond to an y dependent variable and each col to an order */
  tensor *r2y_validation;  /* each order correspond to a principal component; each row correspond to an y dependent variable and each col to an order */
  tensor *q2y;
  tensor *sdep; /* Standard Deviation over Prediction */
  tensor *sdec; /* Standard Deviation over Recalculating */
  tensor *recalculated_y;
  tensor *predicted_y;
  tensor *recalc_residuals;
  tensor *pred_residuals;
  tensor *q2y_yscrambling;
  tensor *sdep_yscrambling;
} UPLSMODEL;

/**
 * Allocate an empty UPLS model.
 *
 * Args:
 *   m: [out] Address of the pointer that receives the new model.
 *
 * Note:
 *   Release the model with :c:func:`DelUPLSModel`.
 */
void NewUPLSModel(UPLSMODEL **m);

/**
 * Free a UPLS model and all its members, and set the pointer to ``NULL``.
 *
 * Args:
 *   m: [in,out] Address of the model pointer created with :c:func:`NewUPLSModel`.
 */
void DelUPLSModel(UPLSMODEL **m);

/**
 * Fit a multi-way PLS regression model (Wold, Geladi and Esbensen,
 * J. Chemometrics 1, 41-56, 1987).
 *
 * Every column of every X and Y layer is preprocessed, then the latent
 * variables are extracted one at a time with NIPALS (up to 1000 iterations,
 * relative tolerance ``UPLSCONVERGENCE`` = 1e-3). When Y is fully explained
 * the remaining latent variables are skipped (a message is printed on
 * standard error), so the model can have fewer than ``npc`` latent
 * variables. Fitting an already fitted model replaces its results.
 *
 * Args:
 *   X_: [in] X layers, all ``n x f``.
 *   Y_: [in] Y layers, all ``n x ny``.
 *   npc: [in] Number of latent variables; capped to ``f``.
 *   xautoscaling: [in] Column preprocessing of X: 0 = mean centering only,
 *     1 = standard deviation, 2 = root mean square, 3 = Pareto, 4 = range,
 *     5 = level.
 *   yautoscaling: [in] Column preprocessing of Y, same options.
 *   m: [out] Model created with :c:func:`NewUPLSModel`; filled with the
 *     results (see :c:type:`UPLSMODEL`).
 *   s: [in] Unused, pass ``NULL``.
 *
 * Returns:
 *   Nothing, the results are written into ``m``.
 */
void UPLS(tensor *X_,
          tensor *Y_,
          size_t npc,
          size_t xautoscaling,
          size_t yautoscaling,
          UPLSMODEL *m,
          scisignal *s);

/**
 * Compute the X scores of new objects.
 *
 * Args:
 *   X_: [in] X layers of the new objects (same layers and features as the
 *     training X, ``m`` rows each).
 *   m: [in] Model fitted with :c:func:`UPLS`.
 *   npc: [in] Number of latent variables to use.
 *   ptscores: [out] Matrix created with ``initMatrix``; resized to
 *     ``m x npc``.
 *
 * Returns:
 *   Nothing, the scores are written into ``ptscores``.
 */
void UPLSScorePredictor(tensor *X_,
                        UPLSMODEL *m,
                        size_t npc,
                        matrix *ptscores);

/**
 * Compute the responses from X scores.
 *
 * Args:
 *   tscores: [in] X scores (from :c:func:`UPLSScorePredictor`).
 *   m: [in] Model fitted with :c:func:`UPLS`.
 *   npc: [in] Number of latent variables to use.
 *   py: [out] Tensor created with ``initTensor``; receives one ``m x ny``
 *     matrix per Y layer, in the original units of Y.
 *
 * Returns:
 *   Nothing, the responses are written into ``py``.
 */
void UPLSYPredictor(matrix *tscores,
                    UPLSMODEL *m,
                    size_t npc,
                    tensor *py);

/**
 * Compute R² of X and Y and the SDEC for 1..npc latent variables.
 *
 * Args:
 *   X_: [in] X layers (training or test objects).
 *   Y_: [in] True Y layers of the same objects.
 *   m: [in] Model fitted with :c:func:`UPLS`.
 *   npc: [in] Maximum number of latent variables.
 *   r2x: [out] Vector created with ``initDVector``; one cumulative R² of X
 *     per number of latent variables.
 *   r2y: [out] Tensor created with ``initTensor``; one ``ny x Ky`` matrix of
 *     R² values per number of latent variables.
 *   sdec: [out] Tensor; standard deviation of the errors, same layout.
 *
 * Returns:
 *   Nothing, the statistics are written into the outputs.
 */
void UPLSRSquared(tensor *X_,
                  tensor *Y_,
                  UPLSMODEL *m,
                  size_t npc,
                  dvector *r2x,
                  tensor *r2y,
                  tensor *sdec);

/**
 * Compute the residual and total sums of squares of X and Y for 1..npc
 * latent variables (used by the cross-validation functions).
 *
 * Args:
 *   X_: [in] X layers.
 *   Y_: [in] True Y layers.
 *   m: [in] Model fitted with :c:func:`UPLS`.
 *   npc: [in] Maximum number of latent variables.
 *   xss_err: [out] Residual sum of squares of X, one value per number of
 *     latent variables.
 *   xss_tot: [out] Total sum of squares of X, same size.
 *   yss_err: [out] Residual sums of squares of Y: one ``ny x Ky`` matrix per
 *     number of latent variables.
 *   yss_tot: [out] Total sums of squares of Y, same layout.
 *   predicted_y: [out] Predicted Y layers. ``NULL`` to skip.
 *
 * Returns:
 *   Nothing, the results are written into the outputs.
 */
void UPLSRSquared_SSErr_SStot(tensor *X_,
                              tensor *Y_,
                              UPLSMODEL *m,
                              size_t npc,
                              dvector *xss_err,
                              dvector *xss_tot,
                              tensor *yss_err,
                              tensor *yss_tot,
                              tensor *predicted_y);

/**
 * Y-scrambling test: refit models on permuted Y to detect chance
 * correlation (overfitting).
 *
 * Args:
 *   X_: [in] X layers.
 *   Y_: [in] Y layers.
 *   xautoscaling: [in] X preprocessing (see :c:func:`UPLS`).
 *   yautoscaling: [in] Y preprocessing.
 *   npc: [in] Number of latent variables.
 *   block: [in] Number of scrambling blocks.
 *   valtype: [in] Validation of each scrambled model: 0 = leave-one-out,
 *     otherwise random groups.
 *   rgcv_groups: [in] Number of groups for random-group validation.
 *   rgcv_iterations: [in] Number of iterations for random-group validation.
 *   q2y: [out] Address of a tensor receiving the Q² of the scrambled models.
 *   sdep: [out] Address of a tensor receiving their SDEP.
 *   s: [in] Optional stop signal; ``NULL`` if unused.
 *
 * Returns:
 *   Nothing, the results are written into ``q2y`` and ``sdep``.
 */
void UPLSYScrambling(tensor *X_, tensor *Y_,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t npc, size_t block,
                        size_t valtype, size_t rgcv_groups, size_t rgcv_iterations,
                        tensor **q2y, tensor **sdep, scisignal *s);

/**
 * Estimate the predictive ability of a UPLS model by random-group
 * cross-validation.
 *
 * Args:
 *   X_: [in] X layers.
 *   Y_: [in] Y layers.
 *   xautoscaling: [in] X preprocessing (see :c:func:`UPLS`).
 *   yautoscaling: [in] Y preprocessing.
 *   npc: [in] Maximum number of latent variables.
 *   group: [in] Number of cross-validation groups.
 *   iterations: [in] Number of random repetitions.
 *   r2x: [out] Cross-validated R² of X for 1..npc latent variables.
 *   q2y: [out] Q² of Y: one ``ny x Ky`` matrix per number of latent variables.
 *   sdep: [out] Standard deviation of the prediction errors, same layout.
 *   predicted_y: [out] Cross-validated predictions.
 *   pred_residuals: [out] Cross-validated residuals.
 *   s: [in] Optional stop signal; ``NULL`` if unused.
 *
 * Returns:
 *   Nothing, the results are written into the outputs.
 */
void UPLSRandomGroupsCV(tensor *X_, tensor *Y_,
                         size_t xautoscaling, size_t yautoscaling,
                         size_t npc, size_t group, size_t iterations,
                         dvector **r2x, tensor **q2y, tensor **sdep, tensor **predicted_y, tensor **pred_residuals, scisignal *s);

/**
 * Estimate the predictive ability of a UPLS model by leave-one-out
 * cross-validation.
 *
 * Args:
 *   X_: [in] X layers.
 *   Y_: [in] Y layers.
 *   xautoscaling: [in] X preprocessing (see :c:func:`UPLS`).
 *   yautoscaling: [in] Y preprocessing.
 *   npc: [in] Maximum number of latent variables.
 *   r2x: [out] Cross-validated R² of X for 1..npc latent variables.
 *   q2y: [out] Q² of Y: one ``ny x Ky`` matrix per number of latent variables.
 *   sdep: [out] Standard deviation of the prediction errors, same layout.
 *   predicted_y: [out] Cross-validated predictions.
 *   pred_residuals: [out] Cross-validated residuals.
 *   s: [in] Optional stop signal; ``NULL`` if unused.
 *
 * Returns:
 *   Nothing, the results are written into the outputs.
 */
void UPLSLOOCV(tensor *X_, tensor *Y_,
                         size_t xautoscaling, size_t yautoscaling,
                         size_t npc,
                         dvector **r2x, tensor **q2y, tensor **sdep, tensor **predicted_y, tensor **pred_residuals, scisignal *s);

/**
 * Choose the number of latent variables from a table of R² or Q² values.
 *
 * Looks for the best value and stops when it drops.
 *
 * Args:
 *   rq2: [in] R² or Q² values: one matrix per number of latent variables
 *     (e.g. ``model->q2y``).
 *
 * Returns:
 *   The **0-based** index of the selected latent variable: the number of
 *   latent variables to use is the return value + 1.
 */
size_t UPLSGetPCModelCutOff(tensor *rq2);

/**
 * Print the content of a UPLS model to standard output.
 *
 * Args:
 *   m: [in] The model to print.
 */
void PrintUPLSModel(UPLSMODEL *m);

#endif

/* Implements Partial Least Squares (PLS) regression.
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

#ifndef PLS_H
#define PLS_H
#include "tensor.h"
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define PLSCONVERGENCE 1e-8

/**
 * Partial Least Squares regression model.
 *
 * Create it with :c:func:`NewPLSModel`, fill it with :c:func:`PLS` and free
 * it with :c:func:`DelPLSModel`. All the members are owned by the model.
 * With ``n`` objects, ``f`` features, ``ny`` responses and ``nlv`` latent
 * variables:
 *
 * Members:
 *   xscores: X scores ``T``, ``n x nlv``.
 *   xloadings: X loadings ``P``, ``f x nlv``.
 *   xweights: X weights ``W``, ``f x nlv`` (used to compute the scores of
 *     new objects; not unit length).
 *   yscores: Y scores ``U``, ``n x nlv``.
 *   yloadings: Y loadings ``Q``, ``ny x nlv`` (all ``1`` for one response).
 *   b: Inner regression coefficient of each latent variable (size ``nlv``).
 *   xvarexp: Percentage of the X variance explained by each latent variable.
 *   xcolaverage: Column means of X (size ``f``).
 *   xcolscaling: Column scalings of X (size ``f``).
 *   ycolaverage: Column means of Y (size ``ny``).
 *   ycolscaling: Column scalings of Y (size ``ny``).
 *   recalculated_y: Fitted Y of the training objects, ``n x (ny*nlv)``:
 *     column ``ny*lv + j`` is response ``j`` with ``lv+1`` latent variables.
 *   recalc_residuals: ``recalculated_y`` minus the true Y, same layout.
 *   predicted_y: Cross-validated predictions (filled by the validation
 *     functions of ``modelvalidation.h``), same layout.
 *   pred_residuals: Cross-validated residuals, same layout.
 *   r2y_recalculated: R² of the fitted Y, ``nlv x ny``.
 *   r2y_validation: R² of the cross-validated Y, ``nlv x ny``.
 *   q2y: Q² (cross-validated R²), ``nlv x ny``.
 *   sdep: Standard deviation of the prediction errors, ``nlv x ny``.
 *   sdec: Standard deviation of the recalculation errors, ``nlv x ny``.
 *   bias: Prediction bias, ``nlv x ny``.
 *   roc_recalculated: PLS-DA ROC curves of the fitted Y (one matrix per
 *     latent variable).
 *   roc_validation: PLS-DA ROC curves of the cross-validated Y.
 *   roc_auc_recalculated: ROC AUC of the fitted Y, ``nlv x ny``.
 *   roc_auc_validation: ROC AUC of the cross-validated Y, ``nlv x ny``.
 *   precision_recall_recalculated: Precision-recall curves of the fitted Y.
 *   precision_recall_validation: Precision-recall curves of the
 *     cross-validated Y.
 *   precision_recall_ap_recalculated: Average precision of the fitted Y.
 *   precision_recall_ap_validation: Average precision of the
 *     cross-validated Y.
 *   yscrambling: Y-scrambling results (R² and Q² of models fitted on
 *     permuted Y).
 */
typedef struct{
  matrix *xscores;
  matrix *xloadings;
  matrix *xweights;
  matrix *yscores;
  matrix *yloadings;
  dvector *b;
  dvector *xvarexp;
  dvector *xcolaverage;
  dvector *xcolscaling;
  dvector *ycolaverage;
  dvector *ycolscaling;
  matrix *recalculated_y;
  matrix *recalc_residuals;
  matrix *predicted_y;
  matrix *pred_residuals;
  /* Regression variables */
  matrix *r2y_recalculated; /* each column correspond to an y dependent variable and each row correspond to a principal component*/
  matrix *r2y_validation;
  matrix *q2y;
  matrix *sdep; /* Standard Deviation over Prediction */
  matrix *sdec; /* Standard Deviation over Recalculating */
  matrix *bias;
  /* Discriminant Analyisis variables */
  tensor *roc_recalculated;
  tensor *roc_validation;
  matrix *roc_auc_recalculated;
  matrix *roc_auc_validation;
  tensor *precision_recall_recalculated;
  tensor *precision_recall_validation;
  matrix *precision_recall_ap_recalculated;
  matrix *precision_recall_ap_validation;
  matrix *yscrambling;
} PLSMODEL;

/**
 * Allocate an empty PLS model.
 *
 * Args:
 *   m: [out] Address of the pointer that receives the new model.
 *
 * Note:
 *   Release the model with :c:func:`DelPLSModel`.
 */
void NewPLSModel(PLSMODEL **m);

/**
 * Free a PLS model and all its members, and set the pointer to ``NULL``.
 *
 * Args:
 *   m: [in,out] Address of the model pointer created with :c:func:`NewPLSModel`.
 */
void DelPLSModel(PLSMODEL **m);

/*
 * Internal helper: compute one latent variable with NIPALS (Geladi &
 * Kowalski, Anal. Chim. Acta 185, 1986), returning the scores t, u, the
 * loadings p, q, the weights w and the inner coefficient bcoef.
 * All the vectors must be allocated. Used by PLS().
 */
void LVCalc(matrix *X, matrix *Y, dvector *t, dvector *u, dvector *p, dvector *q, dvector *w, double *bcoef);


/**
 * Fit a PLS regression model with the NIPALS algorithm.
 *
 * X and Y are preprocessed, then ``nlv`` latent variables are extracted one
 * at a time (up to 1000 iterations each). The model also stores the fitted Y
 * of the training objects for every number of latent variables.
 *
 * Args:
 *   mx: [in] Predictors, ``n x f``.
 *   my: [in] Responses, ``n x ny`` (one column per response).
 *   nlv: [in] Number of latent variables; capped to ``f``.
 *   xautoscaling: [in] Column preprocessing of X: 0 = mean centering only,
 *     1 = standard deviation, 2 = root mean square, 3 = Pareto, 4 = range,
 *     5 = level (all options mean-center first).
 *   yautoscaling: [in] Column preprocessing of Y, same options.
 *   model: [out] Model created with :c:func:`NewPLSModel`; filled with the
 *     results (see :c:type:`PLSMODEL`).
 *   s: [in] Optional stop signal (``*s == SIGSCIENTIFICSTOP`` stops the
 *     computation); ``NULL`` if unused.
 *
 * Returns:
 *   Nothing, the results are written into ``model``.
 *
 * Example:
 *   .. code-block:: c
 *
 *      PLSMODEL *model;
 *      matrix *ypred, *tscores;
 *      NewPLSModel(&model);
 *      PLS(x, y, 3, 1, 0, model, NULL);   // 3 LVs, autoscale X, center Y
 *      initMatrix(&tscores);
 *      initMatrix(&ypred);
 *      PLSYPredictorAllLV(xnew, model, tscores, ypred);
 *      // ypred->data[i][ny*lv + j]: response j of object i with lv+1 LVs
 *      DelMatrix(&ypred);
 *      DelMatrix(&tscores);
 *      DelPLSModel(&model);
 */
void PLS(matrix *mx,
         matrix *my,
         size_t nlv,
         int xautoscaling,
         int yautoscaling,
         PLSMODEL *model,
         scisignal *s);

/**
 * Compute the regression coefficients ``B = W (P' W)^-1 b`` with ``nlv``
 * latent variables.
 *
 * They express the model as a linear equation on the preprocessed data:
 * ``(y - ycolaverage) / ycolscaling = sum_j B[j] (x[j] - xcolaverage[j]) / xcolscaling[j]``.
 *
 * Args:
 *   model: [in] Model fitted with :c:func:`PLS` on a **single** response.
 *   nlv: [in] Number of latent variables to use.
 *   betas: [out] Vector created with ``initDVector``; filled with ``f``
 *     coefficients.
 *
 * Returns:
 *   Nothing, the coefficients are written into ``betas``.
 *
 * Warning:
 *   Only valid for models with one response: the Y loadings are not used,
 *   so the result is meaningless when ``ny > 1``.
 */
void PLSBetasCoeff(PLSMODEL *model,
                   size_t nlv,
                   dvector *betas);

/**
 * Compute the X scores of new objects.
 *
 * Args:
 *   mx: [in] Objects, ``m x f``, same features as the training X.
 *   model: [in] Model fitted with :c:func:`PLS`.
 *   nlv: [in] Number of latent variables to use.
 *   xscores: [out] Matrix created with ``initMatrix``; resized to ``m x nlv``.
 *
 * Returns:
 *   Nothing, the scores are written into ``xscores``.
 */
void PLSScorePredictor(matrix *mx,
                       PLSMODEL *model,
                       size_t nlv,
                       matrix *xscores);

/**
 * Compute the responses from X scores with ``nlv`` latent variables.
 *
 * Args:
 *   tscore: [in] X scores, ``m x a`` with ``a >= nlv`` (from
 *     :c:func:`PLSScorePredictor` or ``model->xscores``).
 *   model: [in] Model fitted with :c:func:`PLS`.
 *   nlv: [in] Number of latent variables to use.
 *   y: [out] Matrix created with ``initMatrix``; resized to ``m x ny`` and
 *     filled with the responses in the original units of Y.
 *
 * Returns:
 *   Nothing, the responses are written into ``y``.
 */
void PLSYPredictor(matrix *tscore,
                   PLSMODEL *model,
                   size_t nlv,
                   matrix *y);

/**
 * Predict the responses of new objects for every number of latent variables.
 *
 * Args:
 *   mx: [in] Objects, ``m x f``, same features as the training X.
 *   model: [in] Model fitted with :c:func:`PLS`.
 *   tscores: [out] Matrix created with ``initMatrix``; receives the X
 *     scores, ``m x nlv``. Pass ``NULL`` if not needed.
 *   y: [out] Matrix created with ``initMatrix``; resized to
 *     ``m x (ny*nlv)``: column ``ny*lv + j`` is response ``j`` with ``lv+1``
 *     latent variables.
 *
 * Returns:
 *   Nothing, the results are written into ``tscores`` and ``y``.
 */
void PLSYPredictorAllLV(matrix *mx,
                        PLSMODEL *model,
                        matrix *tscores,
                        matrix *y);


/**
 * Compute R², RMSE and bias of regression predictions for every number of
 * latent variables.
 *
 * Args:
 *   my_true: [in] True responses, ``n x ny``.
 *   my_pred: [in] Predictions in the layout of :c:func:`PLSYPredictorAllLV`,
 *     ``n x (ny*nlv)`` (fitted values give R², cross-validated ones Q²).
 *   ccoeff: [out] Matrix created with ``initMatrix``; resized to
 *     ``nlv x ny``: element ``[lv][j]`` is R² of response ``j`` with
 *     ``lv+1`` latent variables. ``NULL`` to skip.
 *   rmse: [out] Root mean squared error, ``nlv x ny``. ``NULL`` to skip.
 *   bias: [out] Bias of the predictions, ``nlv x ny``. ``NULL`` to skip.
 *
 * Returns:
 *   Nothing, the statistics are written into the output matrices.
 */
void PLSRegressionStatistics(matrix *my_true,
                             matrix *my_pred,
                             matrix *ccoeff,
                             matrix *rmse,
                             matrix *bias);

/**
 * Compute ROC and precision-recall curves of a PLS discriminant analysis
 * (PLS-DA) for every number of latent variables.
 *
 * Args:
 *   my_true: [in] True classes, ``n x ny`` (one 0/1 column per class).
 *   my_score: [in] Predicted class scores in the layout of
 *     :c:func:`PLSYPredictorAllLV`, ``n x (ny*nlv)``.
 *   roc: [out] Tensor created with ``initTensor``; receives one ``n x 2ny``
 *     matrix per latent variable (columns ``2j``, ``2j+1``: false and true
 *     positive rates of class ``j``). ``NULL`` to skip.
 *   roc_auc: [out] Matrix created with ``initMatrix``; one row per latent
 *     variable, one column per class (area under the ROC curve).
 *     ``NULL`` to skip.
 *   precision_recall: [out] Tensor created with ``initTensor``; one
 *     ``n x 2ny`` matrix per latent variable (columns ``2j``, ``2j+1``:
 *     recall and precision of class ``j``). ``NULL`` to skip.
 *   precision_recall_ap: [out] Matrix; average precision, same layout as
 *     ``roc_auc``. ``NULL`` to skip.
 *
 * Returns:
 *   Nothing, the results are written into the outputs.
 */
void PLSDiscriminantAnalysisStatistics(matrix *my_true,
                                       matrix *my_score,
                                       tensor *roc,
                                       matrix *roc_auc,
                                       tensor *precision_recall,
                                       matrix *precision_recall_ap);

/**
 * Compute the Variable Importance in Projection (VIP) of each feature.
 *
 * Uses S. Wold's definition over all the latent variables of the model:
 * ``VIP_j = sqrt(f * sum_a (w_ja² / ||w_a||²) SSY_a / SSY)``. The average of
 * the squared VIP values is 1; features with ``VIP > 1`` are usually
 * considered important.
 *
 * Args:
 *   model: [in] Model fitted with :c:func:`PLS`.
 *   vip: [out] Matrix created with ``initMatrix``; resized to ``f x 1``.
 *
 * Returns:
 *   Nothing, the values are written into ``vip``.
 */
void PLSVIP(PLSMODEL *model, matrix *vip);

/**
 * Choose the number of latent variables from a table of R²/Q² or AUC values.
 *
 * The values of each row are summed over the responses; the function walks
 * the rows while the sum increases by more than 3% and stops at the first
 * smaller improvement.
 *
 * Args:
 *   coeff: [in] Statistics, ``nlv x ny`` (e.g. ``model->q2y`` or the
 *     ``ccoeff`` output of :c:func:`PLSRegressionStatistics`).
 *
 * Returns:
 *   The **0-based** index of the selected row: the number of latent
 *   variables to use is the return value + 1.
 */
int GetLVCCutoff(matrix* coeff);

/**
 * Print the content of a PLS model to standard output.
 *
 * Args:
 *   model: [in] The model to print.
 */
void PrintPLSModel(PLSMODEL *model);

#endif

/* Implements Principal Component Analysis (PCA).
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

#ifndef PCA_H
#define PCA_H
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define PCACONVERGENCE 1e-10

/**
 * Principal Component Analysis model.
 *
 * Create it with :c:func:`NewPCAModel`, fill it with :c:func:`PCA` and free
 * it with :c:func:`DelPCAModel`. All the members are allocated and owned by
 * the model: read them, but do not free or resize them.
 *
 * With ``n`` objects (rows), ``f`` features (columns) and ``npc`` components:
 *
 * Members:
 *   scores: Scores ``T``, ``n x npc``: coordinates of each object in the
 *     component space.
 *   loadings: Loadings ``P``, ``f x npc``: weight of each feature in each
 *     component. Each column has unit length.
 *   dmodx: Distance of each object to the model, ``n x npc``: element
 *     ``[i][k]`` is the norm of the residual of object ``i`` after removing
 *     components ``0..k``.
 *   varexp: Variance explained by each component, in percent of the total
 *     variance of the preprocessed matrix (size ``npc``, not cumulative).
 *   colaverage: Column means of the training matrix (size ``f``).
 *   colscaling: Column scaling factors (size ``f``; all ``1`` for scaling 0).
 */
typedef struct{
  matrix *scores;
  matrix *loadings;
  matrix *dmodx;
  dvector *varexp;
  dvector *colaverage;
  dvector *colscaling;
} PCAMODEL;

/**
 * Allocate an empty PCA model.
 *
 * Args:
 *   m: [out] Address of the pointer that receives the new model.
 *
 * Note:
 *   Release the model with :c:func:`DelPCAModel`.
 */
void NewPCAModel(PCAMODEL **m);

/**
 * Free a PCA model and all its members, and set the pointer to ``NULL``.
 *
 * Args:
 *   m: [in,out] Address of the model pointer created with :c:func:`NewPCAModel`.
 */
void DelPCAModel(PCAMODEL **m);

/*
 * Internal helper (used by PCA, PLS and CPCA): fill varexp with the
 * percentage of the total sum of squares ss explained by each eigenvalue.
 */
void calcVarExpressed(double ss, dvector *eval, dvector *varexp);

/*
 * Internal helper: object distance used inside the PCA algorithm.
 */
double calcObjectDistance(matrix *m);

/*
 * Internal helper: NIPALS convergence criterion
 * (DOI: 10.1002/cem.1180010107, page 51), used by PCA/CPCA/PLS.
 */
double calcConvergence(dvector *t_new, dvector *t_old);


/**
 * Fit a PCA model with the NIPALS algorithm.
 *
 * The columns of ``mx`` are preprocessed according to ``scaling``, then the
 * ``npc`` components are extracted one at a time. Each component gets up to
 * 1000 NIPALS iterations; a component that does not converge keeps its best
 * estimate. Missing values (``MISSING``, NaN, inf) are skipped.
 *
 * Args:
 *   mx: [in] Training matrix, ``n x f`` (objects x features). Not modified,
 *     except that NaN/inf values are replaced by ``MISSING``.
 *   scaling: [in] Column preprocessing. All the options mean-center the
 *     columns, then: 0 = no scaling, 1 = standard deviation (autoscaling),
 *     2 = root mean square, 3 = Pareto (square root of the standard
 *     deviation), 4 = range (max - min), 5 = level (column mean).
 *   npc: [in] Number of components to compute; capped to ``f``.
 *   model: [out] Model created with :c:func:`NewPCAModel`. Its members are
 *     resized and filled with the results (see :c:type:`PCAMODEL`).
 *   s: [in] Optional stop signal: when ``*s == SIGSCIENTIFICSTOP`` the
 *     computation stops after the current component. Pass ``NULL`` if unused.
 *
 * Returns:
 *   Nothing, the results are written into ``model``.
 *
 * Example:
 *   .. code-block:: c
 *
 *      matrix *x;
 *      PCAMODEL *model;
 *      NewMatrix(&x, 5, 3);          // 5 objects, 3 features
 *      // ... fill x->data[i][j] ...
 *      NewPCAModel(&model);
 *      PCA(x, 1, 2, model, NULL);    // autoscaling, 2 components
 *      // model->scores: 5 x 2, model->loadings: 3 x 2,
 *      // model->varexp->data[k]: % of variance of component k
 *      DelPCAModel(&model);
 *      DelMatrix(&x);
 */
void PCA(matrix *mx,
         int scaling,
         size_t npc,
         PCAMODEL *model,
         scisignal *s);


/**
 * Project new objects into a fitted PCA model.
 *
 * The objects are preprocessed with the column means and scalings of the
 * training data, then projected onto the loadings. Projecting the training
 * matrix gives the same values as ``model->scores``.
 *
 * Args:
 *   mx: [in] Objects to project, ``m x f``, with the same features (same
 *     order) as the training matrix.
 *   model: [in] Model fitted with :c:func:`PCA`.
 *   npc: [in] Number of components to use; capped to the model components.
 *   pscores: [out] Matrix created with ``initMatrix`` or ``NewMatrix``. It is
 *     resized and filled with the predicted scores, ``m x npc``. The caller
 *     frees it with ``DelMatrix``.
 *
 * Returns:
 *   Nothing, the scores are written into ``pscores``.
 *
 * Example:
 *   .. code-block:: c
 *
 *      matrix *pscores;
 *      initMatrix(&pscores);
 *      PCAScorePredictor(xnew, model, 2, pscores);   // xnew: m x f
 *      // pscores->data[i][k]: score of new object i on component k
 *      DelMatrix(&pscores);
 */
void PCAScorePredictor(matrix *mx,
                       PCAMODEL *model,
                       size_t npc,
                       matrix *pscores);

/**
 * Rebuild data in the original feature space from scores: ``X = T P'``.
 *
 * The product of the first ``npc`` scores and loadings is computed, then the
 * scaling and centering are reverted, so the result is in the units of the
 * training data. With all the components the training matrix is rebuilt
 * exactly.
 *
 * Args:
 *   t: [in] Scores, ``n x a`` with ``a >= npc`` (e.g. ``model->scores`` or the
 *     output of :c:func:`PCAScorePredictor`).
 *   p: [in] Loadings, ``f x a`` (``model->loadings``).
 *   colaverage: [in] Column means, size ``f`` (``model->colaverage``).
 *   colscaling: [in] Column scalings, size ``f`` (``model->colscaling``).
 *   npc: [in] Number of components to use; capped to the columns of ``t``.
 *   mx: [out] Matrix created with ``initMatrix``; resized and filled with the
 *     reconstructed data, ``n x f``. The caller frees it with ``DelMatrix``.
 *
 * Returns:
 *   Nothing, the result is written into ``mx``.
 */
void PCAIndVarPredictor(matrix *t,
                        matrix *p,
                        dvector *colaverage,
                        dvector *colscaling,
                        size_t npc,
                        matrix *mx);

/**
 * Compute how well a model reconstructs a matrix (R²) for 1..npc components.
 *
 * For each number of components ``a = 1..npc`` the objects of ``mx`` are
 * projected and reconstructed with ``a`` components, and
 * ``R²(a) = 1 - SS(residual) / SS(mx centered on the model column means)`` is
 * appended to ``r2``. The values are cumulative: ``r2[a-1]`` refers to the
 * first ``a`` components together.
 *
 * Args:
 *   mx: [in] Matrix to reconstruct, ``n x f`` (training data or new objects).
 *   model: [in] Model fitted with :c:func:`PCA`.
 *   npc: [in] Maximum number of components; capped to the model components.
 *   r2: [out] Empty vector created with ``initDVector``: ``npc`` values are
 *     appended to it. The caller frees it with ``DelDVector``.
 *
 * Returns:
 *   Nothing, the values are appended to ``r2``.
 */
void PCARSquared(matrix *mx,
                 PCAMODEL *model,
                 size_t npc,
                 dvector *r2);

/**
 * Estimate the number of significant components by cross-validation.
 *
 * The objects are split at random into ``group`` groups. Each group in turn
 * is left out, a PCA model is fitted on the other objects and
 * :c:func:`PCARSquared` is computed on the left-out group. The procedure is
 * repeated ``iterations`` times and the results are averaged, giving the
 * cross-validated R² (Q²) for 1..npc components (S. Wold, Technometrics 20,
 * 397-405, 1978). Choose the number of components where Q² stops increasing.
 * The random split is deterministic for a given input size.
 *
 * Args:
 *   mx: [in] Data matrix, ``n x f``.
 *   npc: [in] Maximum number of components to evaluate.
 *   scaling: [in] Column preprocessing, as in :c:func:`PCA`.
 *   group: [in] Number of cross-validation groups (e.g. 5 or 7).
 *   iterations: [in] Number of random repetitions of the split.
 *   r2: [out] Empty vector created with ``initDVector``: receives ``npc``
 *     values, element ``a-1`` being Q² with ``a`` components.
 *   s: [in] Optional stop signal (see :c:func:`PCA`); ``NULL`` if unused.
 *
 * Returns:
 *   Nothing, the values are written into ``r2``.
 */
void PCARankValidation(matrix *mx,
                       size_t npc,
                       size_t scaling,
                       size_t group,
                       size_t iterations,
                       dvector *r2,
                       scisignal *s);

/**
 * Compute the squared prediction error (SPE) of each object and its
 * per-feature contributions.
 *
 * Despite the name, the statistic is the SPE (Q residual), not Hotelling's
 * T². Each object is projected and reconstructed with ``npc`` components;
 * the residuals are computed in the preprocessed (centered and scaled)
 * space. ``contributions[i][j]`` is the squared residual of feature ``j`` for
 * object ``i`` and ``spe[i]`` is their sum. Large values identify objects
 * (and the features responsible) that the model does not describe well.
 *
 * Args:
 *   x: [in] Objects to evaluate, ``n x f``.
 *   model: [in] Model fitted with :c:func:`PCA`.
 *   npc: [in] Number of components of the model to use.
 *   spe: [out] Vector created with ``initDVector``; resized to ``n``.
 *   contributions: [out] Matrix created with ``initMatrix``; resized to
 *     ``n x f``.
 *
 * Returns:
 *   Nothing, the results are written into ``spe`` and ``contributions``.
 */
void PCATsqContributions(
  matrix *x,
  PCAMODEL *model,
  size_t npc,
  dvector *spe,
  matrix *contributions);

/**
 * Compute the residual matrix of the training data after ``pc`` components.
 *
 * Returns ``E = X_preprocessed - T P'`` using the first ``pc`` components.
 * It uses the training scores stored in the model, so ``mx`` must be the
 * same matrix used to fit the model; for new objects use
 * :c:func:`PCAScorePredictor` and :c:func:`PCAIndVarPredictor`.
 *
 * Args:
 *   mx: [in] The training matrix, ``n x f``.
 *   model: [in] Model fitted on ``mx`` with :c:func:`PCA`.
 *   pc: [in] Number of components to remove; capped to ``f``.
 *   rmx: [out] Matrix created with ``initMatrix``; resized and filled with
 *     the residuals in the preprocessed (centered and scaled) space, ``n x f``.
 *
 * Returns:
 *   Nothing, the residuals are written into ``rmx``.
 */
void GetResidualMatrix(matrix *mx,
                       PCAMODEL *model,
                       size_t pc,
                       matrix *rmx);

/**
 * Print the explained variance, scores, loadings and DModX of a model to
 * standard output.
 *
 * Args:
 *   m: [in] The model to print.
 */
void PrintPCA(PCAMODEL *m);

#endif

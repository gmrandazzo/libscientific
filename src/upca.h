/* Implements UPCA algorithms.
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

#ifndef UPCA_H
#define UPCA_H
#include "list.h"
#include "tensor.h"
#include "matrix.h"
#include "vector.h"
#include "preprocessing.h"
#include "scientificinfo.h"

#define UPCACONVERGENCE 1e-10

/**
 * Multi-way (unfolded) Principal Component Analysis model.
 *
 * The input is a tensor of ``K`` layers (``X->m[k]``), all ``n x f``: ``n``
 * objects, ``f`` features measured in ``K`` conditions. Create the model with
 * :c:func:`NewUPCAModel`, fill it with :c:func:`UPCA` and free it with
 * :c:func:`DelUPCAModel`.
 *
 * Members:
 *   scores: Scores, ``n x npc``.
 *   loadings: One ``f x K`` loading matrix per component
 *     (``loadings->m[pc]``), constrained to rank 1 (outer product of a
 *     feature profile and a layer profile).
 *   varexp: Percentage of variance explained by each component.
 *   colaverage: Column means of each layer (one vector of size ``f`` per layer).
 *   colscaling: Column scalings of each layer.
 */
typedef struct{
  matrix *scores;
  tensor *loadings;
  dvector *varexp;
  dvectorlist *colaverage;
  dvectorlist *colscaling;
} UPCAMODEL;

/**
 * Allocate an empty UPCA model.
 *
 * Args:
 *   m: [out] Address of the pointer that receives the new model.
 *
 * Note:
 *   Release the model with :c:func:`DelUPCAModel`.
 */
void NewUPCAModel(UPCAMODEL **m);

/**
 * Free a UPCA model and all its members, and set the pointer to ``NULL``.
 *
 * Args:
 *   m: [in,out] Address of the model pointer created with :c:func:`NewUPCAModel`.
 */
void DelUPCAModel(UPCAMODEL **m);

/**
 * Check that all the layers of a tensor have the same shape.
 *
 * Args:
 *   X: [in] The tensor to check.
 *
 * Returns:
 *   0 if every layer has the same number of rows and columns, -1 otherwise.
 */
int CheckTensor(tensor *X);

/**
 * Fit a multi-way PCA model (Wold, Geladi and Esbensen, J. Chemometrics 1,
 * 41-56, 1987).
 *
 * Every column of every layer is preprocessed, then the components are
 * extracted one at a time with NIPALS (up to 1000 iterations each). The
 * loading matrix of each component is constrained to rank 1, so the
 * components are not exactly orthogonal: with several components the
 * explained variances can add up to slightly more than 100%. Fitting an
 * already fitted model replaces its results.
 *
 * Args:
 *   X: [in] Layers, all ``n x f`` (see :c:func:`CheckTensor`).
 *   npc: [in] Number of components; capped to ``f``.
 *   autoscaling: [in] Column preprocessing: 0 = mean centering only,
 *     1 = standard deviation, 2 = root mean square, 3 = Pareto, 4 = range,
 *     5 = level.
 *   m: [out] Model created with :c:func:`NewUPCAModel`; filled with the
 *     results (see :c:type:`UPCAMODEL`).
 *   s: [in] Optional stop signal; ``NULL`` if unused.
 *
 * Returns:
 *   Nothing, the results are written into ``m``. Nothing is computed if the
 *   layers have different shapes.
 */
void UPCA(tensor *X,
          size_t npc,
          size_t autoscaling,
          UPCAMODEL *m,
          scisignal *s);

/**
 * Project new objects into a UPCA model.
 *
 * Projecting the training tensor reproduces ``model->scores``.
 *
 * Args:
 *   X: [in] Layers of the new objects: same number of layers and features as
 *     the training data, ``m`` rows each.
 *   model: [in] Model fitted with :c:func:`UPCA`.
 *   npc: [in] Number of components to use; capped to the model components.
 *   pscores: [out] Matrix created with ``initMatrix``; resized to ``m x npc``.
 *
 * Returns:
 *   Nothing, the scores are written into ``pscores``.
 */
void UPCAScorePredictor(tensor *X,
                        UPCAMODEL *model,
                        size_t npc,
                        matrix *pscores);

/**
 * Rebuild layers in the original units from scores and loadings.
 *
 * Computes ``X[k][i][j] = sum_pc T[i][pc] P[pc][j][k]`` with the first
 * ``npc`` components, then reverts the scaling and centering of every layer.
 *
 * Args:
 *   T: [in] Scores, ``n x a`` (``model->scores`` or the output of
 *     :c:func:`UPCAScorePredictor`).
 *   P: [in] Loadings (``model->loadings``).
 *   colaverage: [in] Column means of each layer (``model->colaverage``).
 *   colscaling: [in] Column scalings of each layer (``model->colscaling``).
 *   npc: [in] Number of components to use; capped to the available ones.
 *   X: [out] Tensor created with ``initTensor``; receives one ``n x f``
 *     matrix per layer.
 *
 * Returns:
 *   Nothing, the result is written into ``X``.
 */
void UPCAIndVarPredictor(matrix *T,
                         tensor *P,
                         dvectorlist *colaverage,
                         dvectorlist *colscaling,
                         size_t npc,
                         tensor *X);

/**
 * Print the content of a UPCA model to standard output.
 *
 * Args:
 *   m: [in] The model to print.
 */
void PrintUPCAModel(UPCAMODEL *m);

#endif

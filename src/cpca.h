/* Implements Consensus Principal Component Analysis (CPCA).
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

#ifndef CPCA_H
#define CPCA_H
#include <stdio.h>
#include "tensor.h"
#include "list.h"
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define CPCACONVERGENCE 1e-10

/**
 * Consensus Principal Component Analysis model (multi-block PCA).
 *
 * Create it with :c:func:`NewCPCAModel`, fill it with :c:func:`CPCA` and free
 * it with :c:func:`DelCPCAModel`. All the members are owned by the model.
 * With ``n`` objects, ``B`` blocks (block ``b`` has ``f_b`` columns) and
 * ``npc`` components:
 *
 * Members:
 *   block_scores: One ``n x B`` matrix per component (``block_scores->m[k]``):
 *     column ``b`` holds the scores of block ``b`` on component ``k``.
 *   block_loadings: One ``f_b x npc`` matrix per block
 *     (``block_loadings->m[b]``).
 *   super_scores: Consensus scores, ``n x npc``. They equal the PCA scores
 *     of the concatenated blocks, each divided by ``sqrt(f_b)``.
 *   super_weights: Contribution of each block to each component, ``B x npc``
 *     (unit-length columns).
 *   scaling_factor: Block scaling ``sqrt(f_b)`` of each block (size ``B``).
 *   total_expvar: Percentage of the total variance explained by each
 *     component (size ``npc``, not cumulative).
 *   block_expvar: One vector of size ``B`` per component: the cumulative
 *     percentage of the variance of each block explained so far.
 *   colaverage: Column means of each block (one vector per block).
 *   colscaling: Column scalings of each block (one vector per block).
 */
typedef struct {
  tensor *block_scores;
  tensor *block_loadings;
  matrix *super_scores;
  matrix *super_weights;
  dvector *scaling_factor;
  dvector *total_expvar;
  dvectorlist *block_expvar;
  dvectorlist *colaverage;
  dvectorlist *colscaling;
} CPCAMODEL;

/**
 * Allocate an empty CPCA model.
 *
 * Args:
 *   m: [out] Address of the pointer that receives the new model.
 *
 * Note:
 *   Release the model with :c:func:`DelCPCAModel`.
 */
void NewCPCAModel(CPCAMODEL **m);

/**
 * Free a CPCA model and all its members, and set the pointer to ``NULL``.
 *
 * Args:
 *   m: [in,out] Address of the model pointer created with :c:func:`NewCPCAModel`.
 */
void DelCPCAModel(CPCAMODEL **m);

/**
 * Fit a Consensus PCA model on a set of blocks.
 *
 * Every block is preprocessed column by column according to ``scaling`` and
 * divided by ``sqrt(f_b)``; the consensus components are then extracted one
 * at a time (Westerhuis, Kourti and MacGregor, J. Chemometrics 12, 301-321,
 * 1998), with up to 1000 iterations per component. Fitting an already
 * fitted model replaces its results.
 *
 * Args:
 *   x: [in] The blocks: one tensor matrix per block (``x->m[b]``, ``n x f_b``),
 *     all with the same rows in the same object order.
 *   scaling: [in] Column preprocessing: 0 = mean centering only, 1 = standard
 *     deviation, 2 = root mean square, 3 = Pareto, 4 = range, 5 = level.
 *   npc: [in] Number of components; capped to the columns of the smallest
 *     block.
 *   model: [out] Model created with :c:func:`NewCPCAModel`; filled with the
 *     results (see :c:type:`CPCAMODEL`).
 *
 * Returns:
 *   Nothing, the results are written into ``model``.
 *
 * Example:
 *   .. code-block:: c
 *
 *      tensor *blocks;
 *      CPCAMODEL *model;
 *      NewTensor(&blocks, 2);
 *      NewTensorMatrix(blocks, 0, 20, 10);   // block 0: 20 objects x 10 columns
 *      NewTensorMatrix(blocks, 1, 20, 7);    // block 1: 20 objects x 7 columns
 *      // ... fill blocks->m[b]->data[i][j] ...
 *      NewCPCAModel(&model);
 *      CPCA(blocks, 1, 3, model);
 *      // model->super_scores: 20 x 3, model->super_weights: 2 x 3
 *      DelCPCAModel(&model);
 *      DelTensor(&blocks);
 */
void CPCA(tensor *x, int scaling, size_t npc, CPCAMODEL *model);


/**
 * Project new objects into a CPCA model.
 *
 * Projecting the training blocks reproduces ``model->super_scores`` within
 * the convergence tolerance of the algorithm.
 *
 * Args:
 *   x: [in] The blocks of the new objects: same number of blocks and
 *     columns, in the same order, as the training data (``m`` rows each).
 *   model: [in] Model fitted with :c:func:`CPCA`.
 *   npc: [in] Number of components to use; capped to the model components.
 *   p_super_scores: [out] Matrix created with ``initMatrix``; resized to
 *     ``m x npc``.
 *   p_block_scores: [out] Tensor created with ``initTensor``; receives one
 *     ``m x B`` matrix per component (same layout as
 *     ``model->block_scores``).
 *
 * Returns:
 *   Nothing, the scores are written into the outputs.
 */
void CPCAScorePredictor(tensor *x,
                        CPCAMODEL *model,
                        size_t npc,
                        matrix *p_super_scores,
                        tensor *p_block_scores);

/**
 * Print the content of a CPCA model to standard output.
 *
 * Args:
 *   m: [in] The model to print.
 */
void PrintCPCA(CPCAMODEL *m);

#endif

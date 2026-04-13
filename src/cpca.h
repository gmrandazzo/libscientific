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

#define CPCACONVERGENCE 1e-18

/**
 * CPCA model data structure.
 * 
 * - **block_scores** matrix of scores
 * - **block_loadings** matrix of loadings
 * - **super_scores** matrix of super scores
 * - **super_weights** matrix of super weigths
 * - **scaling_factor** dvector of scaling factors
 * - **total_expvar** dvector of total explained variance
 * - **block_expvar** dvector list of block explained variance
 * - **colaverage** dvector list of column average
 * - **colscaling** dvector list of column scaling
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
 * Inizialize a new CPCA Model
 */
void NewCPCAModel(CPCAMODEL **m);

/**
 * Delete a new CPCA Model
 */
void DelCPCAModel(CPCAMODEL **m);

/**
 * Consensus Principal Component Analysis
 * 
 * @param [in] x libscientific tensor data input
 * @param [in] scaling scaling type expressed as unsigned int type
 * @param [in] npc number of desired principal components
 * @param [out] model initialized model using NewCPCAModel(...). The datastructure will be populated with results
 * 
 * 
 * Available scalings:
 *
 * - 0: No scaling. Only mean centering
 *
 * - 1: Mean centering and STDEV scaling
 *
 * - 2: Mean centering and Root-Mean-Square column scaling
 *
 * - 3: Mean centering and Pareto scaling
 *
 * - 4: Mean centering and min-max range scaling
 *
 * - 5: Mean centering and level scaling
 *
 */
void CPCA(tensor *x, int scaling, size_t npc, CPCAMODEL *model);


/**
 * Project objects in a CPCA model.
 * 
 * @param [in] x libscientific tensor data input
 * @param [in] model CPCA model
 * @param [in] npc number of desired principal components
 * @param [out] p_super_scores predicted super scores
 * @param [out] p_block_scores predicted block of scores
 */
void CPCAScorePredictor(tensor *x,
                        CPCAMODEL *model,
                        size_t npc,
                        matrix *p_super_scores,
                        tensor *p_block_scores);

/**
 * @brief Print CPCAMODEL to video.
 *
 * @param [in] m computed cpca model
 *
 * @par Returns
 *    Nothing.
 */
void PrintCPCA(CPCAMODEL *m);

#endif

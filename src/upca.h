/* upca.h
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
 * UPCA model data structure.
 * 
 * - **scores** matrix of scores
 * - **loadings** tesnor of loadings
 * - **varexp** vector of explained variance by every component 
 * - **colaverage** vector list of column average
 * - **colaverage** vector list of column scaling
 */
typedef struct{
  matrix *scores;
  tensor *loadings;
  dvector *varexp;
  dvectorlist *colaverage;
  dvectorlist *colscaling;
} UPCAMODEL;

/**
 *  Initialize an empty UPCAMODEL 
 */
void NewUPCAModel(UPCAMODEL **m);

/**
 * Delete an UPCAMODEL 
 */
void DelUPCAModel(UPCAMODEL **m);

int CheckTensor(tensor *X_);

/**
 * Unfolded Principal Component Analysis
 * 
 * @param [in] X_ input tensor
 * @param [in] npc number of desired principal components
 * @param [in] autoscaling
 * @param [out] m initialized model using NewUPCAModel(...). The datastructure will be populated with results
 * @param [in] ssignal libscientific signal. Default value is NULL
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
void UPCA(tensor *X_,
          size_t npc,
          size_t autoscaling,
          UPCAMODEL *m,
          ssignal *s);

/**
 *  Predict scores given an unfolded principal component analysis and a tensor as input. 
 * 
 * @param [in] X_ input tensor
 * @param [in] model computed UPCAMODEL
 * @param [in] npc number of desired principal components
 * @param [out] pscores predicted scores
 * 
 */
void UPCAScorePredictor(tensor *X_,
                        UPCAMODEL *model,
                        size_t npc,
                        matrix *pscores);

/**
 * @brief Reconstruct the original input tensor from UPCA model using scores and loadings.
 *
 * @param [in] t upca scores with size #objects x npc
 * @param [in] p upca loadings with size order x npc x #features
 * @param [in] colaverage input list of column average with size #features
 * @param [in] colscaling input list of column scaling with size #features
 * @param [in] npc desired principal components to use for the tensor reconstruction
 * @param [in] X ouptut reconstructed tensor
 *
 * @par Returns
 *    Nothing.
 */
void UPCAIndVarPredictor(matrix *T,
                         tensor *P,
                         dvectorlist *colaverage,
                         dvectorlist *colscaling,
                         size_t npc,
                         tensor *X);

/**
 * @brief Print UPCAMODEL to video.
 *
 * @param [in] m computed upca model
 *
 * @par Returns
 *    Nothing.
 */
void PrintUPCAModel(UPCAMODEL *m);

#endif

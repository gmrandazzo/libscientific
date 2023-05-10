/* pca.h
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

#ifndef PCA_H
#define PCA_H
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define PCACONVERGENCE 1e-10

/**
 * PCA model data structure.
 * 
 * - **scores** matrix of scores
 * - **loadings** matrix of loadings
 * - **varexp** vector of explained variance by every component 
 * - **colaverage** input matrix column average
 * - **colaverage** input matrix column scaling
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
 * Initialize an empty PCAMODEL 
 */
void NewPCAModel(PCAMODEL **m);

/**
 * Delete a PCAMODEL 
 */
void DelPCAModel(PCAMODEL **m);

void calcVarExpressed(double ss, dvector *eval, dvector *varexp);

/*
 * Calculate the object distance.
 * THIS METHOD IS USED INSIDE THE PCA ALGORITHM
 */
double calcObjectDistance(matrix *m);

/*
 * Calculate convergence criteria as described in 
 * DOI: 10.1002/cem.1180010107 page 51 d = ...
 * THIS METHOD IS USED INSIDE THE PCA/CPCA/PLS ALGORITHM
 */
double calcConvergence(dvector *t_new, dvector *t_old);


/**
 * @brief Calculate a principal component analysis using the NIPALS algorithm.
 *
 * @param [in] mx libscientific matrix data input 
 * @param [in] scaling scaling type expressed as unsigned int type
 * @param [in] npc number of desired principal components
 * @param [out] PCAMODEL initialized model using NewPCAModel(...). The datastructure will be populated with results
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
 * @par Returns
 *    Nothing.
 */
void PCA(matrix *mx,
         int scaling,
         size_t npc,
         PCAMODEL *model,
         ssignal *s);


/**
 * @brief Predict scores given a principal component analysis and a matrix as input.
 *
 * @param [in] mx libscientific matrix data input 
 * @param [in] PCAMODEL pca model input
 * @param [in] npc number of desired principal components
 * @param [in] pscores predicted scores 
 *
 * @par Returns
 *    Nothing.
 */
void PCAScorePredictor(matrix *mx,
                       PCAMODEL *model,
                       size_t npc,
                       matrix *pscores);

/**
 * @brief Reconstruct the original input matrix from PCA model using scores and loadings.
 *
 * @param [in] t pca scores with size #objects x npc
 * @param [in] p pca loadings with size npc x #features
 * @param [in] colaverage input column average with size #features
 * @param [in] colscaling input column scaling with size #features
 * @param [in] npc desired principal components to use for the matrix reconstruction
 * @param [in] mx ouptut reconstructed matrix
 *
 * @par Returns
 *    Nothing.
 */
void PCAIndVarPredictor(matrix *t,
                        matrix *p,
                        dvector *colaverage,
                        dvector *colscaling,
                        size_t npc,
                        matrix *mx);

void PCARSquared(matrix *mx,
                 PCAMODEL *model,
                 size_t npc,
                 dvector *r2);

void PCARankValidation(matrix *mx,
                       size_t npc,
                       size_t scaling,
                       size_t group,
                       size_t iterations,
                       dvector *r2,
                       ssignal *s);

/**
 * @brief Compute the residual matrix for a specific number of component.
 *
 * @param [in] mx original input matrix
 * @param [in] model computed pca model
 * @param [in] pc max component to extract the residual matrix
 * @param [in] rmx residual matrix of output
 *
 * @par Returns
 *    Nothing.
 */
void GetResidualMatrix(matrix *mx,
                       PCAMODEL *model,
                       size_t pc,
                       matrix *rmx);

/**
 * @brief Print PCAMODEL to video.
 *
 * @param [in] m computed pca model
 *
 * @par Returns
 *    Nothing.
 */
void PrintPCA(PCAMODEL *m);

#endif

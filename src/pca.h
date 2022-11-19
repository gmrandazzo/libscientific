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

#define PCACONVERGENCE 1e-8

typedef struct{
  matrix *scores;
  matrix *loadings;
  matrix *dmodx;
  dvector *varexp;
  dvector *colaverage;
  dvector *colscaling;
} PCAMODEL;

void NewPCAModel(PCAMODEL **m);
void DelPCAModel(PCAMODEL **m);

void calcVarExpressed(double ss, dvector *eval, dvector *varexp);
double calcObjectDistance(matrix *m);

/*
 * Calculate a principal component analysis
 */
void PCA(matrix *mx,
         int scaling,
         size_t npc,
         PCAMODEL *model,
         ssignal *s);

/*
 * Predict scores given a principal component analysis and a matrix as input
 */
void PCAScorePredictor(matrix *mx,
                       PCAMODEL *model,
                       size_t npc,
                       matrix *pscores);

/*
 * Reconstruct the original matrix from PCA model, scores and loadings
 */
void PCAIndVarPredictor(matrix *t,
                        matrix *p,
                        dvector *colaverage,
                        dvector *colsdev,
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
                       dvector *r,
                       ssignal *s);

/*Compute the residual matrix for a specific number of component.
 * mx = matrix of origin model
 * model = PCA model where are stored scores and loadings
 * pc = max component to extract the residual matrix
 * rmx = residual matrix of output. must be initialized with initMatrix(&rmx)
 */
void GetResidualMatrix(matrix *mx,
                       PCAMODEL *model,
                       size_t pc,
                       matrix *rmx);

void PrintPCA(PCAMODEL *m);

#endif

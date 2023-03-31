/* cpca.h
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

#ifndef CPCA_H
#define CPCA_H
#include <stdio.h>
#include "tensor.h"
#include "list.h"
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define CPCACONVERGENCE 1e-18

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

void NewCPCAModel(CPCAMODEL **m);
void DelCPCAModel(CPCAMODEL **m);

/*
 * Consensus Principal Component Analysis
 *
 * ANALYSIS OF MULTIBLOCK AND HIERARCHICAL PCA AND PLS MODELS
 * JOHAN A. WESTERHUIS, THEODORA KOURTI* AND JOHN F. MACGREGOR
 * J. Chemometrics 12, 301–321 (1998)
 *
 * N.B.: The superscores of CPCA are the scores of a PCA!!
 */
void CPCA(tensor *x, int scaling, size_t npc, CPCAMODEL *model);


/*
 * Project objects in a CPCA model.
 */
void CPCAScorePredictor(tensor *x,
                        CPCAMODEL *model,
                        size_t npc,
                        matrix *p_super_scores,
                        tensor *p_block_scores);

void PrintCPCA(CPCAMODEL *m);

#endif

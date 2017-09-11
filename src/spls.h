/* spls.c
* Serial PLS
* Copyright (C) <2017>  Giuseppe Marco Randazzo
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

#ifndef SPLS_H
#define SPLS_H
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

typedef struct{
  matrix *xscores;
  matrix *xloadings;
  matrix *xweights;
  matrix *yscores;
  matrix *yloadings;
  matrix *cweights;
  dvector *b;
  dvector *xvarexp;
  dvector *xcolaverage;
  dvector *xcolscaling;
  dvector *ycolaverage;
  dvector *ycolscaling;
  matrix *r2y_model; /* each column correspond to an y dependent variable and each row correspond to a principal component*/
  matrix *r2y_validation;
  matrix *recalc_residuals;
  matrix *q2y;
  matrix *sdep; /* Standard Deviation over Prediction */
  matrix *sdec; /* Standard Deviation over Recalculating */
  matrix *bias;
  matrix *recalculated_y;
  matrix *predicted_y;
  matrix *pred_residuals;
  matrix *r2q2scrambling;
  matrix *q2_sample_validation;
  matrix *sdep_sample_validation;
  matrix *q2_sample_validation_surface;
  matrix *sdep_sample_validation_surface;
} SPLSMODEL;

/*
 * Description: Create a new data structure
 * to define the number of blocks of the x matrix
 */
typedef struct{
  size_t from;
  size_t to;
} BLOCKS;

/*
 * Description: Create n BLOCKS
 */
void SetNBlocks(BLOCKS **b, size_t n);

/*
 * Description: Destroy n BLOCKS
 */
void DelBlocks( BLOCKS **b);

/*
 * Description: Create an SPLSMODEL
 */
void NewSPLSModel(SPLSMODEL **m);

/*
 * Description: Destroy an SPLSMODEL
 */
void DelSPLSModel(SPLSMODEL **m);

/*
 * Description: Compute the SPLS
 *
 * Algorithm for n blocks X1, X2, ..., Xn
 * X1 = T1 x P1' + E1
 * X2 = T2 x P2' + E2
 * Xn = Tn x Pn' + En
 * Y = T1 x C1' + T2 x C2' + ... + Tn x Cn' + F
 * The regression model can be divided into n parts...
 *
 * Y = X1 x B1 + X2 x B2 + ... + Xn x Bn + F
 * N.B.: B1 and B2 are regression matrices...
 */
void SPLS(matrix *mx, matrix *my, BLOCKS *b, size_t nlv, size_t xautoscaling, size_t yautoscaling, SPLSMODEL *m, ssignal *s);

#endif

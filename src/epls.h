/* epls.h
*
* Copyright (C) <2018>  Giuseppe Marco Randazzo
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

#ifndef EPLS_H
#define EPLS_H
#include "matrix.h"
#include "pls.h"
#include "vector.h"
#include "scientificinfo.h"

typedef struct{
  PLSMODEL **models;
  matrix **q2;
  matrix **sdep;
  size_t n_models;
  size_t nlv;
  size_t ny;
} EPLSMODEL;

/*
 * Description: Create a new PLSMODEL
 */
void NewEPLSModel(EPLSMODEL **m);

/*
 * Description: Delete a PLSMODEL
 */
void DelEPLSModel(EPLSMODEL **m);

typedef enum{
    Bagging = 0,
    RandomSubspaceMethod = 1
} ELearningMethod;

void EPLS(matrix *mx, matrix *my, size_t nlv, size_t xautoscaling, size_t yautoscaling, size_t n_models, double testsize, EPLSMODEL *model, ELearningMethod agtype, ssignal *s);

typedef enum{
  Mean = 0,
  Median = 1,
  WeightedMedian = 2 // weight according the error SDEP of the model
} CombinationRule;

void EPLSGetSXScore(EPLSMODEL *m, CombinationRule crule, matrix *sxscores);
void EPLSGetSXLoadings(EPLSMODEL *m, CombinationRule crule, matrix *sxloadings);
void EPLSGetSYScore(EPLSMODEL *m, CombinationRule crule, matrix *syscores);
void EPLSGetSYLoadings(EPLSMODEL *m, CombinationRule crule, matrix *syloadings);
void EPLSGetSWeights(EPLSMODEL *m, CombinationRule crule, matrix *sweights);
void EPLSGetSBetaCoefficients(EPLSMODEL *m, CombinationRule crule, matrix *sbetas);

void EPLSYPrediction(matrix *mx, EPLSMODEL *m, CombinationRule crule, matrix *py);

void PrintEPLSModel(EPLSMODEL *m);

#endif

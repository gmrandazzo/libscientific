/* lda.h
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

#ifndef LDA_H
#define LDA_H
#include <stdio.h>
#include "matrix.h"
#include "vector.h"
#include "tensor.h"
#include "scientificinfo.h"


typedef struct{
  matrix *inv_cov;
  tensor *features;
  tensor *mnpdf;
  matrix *evect;
  matrix *mu;
  matrix *fmean;
  matrix *fsdev;
  dvector *eval;
  dvector *pprob;
  uivector *classid;
  size_t nclass;
  size_t class_start;
  dvector *sens, *spec, *ppv, *npv, *acc;
} LDAMODEL;

void NewLDAModel(LDAMODEL **m);
void DelLDAModel(LDAMODEL **m);
void PrintLDAModel(LDAMODEL *m);

void LDA(matrix *mx, uivector *y, LDAMODEL *lda);
/*prediction
 * OUTPUT:
 *  - predicted features
 *  - probability
 *  - multivariate normal probability distribution of features
 *  - class prediction
 */
void LDAPrediction(matrix *mx, LDAMODEL *lda, matrix **pfeatures, matrix **probability, matrix **mnpdf, uivector **prediction);

void LDAStatistics(dvector *y_true, dvector *y_score, matrix **roc, double *roc_auc, matrix **precision_recal, double *pr_auc);

void LDARandomGroupsCV(matrix *mx, uivector *my, size_t group, size_t iterations, dvector **sens, dvector **spec, dvector **ppv, dvector **npv, dvector **acc, size_t nthreads, ssignal *s);

void LDALOOCV(matrix* mx, uivector* my, dvector** sens, dvector** spec, dvector** ppv, dvector** npv, dvector **acc, size_t nthreads, ssignal *s);

#endif

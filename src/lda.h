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
  matrix *yscrambling;
  matrix *recalculated_y;
  matrix *recalculated_residuals;
  matrix *predicted_y;
  matrix *predicted_residuals;
  dvector *eval;
  dvector *pprob;
  uivector *classid;
  size_t nclass;
  size_t class_start;
  tensor *roc;
  tensor *pr;
  dvector *roc_aucs;
  dvector *pr_aucs;
} LDAMODEL;

/**
 * Initialize and LDAMODEL
 */
void NewLDAModel(LDAMODEL **m);

/**
 * Delete and LDAMODEL
 */
void DelLDAModel(LDAMODEL **m);

/**
 * @brief Print LDAMODEL to video.
 *
 * @param [in] m computed lda model
 *
 * @par Returns
 *    Nothing.
 */
void PrintLDAModel(LDAMODEL *m);

/**
 * Calculate and LDA Fisher discriminant analyisis model
 * 
 * @param [in] mx input matrix of independent variables
 * @param [in] my input matrix of classes
 * @param [out] lda initialized model using NewLDAModel(...). The datastructure will be populated with results
 */
void LDA(matrix *mx, matrix *my, LDAMODEL *lda);

/** Predict classes using an LDA fisher discriminant analysis model
 * 
 * @param [in] mx input matrix of independent variables
 * @param [in] lda computed LDAMODEL
 * @param [out] pfeatures predicted features
 * @param [out] probabilities
 * @param [out] multivariate normal probability distribution of features
 * @param [out] predicted classes
 * 
 */
void LDAPrediction(matrix *mx,
                   LDAMODEL *lda,
                   matrix *pfeatures,
                   matrix *probability,
                   matrix *mnpdf,
                   matrix *prediction);

/**
 * Binary statistics
 */
void LDAStatistics(dvector *y_true,
                   dvector *y_pred,
                   matrix *roc,
                   double *roc_auc,
                   matrix *precision_recal,
                   double *pr_auc);

/**
 *  Multiclass statistics
 */
void LDAMulticlassStatistics(matrix *y_true,
                             matrix *y_pred,
                             tensor *roc,
                             dvector *roc_aucs,
                             tensor *precision_recals,
                             dvector *pr_aucs);

int getNClasses(matrix *my);

#endif

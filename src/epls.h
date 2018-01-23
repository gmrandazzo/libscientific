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
  uivector **model_feature_ids;
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

/* Description: Extract a subspace matrix giving the featureids */
void SubspaceMatrix(matrix *mx, uivector *featureids, matrix **x_subspace);

typedef enum{
    Bagging = 0,
    FixedRandomSubspaceMethod = 1,
    BaggingRandomSubspaceMethod = 2
} ELearningMethod;

typedef struct{
  ELearningMethod algorithm;
  size_t n_models;
  double trainsize;
  size_t r_fix;
} ELearningParameters;

void EPLS(matrix *mx, matrix *my, size_t nlv, size_t xautoscaling, size_t yautoscaling, EPLSMODEL *m, ELearningParameters eparm, ssignal *s);

typedef enum{
  Averaging = 0, /* Hard voting */
  Median = 1 /* Hard voting */
  /*WeightedMedian = 2  Soft voting */
} CombinationRule;

void EPLSGetSXScore(EPLSMODEL *m, CombinationRule crule, matrix *sxscores);
void EPLSGetSXLoadings(EPLSMODEL *m, CombinationRule crule, matrix *sxloadings);
void EPLSGetSYScore(EPLSMODEL *m, CombinationRule crule, matrix *syscores);
void EPLSGetSYLoadings(EPLSMODEL *m, CombinationRule crule, matrix *syloadings);
void EPLSGetSWeights(EPLSMODEL *m, CombinationRule crule, matrix *sweights);
void EPLSGetSBetaCoefficients(EPLSMODEL *m, CombinationRule crule, matrix *sbetas);

/* Description: */
void EPLSYPRedictorAllLV(matrix *mx, EPLSMODEL *m, CombinationRule crule, array **tscores, matrix **py);

/*
 * Description: Calculate the correlation coefficient (ccoeff),
 *              the standard deviation of the prediction (stdev),
 *              the bias of the prediction (bias) in a regression model.
 *              mx and my could be the training or the test datasets.
 */
void EPLSRegressionStatistics(matrix *my_true, matrix *my_pred, matrix** ccoeff, matrix **stdev, matrix **bias);

/*
 * Description: Calculate the roc curve, the auc, the precision recall curve,
 *              the precision_recall_auc of a classification model.
 *              mx and my could be the training or the test datasets.
 */
void EPLSDiscriminantAnalysisStatistics(matrix *my_true, matrix *my_score, array **roc, matrix **roc_auc, array **precision_recall, matrix **precision_recall_ap);

void PrintEPLSModel(EPLSMODEL *m);

#endif

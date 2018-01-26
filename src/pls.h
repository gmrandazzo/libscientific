/* pls.h
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

#ifndef PLS_H
#define PLS_H
#include "array.h"
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define PLSCONVERGENCE 1e-8

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
  matrix *recalculated_y;
  matrix *recalc_residuals;
  matrix *predicted_y;
  matrix *pred_residuals;
  /* Regression variables */
  matrix *r2y_model; /* each column correspond to an y dependent variable and each row correspond to a principal component*/
  matrix *r2y_validation;
  matrix *q2y;
  matrix *sdep; /* Standard Deviation over Prediction */
  matrix *sdec; /* Standard Deviation over Recalculating */
  matrix *bias;

  /* Discriminant Analyisis variables */
  array *roc_model;
  array *roc_validation;
  matrix *roc_auc_model;
  matrix *roc_auc_validation;
  array *precision_recall_model;
  array *precision_recall_validation;
  matrix *precision_recall_ap_model;
  matrix *precision_recall_ap_validation;

  matrix *yscrambling;
} PLSMODEL;

/*
 * Description: Create a new PLSMODEL
 */
void NewPLSModel(PLSMODEL **m);

/*
 * Description: Delete a PLSMODEL
 */
void DelPLSModel(PLSMODEL **m);

/*
 * Calculate the latent variables according the NIPALS algorithm
 * The function use X, Y and return
 * - deflated X
 * - deflated Y
 * - scores t, u
 * - loadings p, q
 * - weights w (or c)
 * - the beta coefficient bcoef
 *
 * N.B.: all the dvector must be initialised.
 * See Geladi ref. for details
 */
void LVCalc(matrix **X, matrix **Y, dvector **t, dvector **u, dvector **p, dvector **q, dvector **w, double *bcoef);

/*
 * Description PLS calculation from P. Geladi algorithm
 */
void PLS(matrix *mx, matrix *my, size_t nlv, size_t xautoscaling, size_t yautoscaling, PLSMODEL *model, ssignal *s);

/*
 * Description: Calculate betas coefficients from a pls model at specific nlv latent variables
 */
void PLSBetasCoeff(PLSMODEL *model, size_t nlv, dvector **betas);

/*
 * Description: Project a matrix and predict the scores into the new space.
 * This function is used before predict the Y values
 */
void PLSScorePredictor(matrix *mx, PLSMODEL *model, size_t nlv, matrix **xscores);

/*
 * Description: Calculate the Y values at a specific lv number (nlv).
 * Output: y shape (y->row, y-col)
 */
void PLSYPredictor(matrix *tscore, PLSMODEL *model, size_t nlv, matrix **y);

/*
 * Description: Calculate the Y values at all the lv.
 * Output: y shape (y->row, y->col*nlv)
 */
void PLSYPredictorAllLV(matrix *mx, PLSMODEL *model, matrix **tscores, matrix **y);


/*
 * Description: Calculate the correlation coefficient (ccoeff),
 *              the standard deviation of the prediction (stdev),
 *              the bias of the prediction (bias) in a regression model.
 *              mx and my could be the training or the test datasets.
 */
void PLSRegressionStatistics(matrix *my_true, matrix *my_pred, matrix** ccoeff, matrix **stdev, matrix **bias);

/*
 * Description: Calculate the roc curve, the auc, the precision recall curve,
 *              the precision_recall_auc of a classification model.
 *              mx and my could be the training or the test datasets.
 */
void PLSDiscriminantAnalysisStatistics(matrix *my_true, matrix *my_score, array **roc, matrix **roc_auc, array **precision_recall, matrix **precision_recall_ap);

/*
 * Description: Calculate the ROC curve with AUC and the Precision-Recall crurve with his
 * area usefull in PLS-DA case

void PLSBinaryClassificationScores(matrix *mx, matrix *my, PLSMODEL *model, size_t nlv, dvector** roc_auc, dvector **pr_ap, matrix **roc_curve, matrix **pr_curve);
*/

/*
 * Description: Calculate the PLS Very important variables
 */
void PLSVIP(PLSMODEL *model, matrix **vip);

/*
 * Description: Get the Cutoff based on the grow of r2, q2 in case of regression
 * or auc in case of classification.
 */
int GetLVCCutoff(matrix* coeff);

void PrintPLSModel(PLSMODEL *model);

#endif

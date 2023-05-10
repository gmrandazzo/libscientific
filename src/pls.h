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
#include "tensor.h"
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define PLSCONVERGENCE 1e-8

/**
 * PLS model data structure
 * 
 * - **xscores** x space scores
 * - **xloadings** x space loadings
 * - **xweights** x space weights
 * - **yscores** y space scores
 * - **yloadings** y space loadings
 * - **b** pls regression coefficients
 * - **xvarexp** variance explained in the x space
 * - **xcolaverage** x independent variable column average
 * - **xcolscaling** x independent variable column scaling
 * - **ycolaverage** y independent variable column average
 * - **ycolscaling** y independent variable column scaling
 * - **recalculated_y** y recalculated
 * - **recalc_residuals** y recalculated residuals
 * - **predicted_y** y predicted 
 * - **pred_residuals** y predicted residuals 
 * - **r2y_recalculated** r squared using y recalculated values
 * - **r2y_validation**
 * - **q2y** q squared using y predicted values
 * - **sdep** standard deviation over prediction using y predictions
 * - **sdec** standard deviation over recalculation using y recalculated
 * - **bias** bias
 * - **roc_recalculated** receiver operating characteristic using y recalculated 
 * - **roc_validation** receiver operating characteristic using y predicted 
 * - **roc_auc_recalculated** receiver operating characteristic area under the curve using y recalculated
 * - **roc_auc_validation** eceiver operating characteristic area under the curve using y predicted
 * - **precision_recall_recalculated** precision-recall curve using y recalculated
 * - **precision_recall_validation** precision-recall curve using y predicted
 * - **precision_recall_ap_recalculated** precision-recall aread under the curve using y recalculated
 * - **precision_recall_ap_validation** precision-recall aread under the curve using y predicted
 * - **yscrambling** y-scrambling r-squared and q-squared
 */
typedef struct{
  matrix *xscores;
  matrix *xloadings;
  matrix *xweights;
  matrix *yscores;
  matrix *yloadings;
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
  matrix *r2y_recalculated; /* each column correspond to an y dependent variable and each row correspond to a principal component*/
  matrix *r2y_validation;
  matrix *q2y;
  matrix *sdep; /* Standard Deviation over Prediction */
  matrix *sdec; /* Standard Deviation over Recalculating */
  matrix *bias;
  /* Discriminant Analyisis variables */
  tensor *roc_recalculated;
  tensor *roc_validation;
  matrix *roc_auc_recalculated;
  matrix *roc_auc_validation;
  tensor *precision_recall_recalculated;
  tensor *precision_recall_validation;
  matrix *precision_recall_ap_recalculated;
  matrix *precision_recall_ap_validation;
  matrix *yscrambling;
} PLSMODEL;

/**
 * Create a new PLSMODEL
 */
void NewPLSModel(PLSMODEL **m);

/**
 * Delete a PLSMODEL
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
void LVCalc(matrix *X, matrix *Y, dvector *t, dvector *u, dvector *p, dvector *q, dvector *w, double *bcoef);


/**
 * @brief Calculate a partial least squares model using the NIPALS algorithm.
 *
 * @param [in] mx libscientific matrix data input: x independent variables
 * @param [in] my libscientific matrix data input: y dependent variables 
 * @param [in] nlv number of desired latent variables
 * @param [in] xautoscaling scaling typeon the x independent variables expressed as unsigned int type
 * @param [in] yautoscaling scaling typeon the y dependent variables expressed as unsigned int type
 * @param [out] PLSMODEL output -initialized model using NewPLSAModel(...). The datastructure will be populated with results
 * @param [in] ssignal libscientific signal. Default value is NULL.
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
void PLS(matrix *mx,
         matrix *my,
         size_t nlv,
         int xautoscaling,
         int yautoscaling,
         PLSMODEL *model,
         ssignal *s);

/**
 * Calculate betas coefficients from a pls model at specific nlv latent variables
 * @param [in] model libscientific matrix data input: x independent variables
 * @param [in] nlv number of desired latent variables
 * @param [out] betas output - initialized libscientific dvector
 */
void PLSBetasCoeff(PLSMODEL *model,
                   size_t nlv,
                   dvector *betas);

/**
 * Project a matrix and predict the scores into the new space.
 * This function is used before predict the Y values
 * @param [in] mx libscientific matrix data input: x independent variables
 * @param [in] model PLSMODEL
 * @param [in] nlv number of desired latent variables
 * @param [out] xscores output - initialized libscientific matrix
 */
void PLSScorePredictor(matrix *mx,
                       PLSMODEL *model,
                       size_t nlv,
                       matrix *xscores);

/**
 * Calculate the Y values at a specific lv number (nlv).
 * 
 * N.B.: The output of y will be: [y->row][y-col]
 * @param [in] tscore input matrix of scores calculated using PLSScorePredictor
 * @param [in] model PLSMODEL
 * @param [in] nlv number of desired latent variables
 * @param [out] y output - initialized libscientific matrix
 */
void PLSYPredictor(matrix *tscore,
                   PLSMODEL *model,
                   size_t nlv,
                   matrix *y);

/**
 * Calculate the Y values at all the lv.
 * 
 * N.B.: The output of y will be: [y->row][y->col*nlv]
 * @param [in] mx input matrix of scores calculated using PLSScorePredictor
 * @param [in] model PLSMODEL
 * @param [out] tscores output - predicted scores
 * @param [out] y output - y predicted
 */
void PLSYPredictorAllLV(matrix *mx,
                        PLSMODEL *model,
                        matrix *tscores,
                        matrix *y);


/**
 * Calculate the correlation coefficient (ccoeff),
 * the root mean square error of the prediction (rmse),
 * the bias of the prediction (bias) in a regression model.
 * mx and my could be the training or the test datasets.
 *
 * @param [in] my_true input matrix of y true values
 * @param [in] my_pred input matrix of y predicted or recalculated values
 * @param [out] ccoeff output - correlation coefficients q-squared or r-squared if my_pred is respectivelly predicted or recalculated 
 * @param [out] rmse output - root mean square error
 * @param [out] bias output - bias defined as how distant are we from the diagonal
 */
void PLSRegressionStatistics(matrix *my_true,
                             matrix *my_pred,
                             matrix *ccoeff,
                             matrix *rmse,
                             matrix *bias);

/**
 * Calculate the roc curve, the auc, the precision recall curve,
 * the precision_recall_auc of a classification model.
 * mx and my could be the training or the test datasets.
 * 
 * @param [in] my_true input matrix of y true values
 * @param [in] my_scores input matrix of y predicted or recalculated class scores/probabilities
 * @param [out] roc output - roc curves
 * @param [out] roc_auc output - roc AUCs
 * @param [out] precision_recall output - precision-recall curves
 * @param [out] precision_recall_ap output - precision-recall AUCs
 */
void PLSDiscriminantAnalysisStatistics(matrix *my_true,
                                       matrix *my_score,
                                       tensor *roc,
                                       matrix *roc_auc,
                                       tensor *precision_recall,
                                       matrix *precision_recall_ap);

/*
 * Description: Calculate the ROC curve with AUC and the Precision-Recall crurve with his
 * area usefull in PLS-DA case

void PLSBinaryClassificationScores(matrix *mx, matrix *my, PLSMODEL *model, size_t nlv, dvector** roc_auc, dvector **pr_ap, matrix **roc_curve, matrix **pr_curve);
*/

/**
 * Calculate the PLS Very important variables
 */
void PLSVIP(PLSMODEL *model, matrix *vip);

/**
 * Get the Cutoff based on the grow of r2, q2 in case of regression
 * or auc in case of classification.
 */
int GetLVCCutoff(matrix* coeff);

/**
 * Print to video a PLSMODEL
 */
void PrintPLSModel(PLSMODEL *model);

#endif

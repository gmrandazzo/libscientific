/* Provides tools for model validation and error estimation.
 * Copyright (C) 2016-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef MODELVALIDATION_H
#define MODELVALIDATION_H
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

typedef enum{
    _PLS_ = 0,
    _PLS_DA_ = 1,
    _EPLS_ = 2,
    _EPLS_DA_ = 3,
    _MLR_ = 4,
    _LDA_ = 5
} AlgorithmType;

typedef enum{
  LOO = 0,
  BootstrapRGCV = 1
} ValidationType;

typedef struct{
  ValidationType vtype;
  size_t rgcv_group;
  size_t rgcv_iterations;
} ValidationArg;

ValidationArg initValidationArg();

typedef struct{
  matrix *mx;
  matrix *my;
  size_t nlv;
  size_t xautoscaling;
  size_t yautoscaling;
} MODELINPUT;

MODELINPUT initModelInput();

/* Description: Generate random groups of train and test set.
 *
 * Usage:
 * matrix *gid, *x_train, *y_train, *x_test, *y_test;
 * initMatrix(&gid);
 * random_kfold_group_generator(&gid, 5 (groups), 10 (objects), 1234567890 (random number));
 *
 * for(g = 0; g < gid->row; g++){ // for each group
 *  initMatrix(&x_train);
 *  initMatrix(&y_train);
 *  initMatrix(&x_test);
 *  initMatrix(&y_test);
 *  kfold_group_train_test_split(arg->mx, arg->my, gid, g, &x_train,&y_train,&x_test, &y_test); // fill in train of groups != g and fill as test the group "g"
 *
 *  do your calculations....
 *
 *  DelMatrix(&x_train);
 *  DelMatrix(&y_train);
 *  DelMatrix(&x_test);
 *  DelMatrix(&y_test);
 * }
 * DelMatrix(&gid);
 */
void random_kfold_group_generator(matrix *gid,
                                  size_t ngroups,
                                  size_t nobj,
                                  unsigned int *srand_init);

void kfold_group_train_test_split(matrix *x,
                                  matrix *y,
                                  matrix *gid,
                                  size_t group_id,
                                  matrix *x_train,
                                  matrix *y_train,
                                  matrix *x_test,
                                  matrix *y_test);

/*
 * Description: Generate random train/test splits
 * 
 * Usage:
 * matrix *x_train, *y_train, *x_test, *y_test;
 * initMatrix(&x_train);
 * initMatrix(&y_train);
 * initMatrix(&x_test);
 * initMatrix(&y_test);
 * train_test_split(mx,
 *                  my,
 *                  0.2,
 *                  &x_train,
 *                  &y_train,
 *                  &x_test,
 *                  &y_test);
 * do your calculations....
 *
 * DelMatrix(&x_train);
 * DelMatrix(&y_train);
 * DelMatrix(&x_test);
 * DelMatrix(&y_test);
 */
void train_test_split(matrix *x,
                      matrix *y,
                      double testsize,
                      matrix *x_train,
                      matrix *y_train,
                      matrix *x_test,
                      matrix *y_test,
                      uivector *testids,
                      unsigned int *srand_init);

/*
 * Description: Calculate the Bootstrap Random group cross validation.
 *              Every method has its model function. 
 * 
 *              Available methods:
 *              - PLS: partial least squares
 *              - PLS_DA: partial least squares discriminant analysis
 *              - EPLS: ensemble partial least squares
 *              - EPLS_DA: ensemble partial least squares discriminant analysis
 *              - MLR: multiple linear regression
 *              - LDA: linear discriminant analysis
 */
void BootstrapRandomGroupsCV(MODELINPUT *input,
                             size_t group,
                             size_t iterations,
                             AlgorithmType algo,
                             matrix *predicted_y,
                             matrix *pred_residuals,
                             size_t nthreads,
                             ssignal *s,
                             int num_arg,
                             ...);

/*
 * Description: Calculate the Leave-One-Out cross validation.
 *              Every method has its model function. 
 * 
 *              Available methods:
 *              - PLS: partial least squares
 *              - PLS_DA: partial least squares discriminant analysis
 *              - EPLS: ensemble partial least squares
 *              - EPLS_DA: ensemble partial least squares discriminant analysis
 *              - MLR: multiple linear regression
 *              - LDA: linear discriminant analysis
 */
void LeaveOneOut(MODELINPUT *input,
                 AlgorithmType algo,
                 matrix *predicted_y,
                 matrix *pred_residuals,
                 size_t nthreads,
                 ssignal *s,
                 int num_arg,
                 ...);

/*
 * Description: Calculate the k-fold cross validation.
 *              Every method has its model function. 
 * 
 *              Available methods:
 *              - PLS: partial least squares
 *              - PLS_DA: partial least squares discriminant analysis
 *              - EPLS: ensemble partial least squares
 *              - EPLS_DA: ensemble partial least squares discriminant analysis
 *              - MLR: multiple linear regression
 *              - LDA: linear discriminant analysis
 */
void KFoldCV(MODELINPUT *input,
             uivector *groups,
             AlgorithmType algo,
             matrix *predicted_y,
             matrix *pred_residuals,
             size_t nthreads,
             ssignal *s,
             int arg,
             ...);

/*
 * Description: Calculate the Y-scrambling cross validation test.
 *              Every method has its model function. 
 * 
 *              Available methods:
 *              - PLS: partial least squares
 *              - PLS_DA: partial least squares discriminant analysis
 *              - EPLS: ensemble partial least squares
 *              - EPLS_DA: ensemble partial least squares discriminant analysis
 *              - MLR: multiple linear regression
 *              - LDA: linear discriminant analysis
 */
void YScrambling(MODELINPUT *input,
                 AlgorithmType algo,
                 ValidationArg varg,
                 size_t iterations,
                 matrix *ccoeff_yscrambling,
                 size_t nthreads,
                 ssignal *s);

#endif

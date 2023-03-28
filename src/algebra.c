/* algebra.c
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

#include <stdio.h>
#include <stdlib.h>

#include "numeric.h"
#include "algebra.h"
#include "matrix.h"
#include "vector.h"

/* Solve Linear System of Equation with the Gauss-Jordan
 */
void SolveLSE(matrix *mx, dvector *solution)
{
  size_t i;
  size_t j;
  size_t k;
  long int l;
  double tmp;
  matrix *X;
  initMatrix(&X);
  MatrixCopy(mx, &X);
  /*
   Problem: Solve the following system:
      x + y + z  = 4
      x - 2y - z = 1
      2x - y - 2z = -1

    Start out by multiplying the
    first X by -1 and add
    it to the second X to
    eliminate x from the second
    (*X).

    -x  - y - z = -4
      x - 2y - z = 1
    ----------------
        -3y - 2z = -3

    Now eliminate x from the third
    X by multiplying the first
    X by -2 and add it to
    the third (*X).

    -2x - 2y - 2z = -8
      2x -  y - 2z = -1
    ------------------
          -3y - 4z = -9

    Next, eliminate y from the third
    X by multiplying the second
    X by -1 and adding it to
    the third (*X).

      3y +  2z = 3
    -3y -  4z = -9
    --------------
          -2z = -6

    Solve the third X for z.

    -2z = -6
      z = 3

    Substitute 3 for z in the
    second X and solve for y.

    -3y - 2z = -3
    -3y - 2(3) = -3
    -3y - 6 = -3
    -3y = 3
    y = -1

    Lastly, substitute -1 for y and
    3 for z in the first X
    and solve for x.

    x + (-1) + 3 = 4
    x + 2 = 4
    x = 2

    The answer is (2, -1, 3).

  */

  /* Reorganize the matrix in order to find the pivot != 0 in the first k row k column */
  for(k = 0; k < X->row; k++){
    if(FLOAT_EQ(X->data[k][k], 0, 1e-4)){
      for(i = 0; i < X->row; i++){
        if(FLOAT_EQ(X->data[i][k], 0, 1e-4) == 0){
          /* move the row i to the null value */
          for(j = 0; j < X->col; j++){
            tmp = X->data[i][j];
            X->data[i][j] = X->data[k][j];
            X->data[k][j] = tmp;
          }
          break;
        }
        else{
          continue;
        }
      }
    }
  }

  /* (*X).row is the number of X, so is equal to the number of unknowns variables */
  for(k = 0; k < X->row; k++){
    for(i = k+1; i < X->row; i++){
      if(FLOAT_EQ(X->data[i][k], 0, 1e-4) == 0){ /* if the value is not 0 */
        if(FLOAT_EQ(X->data[k][k], 0, 1e-4) == 1){
          tmp = 0.f;
        }
        else{
          tmp = X->data[i][k]/X->data[k][k];
        }

        for(j = 0; j < X->col; j++){
          X->data[i][j] = (X->data[k][j] * (-tmp)) + X->data[i][j];
        }
      }
      else
        continue;
    }
  }

  /*Reset solution vector*/
  if(solution->size != X->row){
   DVectorResize(solution, X->row);
  }

  l = (*X).row-1;

  while(l > -1){
    double b = 0.f;
    for(i = 0; i < (*X).col-1; i++){
      if(i != l){
        b += X->data[l][i] * solution->data[i];
      }
      else
        continue;
    }

    if(FLOAT_EQ(X->data[l][l], 0, 1e-4) == 1)
      solution->data[l] = 0.f;
    else
      solution->data[l] = (X->data[l][X->col-1] -b) / X->data[l][l];
    l--;
  }
  DelMatrix(&X);
}



/*
 * Starting from the point that n equations are stocked together and written in vector form as_
 *
 *               y = Xβ + ε
 *
 * where:
 * - y is called the regressand, endogenous variable, response variable, measured variable, or dependent variable.
 * - X is a matrix composed by x_i regressors, exogenous variables, explanatory variables, covariates, input variables, predictor variables, or independent variables.
 * - β is a p-dimensional parameter vector. Its elements are also called effects, or regression coefficients. Statistical estimation and inference in linear regression focuses on β.
 * - is called the error term, disturbance term, or noise. This variable captures all other factors which influence the dependent variable yi other than the regressors xi.
 * -ε is a matrix composed by ε_i parameters that are the relationship between the error term and the regressors, for example whether they are correlated,
 * is a crucial step in formulating a linear regression model, as it will determine the method to use for estimation.
 *
 *
 * Ordinary least squares (OLS) is the simplest and thus most common estimator. It is conceptually simple and computationally straightforward.
 * OLS estimates are commonly used to analyze both experimental and observational data.
 * The OLS method minimizes the sum of squared residuals, and leads to a closed-form expression for the estimated value of the unknown parameter β:
 *
 *              β = (X' * X)^-1 * X'* y
 *
 * The estimator is unbiased and consistent if the errors have finite variance and are uncorrelated with the regressors
 *
 *              E[xiεi] = 0.
 *
 * It is also efficient under the assumption that the errors have finite variance and are homoscedastic, meaning that E[εi2|xi] does not depend on i.
 * The condition that the errors are uncorrelated with the regressors will generally be satisfied in an experiment, but in the case of observational data,
 * it is difficult to exclude the possibility of an omitted covariate z that is related to both the observed covariates and the response variable.
 * The existence of such a covariate will generally lead to a correlation between the regressors and the response variable, and hence to an inconsistent estimator of β.
 * The condition of homoscedasticity can fail with either experimental or observational data.
 * If the goal is either inference or predictive modeling, the performance of OLS estimates can be poor if multicollinearity is present, unless the sample size is large.
 * In simple linear regression, where there is only one regressor (with a constant), the OLS coefficient estimates have a simple form that is closely related to
 * the correlation coefficient between the covariate and the response.
 *
 *
 *
 * In order to find the coefficient value "b" of the polynomial X, the procedure applied is this:
 *
 * b = (Z' * Z')^-1 * Z'* y
 *
 * were :
 * b = vector of coefficients (b0, b1, b2, b...., bn)
 * Z' = is the sign matrix generated for the Centro Simmetric CCD Design
 * y = Dependent value vector. For each experiment we have a dependent value.
 *
 * The standard errors of the varuous estimates are the square roots of the corresponding diagoknal elementso of (Z' * Z')^-1 σ^2
 *
 * For info check the P.Box - Empirical Model-Building and response surfaces  Pag. 304-316
 */

/* x are the indipendent value of the polynomial equation
   y are the dependent value
 */
void OrdinaryLeastSquares(matrix *x, dvector *y, dvector *coefficients)
{
  /*y->size = x->row
   *
   * Z' = n  x m => n = x->col; m = x->row
   *
   * Z'Z = m x n * n x m  =  m x m  => m = x->row
   *
   * Z'y =
   *
   *
   */
  matrix *x_t; /* Z' that is x transposed */
  matrix *x_x_t; /* Z'*Z is the product between x and x_t */
  matrix *x_x_t_i; /* x_x_t_i inverted matrix  (Z' * Z')^-1 */
  dvector *x_t_y; /* Z'* y  = x->col */

  /* transpose x to x_t*/
  NewMatrix(&x_t, x->col, x->row);
  MatrixTranspose(x, x_t);

  /* Z'*Z */
  NewMatrix(&x_x_t, x->col, x->col);
  MatrixDotProduct(x_t, x, x_x_t);

  /* (Z'*Z)^-1 */
  initMatrix(&x_x_t_i);
  MatrixInversion(x_x_t, x_x_t_i);

  /* create a vector for store the results of Z'y */
  NewDVector(&x_t_y, x->col);

  /* Z'y */
  MatrixDVectorDotProduct(x_t, y, x_t_y);

  /* final data coefficient value */
  DVectorResize(coefficients, x->col);
  MatrixDVectorDotProduct(x_x_t_i, x_t_y, coefficients);

  DelDVector(&x_t_y);
  DelMatrix(&x_x_t_i);
  DelMatrix(&x_x_t);
  DelMatrix(&x_t);
}

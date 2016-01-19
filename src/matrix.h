/* matrix.h
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

#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>

#include "vector.h"

typedef struct{
  double **data;
  size_t row, col;
}matrix;

void initMatrix(matrix **m);
void NewMatrix(matrix**, size_t row_, size_t col_);
void ResizeMatrix(matrix **m, size_t row_, size_t col_); /*ResizeMatrix Delete all the data stored inside the matrix */
void DelMatrix(matrix**);

/*check for nan != finite */
void MatrixCheck(matrix *m);
void FindNan(matrix *m);
void PrintMatrix(matrix *m);
int ValInMatrix(matrix *m, double val);
void MatrixSet(matrix *m, double val);
void MatrixCopy(matrix *msrc, matrix **mdst);
void setMatrixValue(matrix *m, size_t row, size_t col, double val);
double getMatrixValue(matrix *m, size_t row, size_t col);
dvector *getMatrixRow(matrix *m, size_t row);
dvector *getMatrixColumn(matrix *m, size_t col);
void MatrixAppendRow(matrix **mx, dvector *row);
void MatrixAppendCol(matrix **mx, dvector *col);

/*  Matrix Operations */

void MatrixDVectorDotProduct(matrix *m, dvector *v, dvector *r);
/*The standard serial is much faster OBSOLETE*/
void MultiThreadMatrixDVectorDotProduct(matrix *mx, dvector *v, dvector *p, size_t nthreads_);
void DVectorMatrixDotProduct(matrix *mx, dvector *v, dvector *p);
/* m  = v1'v2 */
void DVectorTrasposedDVectorDotProduct(dvector *v1, dvector *v2, matrix *m);
/* r = v/mx  = (inv(mx^T)*v^T)^T*/
void DVectorTransposedMatrixDivision(dvector *v, matrix *mx, dvector *r);
void MatrixDotProduct(matrix *m_t, matrix *m, matrix *r);
void RowColOuterProduct(dvector *a, dvector *b, matrix *m);
void MatrixTranspose(matrix *m, matrix *r);
/* Matrix inversion using the LU decomposition */
void MatrixInversion2(matrix *m, matrix **m_inv);
/* Matrix inversion using the Gauss-Jordan algorithm */
void MatrixInversion(matrix *m, matrix **m_inv);
void MatrixPseudoinversion(matrix *m, matrix **m_inv);
void GenIdentityMatrix(matrix **m);
void MeanCenteredMatrix(matrix *mx, matrix *mxc);
void PearsonCorrelMatrix(matrix *mxsrc, matrix *mxdst);
void SpearmanCorrelMatrix(matrix *mxsrc, matrix *mxdst);

/*Calculate the column Average for the matrix mx*/
void MatrixColAverage(matrix *mx, dvector **colaverage);

/*Calculate the row Average for the matrix mx*/
void MatrixRowAverage(matrix *mx, dvector **rowaverage);

/*Calculate the column Standard Deviation for the matrix mx*/
void MatrixColSDEV(matrix *mx, dvector **colsdev);

/*Calculate the column Root Mean Square for the matrix mx*/
void MatrixColRMS(matrix* mx, dvector** colrms);

/*Calculate the column variance for the matrix mx*/
void MatrixColVar(matrix *mx, dvector **colvar);

/*Calculate the covariance matrix*/
void MatrixCovariance(matrix *mx, matrix **cm);

/* Transform a matrix into a logaritmic matrix */
void Matrix2LogMatrix(matrix *mx_in, matrix **mx_out);

/* Transform a matrix into a SQUARE matrix */
void Matrix2SquareMatrix(matrix *mx_in, matrix **mx_out);

/* Transform a matrix into a SQRT matrix */
void Matrix2SQRTMatrix(matrix *mx_in, matrix **mx_out);

/* Transform a matrix into ABS matrix */
void Matrix2ABSMatrix(matrix *mx_in, matrix **mx_out);

/* Develop an interaction factors matrix
 * Es. Use in DOE
 */
void Matrix2IntFactorsMatrix(matrix *mx_in, size_t factors, matrix **mx_out);

/* Transform a matrix into a row centered scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixRowCenterScaling(matrix *mx_in, matrix **mx_out);

/* Transform a matrix into a SVN row scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixSVNScaling(matrix *mx_in, matrix **mx_out);


double Matrixnorm(matrix *mx);
double Matrix1norm(matrix *mx);
void MatrixNorm(matrix *mx, matrix *nmx);
void MatrixColumnMinMax(matrix* mx, size_t col, double* min, double* max);
void MatrixSort(matrix *mx, size_t col_n);
void MatrixGetMaxValue(matrix *mx, size_t *row, size_t *col);
void MatrixGetMinValue(matrix *mx, size_t *row, size_t *col);

/* Singular Value Decomposition
*/
void SVD(matrix* mx, matrix **U, matrix **S, matrix **VT);
void SVDlapack(matrix *mx, matrix **u, matrix **s, matrix **vt);

/* Eigenvectors and Eigenvalues with the QR Method */
void EVectEval(matrix *mx, dvector **eval, matrix **evect);

/* QR Decomposition with the householder method */
void QRDecomposition(matrix *mx, matrix **Q, matrix **R);
void LUDecomposition(matrix *mx, matrix **L, matrix **U);

void HouseReflectorVect(dvector *x, dvector *u);
void HouseholderVectMatrixProduct(dvector *h, matrix *A, matrix *P);
void HouseholderMatrix(dvector *v, matrix *h);

/*These function reduce a matri to an hessenberg form*/
void HouseholderReduction(matrix *mx);
void CholeskyReduction(matrix *m);

#endif

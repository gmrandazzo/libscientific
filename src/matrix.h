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

/* Description:
 * Matrix data structure
 */
typedef struct{
  double **data;
  size_t row, col;
}matrix;

/* Description:
 * Allocate a matrix in memory without a dimension.
 */
void initMatrix(matrix **m);

/* Description:
 * Allocate a matrix in memory with dimension "row" and "col"
 */
void NewMatrix(matrix**, size_t row_, size_t col_);

/* Description:
 * Resize a matrix to new "row" and "col"
 */
void ResizeMatrix(matrix **m, size_t row_, size_t col_); /*ResizeMatrix Delete all the data stored inside the matrix */

/* Description:
 * Destroy a matrix allocation
 */
void DelMatrix(matrix**);

/* Description:
 * Check if some value in matrix are NAN or INFINITE
 */
void MatrixCheck(matrix *m);

/* Description:
 * Find NAN and print the output position to cmdout
 */
void FindNan(matrix *m);

/* Description:
 * Print a matrix to cmdout
 */
void PrintMatrix(matrix *m);

/* Description:
 * Check if a value "val" is in the matrix.
 * if the value is present return 1 else 0.
 */
int ValInMatrix(matrix *m, double val);

/* Description:
 * Set all values of a matrix to "val"
 */
void MatrixSet(matrix *m, double val);


/* Description:
 * Initialize the matrix elements to a random integer
 */
void MatrixInitRandomInt(matrix *m, int low, int high);


/* Description:
 * Initialize the matrix elements to a random float
 */
void MatrixInitRandomFloat(matrix *m, double low, double high);

/* Description:
 * Copy a matrix from a source "msrc" to a destination matrix "mdst"
 */
void MatrixCopy(matrix *msrc, matrix **mdst);

/* Description:
 * Set  matrix value "value" at position "row" and "col"
 */
void setMatrixValue(matrix *m, size_t row, size_t col, double val);

/* Description:
 * Get the matrix value at position "row" and "col"
 */
double getMatrixValue(matrix *m, size_t row, size_t col);

/* Description:
 * Get row matrix with index "row" as double vector.
 */
dvector *getMatrixRow(matrix *m, size_t row);

/* Description:
 * Get column matrix with index "col" as double vector.
 */
dvector *getMatrixColumn(matrix *m, size_t col);

/* Description:
 * Append double vector uivector as row.
 */
void MatrixAppendRow(matrix **mx, dvector *row);

/* Description:
 * Append double vector uivector as column.
 */
void MatrixAppendCol(matrix **mx, dvector *col);
/* Description:
 * Append unsigned int vector uivector as row
 */
void MatrixAppendUIRow(matrix **mx, uivector *row);

/* Description:
 * Append unsigned int vector uivector as column
 */
void MatrixAppendUIRow(matrix **mx, uivector *row);

/*  Matrix Operations */

/* Description:
 * matrix - row double vector product: the result is a row double vector
 * i.e.: X(10x5) * d(5x1) = r(10x1)
 */
void MatrixDVectorDotProduct(matrix *m, dvector *v, dvector *r);

/* Multithread version of MatrixDVectorDotProduct */
void MT_MatrixDVectorDotProduct(matrix *mx, dvector *v, dvector *p);

/* Description:
 * column double vector - matrix product: the result is a column double vector
 * i.e.: d(1x5) * X(5x10) = r(1x10)
 */
void DVectorMatrixDotProduct(matrix *mx, dvector *v, dvector *p);

/* Multithread version of DVectorMatrixDotProduct */
void MT_DVectorMatrixDotProduct(matrix *mx, dvector *v, dvector *p);

/* Description:
 * transposed double vector - double vecrtor product: the result is a matrix
 * i.e.: X  = d'd
 * i.e.: d'(5x1) * d(1x5) = X(1x10)
 */
void DVectorTrasposedDVectorDotProduct(dvector *v1, dvector *v2, matrix *m);

/* r = v/mx  = (inv(mx^T)*v^T)^T*/
void DVectorTransposedMatrixDivision(dvector *v, matrix *mx, dvector *r);

/*
 * Description:
 * Calculate the matrix matrix product
 */
void MatrixDotProduct(matrix *m_t, matrix *m, matrix *r);
void RowColOuterProduct(dvector *a, dvector *b, matrix *m);

/*
 * Description:
 * Generate a transpose matrix of m
 */
void MatrixTranspose(matrix *m, matrix *r);

/*
 * Description:
 * Matrix inversion using the Gauss-Jordan algorithm
 */
void MatrixInversion(matrix *m, matrix **m_inv);

/*
 * Description:
 * Matrix inversion using the LU decompositio
 */
void MatrixLUInversion(matrix *m, matrix **m_inv);

/*
 * Description:
 * Matrix pseudo inversion using the SVD algorithm
 */
void MatrixPseudoinversion(matrix *m, matrix **m_inv);

/*
 * Description:
 * Matrix pseudo inversion using the Moore-Penrose pseudoinverse
 */
void MatrixMoorePenrosePseudoinverse(matrix *m, matrix **inv);

/*
 * Description:
 * Generate the identity matrix
 */
void GenIdentityMatrix(matrix **m);

/*
 * Description:
 * Calculate the mean centered matrix
 */
void MeanCenteredMatrix(matrix *mx, matrix *mxc);

/*
 * Description:
 * Calculate the pearson correlation matrix
 */
void PearsonCorrelMatrix(matrix *mxsrc, matrix *mxdst);

/*
 * Description:
 * Calculate the spearmann correlation matrix
 */
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

/* Calculate the matrix descriptive statistics:
 *  - Column Average
 *  - Column Median
 *  - Column Armonic Average
 *  - Column Variance Population
 *  - Column Variance Sample (Correcter Variance)
 *  - Column Standard Deviation
 *  - Column Standard Deviation Sample (Corrected Standard Deviation)
 *  - Column Max
 *  - Column Min
 *  - Column of missing values
 */
void MatrixColDescStat(matrix *mx, matrix **ds);

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

/*
 * Description:
 * Develop an interaction factors matrix
 * Es. Use in DOE
 */
void Matrix2IntFactorsMatrix(matrix *mx_in, size_t factors, matrix **mx_out);

/*
 * Description:
 * Transform a matrix into a row centered scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixRowCenterScaling(matrix *mx_in, matrix **mx_out);

/*
 * Description:
 * Transform a matrix into a SVN row scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixSVNScaling(matrix *mx_in, matrix **mx_out);

/*
 * Description:
 * calculate the square root of the sum of the squares
 * of all the elements in the matrix
 * ||X|| = double
 */
double Matrixnorm(matrix *mx);
double Matrix1norm(matrix *mx);

/*
 * Description:
 * calculate the determinant of a matrix
 */
double MatrixDeterminant(matrix *mx);

/*
 * Description:
 * Normalize the matrix for the Matrixnorm value.
 * Each value of mx is divided by double Matrixnorm(matrix *mx);
 */
void MatrixNorm(matrix *mx, matrix *nmx);

/*
 * Description:
 * Find the minimum and maximum of a column in matrix
 */
void MatrixColumnMinMax(matrix* mx, size_t col, double* min, double* max);

/*
 * Description:
 * Sort of a matrix ba a column number col_n
 */
void MatrixSort(matrix *mx, size_t col_n);

/*
 * Description:
 * Reverse sort of a matrix by a column number col_n
 */
void MatrixReverseSort(matrix* mx, size_t col_n);

/*
 * Description:
 * find the maximum value in matrix and return the row and col indexes
 */
void MatrixGetMaxValueIndex(matrix *mx, size_t *row, size_t *col);

/*
 * Description:
 * find the minimum value in matrix and return the row and col indexes
 */
void MatrixGetMinValueIndex(matrix *mx, size_t *row, size_t *col);

/*
 * Description:
 * Singular Value Decomposition local implementation
 */
void SVD(matrix* mx, matrix **U, matrix **S, matrix **VT);

/*
 * Description:
 * Singular Value Decomposition lapack implementation
 */
void SVDlapack(matrix *mx, matrix **u, matrix **s, matrix **vt);

/*
 * Description:
 * Eigenvectors and Eigenvalues with the QR Method
 */
void EVectEval(matrix *mx, dvector **eval, matrix **evect);

/*
 * Description:
 * QR Decomposition with the householder method
 */
void QRDecomposition(matrix *mx, matrix **Q, matrix **R);

/*
 * Description:
 * LU Decomposition with the householder method
 */
void LUDecomposition(matrix *mx, matrix **L, matrix **U);

/*
 * Description:
 * Householder Reflector vector
 */
void HouseReflectorVect(dvector *x, dvector *u);

/*
 * Description:
 * Householder Vector Matrix Multiplication
 */
void HouseholderVectMatrixProduct(dvector *h, matrix *A, matrix *P);

/*
 * Description:
 * Householder matrix
 */
void HouseholderMatrix(dvector *v, matrix *h);

/*
 * Description:
 *  Reduce a matrix to an Hessenberg form
 * - Householder implementation
 * - Cholesky implementation
 */
void HouseholderReduction(matrix *mx);
void CholeskyReduction(matrix *m);

#endif

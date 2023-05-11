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

/**
 * Matrix data structure 
 *
 * - **data** two dimensional array of double
 * - **row** number of rows
 * - **col** number of columns
 */
typedef struct{
  double **data;
  size_t row, col;
}matrix;

/**
 * Description:
 * Allocate a matrix in memory without a dimension.
 */
void initMatrix(matrix **m);

/**
 * Allocate a matrix in memory with dimension "row" and "col"
 * 
 * @param [in] m matrix data structure
 * @param [in] row number of rows
 * @param [in] col number of columns
 */
void NewMatrix(matrix** m, size_t row_, size_t col_);

/**
 * Resize a matrix to new "row" and "col"
 */
void ResizeMatrix(matrix *m, size_t row_, size_t col_); /*ResizeMatrix Delete all the data stored inside the matrix */

/**
 * Destroy a matrix allocation
 */
void DelMatrix(matrix**);

/**
 * Check if some value in matrix are NAN or INFINITE
 */
void MatrixCheck(matrix *m);

/**
 * Find NAN and print the output position to cmdout
 */
void FindNan(matrix *m);

/**
 * Print a matrix to cmdout
 */
void PrintMatrix(matrix *m);

/**
 * Check if a value "val" is in the matrix.
 * if the value is present return 1 else 0.
 */
int ValInMatrix(matrix *m, double val);

/**
 * Set all values of a matrix to "val"
 */
void MatrixSet(matrix *m, double val);


/**
 * Initialize the matrix elements to a random integer
 */
void MatrixInitRandomInt(matrix *m, int low, int high);


/**
 * Initialize the matrix elements to a random float
 */
void MatrixInitRandomFloat(matrix *m, double low, double high);

/**
 * Copy a matrix from a source "msrc" to a destination matrix "mdst"
 */
void MatrixCopy(matrix *msrc, matrix **mdst);

/**
 * Set  matrix value "value" at position "row" and "col"
 */
void setMatrixValue(matrix *m, size_t row, size_t col, double val);

/**
 * Get the matrix value at position "row" and "col"
 */
double getMatrixValue(matrix *m, size_t row, size_t col);

/**
 * Get row matrix with index "row" as double vector.
 */
dvector *getMatrixRow(matrix *m, size_t row);

/**
 * Get column matrix with index "col" as double vector.
 */
dvector *getMatrixColumn(matrix *m, size_t col);

/**
 * Append double vector uivector as row.
 */
void MatrixAppendRow(matrix *m, dvector *row);

/**
 * Append double vector uivector as column.
 */
void MatrixAppendCol(matrix *m, dvector *col);

/**
 * Append unsigned int vector uivector as row
 */
void MatrixAppendUIRow(matrix *m, uivector *row);

/**
 * Append unsigned int vector uivector as column
 */
void MatrixAppendUIRow(matrix *m, uivector *row);

/**
 * Delete a specific row in a matrix
 */
void MatrixDeleteRowAt(matrix *m, size_t row);

/**
 * Delete a specific column in a matrix
 */
void MatrixDeleteColAt(matrix *m, size_t col);

/*  Matrix Operations */

/**
 * matrix - row double vector product: the result is a row double vector
 * i.e.: X(10x5) * d(5x1) = r(10x1)
 */
void MatrixDVectorDotProduct(matrix *m, dvector *v, dvector *p);

/** Multithread version of MatrixDVectorDotProduct */
void MT_MatrixDVectorDotProduct(matrix *m, dvector *v, dvector *p);

/**
 * column double vector - matrix product: the result is a column double vector
 * i.e.: d(1x5) * X(5x10) = r(1x10)
 */
void DVectorMatrixDotProduct(matrix *m, dvector *v, dvector *p);

/**
 * Multithread version of DVectorMatrixDotProduct
 */
void MT_DVectorMatrixDotProduct(matrix *m, dvector *v, dvector *p);

/**
 * transposed double vector - double vecrtor product: the result is a matrix
 * i.e.: X  = d'd
 * i.e.: d'(5x1) * d(1x5) = X(1x10)
 */
void DVectorTrasposedDVectorDotProduct(dvector *v1, dvector *v2, matrix *m);

/**
 * Calculate the vector multiplied with a transposed matrix division
 * according the following formula
 * r = v/m  = (inv(m^T)*v^T)^T
 */
void DVectorTransposedMatrixDivision(dvector *v, matrix *m, dvector *r);

/**
 * Calculate the matrix matrix product using the loop unrolling technique
 */
void MatrixDotProduct(matrix *a, matrix *b, matrix *r);
void MatrixDotProduct_(matrix *a, matrix *b, matrix *r);
void MatrixDotProduct_LOOP_UNROLLING(matrix *a, matrix *b, matrix *r);
void RowColOuterProduct(dvector *a, dvector *b, matrix *m);

/**
 * Generate a transpose matrix of m
 */
void MatrixTranspose(matrix *m, matrix *r);

/**
 * Matrix inversion using the Gauss-Jordan algorithm
 */
void MatrixInversion(matrix *m, matrix *m_inv);

/**
 * Matrix inversion using the LU decompositio
 */
void MatrixLUInversion(matrix *m, matrix *m_inv);

/**
 * Matrix pseudo inversion using the SVD algorithm
 */
void MatrixPseudoinversion(matrix *m, matrix *m_inv);

/**
 * Matrix pseudo inversion using the Moore-Penrose pseudoinverse
 */
void MatrixMoorePenrosePseudoinverse(matrix *m, matrix *inv);

/**
 * Return the trace of a square matrix
 */
double MatrixTrace(matrix *m);

/**
 * Generate the identity matrix
 */
void GenIdentityMatrix(matrix *m);

/**
 * Calculate the mean centered matrix
 */
void MeanCenteredMatrix(matrix *m, matrix *mc);

/**
 * Calculate the pearson correlation matrix
 */
void PearsonCorrelMatrix(matrix *msrc, matrix *mdst);

/**
 * Calculate the spearmann correlation matrix
 */
void SpearmanCorrelMatrix(matrix *msrc, matrix *mdst);

/*
 * Calculate the column Average for the matrix m
 */
void MatrixColAverage(matrix *m, dvector *colaverage);

/*
 * Calculate the row Average for the matrix m
 */
void MatrixRowAverage(matrix *m, dvector *rowaverage);

/*
 * Calculate the column Standard Deviation for the matrix m
 */
void MatrixColSDEV(matrix *m, dvector *colsdev);

/*
 * Calculate the column Root Mean Square for the matrix m
 */
void MatrixColRMS(matrix* m, dvector *colrms);

/*
 * Calculate the column variance for the matrix m
 */
void MatrixColVar(matrix *m, dvector *colvar);

/**
 *  Calculate the matrix descriptive statistics:
 * 
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
void MatrixColDescStat(matrix *m, matrix *ds);

/**
 * Calculate the covariance matrix
 */
void MatrixCovariance(matrix *m, matrix *cm);

/**
 * Transform a matrix into a logaritmic matrix
 */
void Matrix2LogMatrix(matrix *m_in, matrix *m_out);

/**
 * Transform a matrix into a SQUARE matrix
 */
void Matrix2SquareMatrix(matrix *m_in, matrix *m_out);

/**
 * Transform a matrix into a SQRT matrix
 */
void Matrix2SQRTMatrix(matrix *m_in, matrix *m_out);

/**
 * Transform a matrix into ABS matrix
 */
void Matrix2ABSMatrix(matrix *m_in, matrix *m_out);

/**
 * Develop an interaction factors matrix for DOE
 */
void Matrix2IntFactorsMatrix(matrix *m_in, size_t factors, matrix *m_out);

/**
 * Transform a matrix into a row centered scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixRowCenterScaling(matrix *m_in, matrix *m_out);

/**
 * Transform a matrix into a SVN row scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixSVNScaling(matrix *m_in, matrix *m_out);

/**
 * calculate the square root of the sum of the squares
 * of all the elements in the matrix
 * ||X|| = double
 */
double Matrixnorm(matrix *m);

/**
 * calculate the determinant of a matrix
 */
double MatrixDeterminant(matrix *m);

/**
 * Normalize the matrix for the Matrixnorm value.
 * Each value of m is divided by double Matrixnorm(matrix *m);
 */
void MatrixNorm(matrix *m, matrix *nm);

/**
 * Find the minimum and maximum of a column in matrix
 */
void MatrixColumnMinMax(matrix* m, size_t col, double* min, double* max);

/**
 * Sort of a matrix ba a column number col_n
 */
void MatrixSort(matrix *m, size_t col_n);

/**
 * Reverse sort of a matrix by a column number col_n
 */
void MatrixReverseSort(matrix* m, size_t col_n);

/**
 * find the maximum value in matrix and return the row and col indexes
 */
void MatrixGetMaxValueIndex(matrix *m, size_t *row, size_t *col);

/**
 * find the minimum value in matrix and return the row and col indexes
 */
void MatrixGetMinValueIndex(matrix *m, size_t *row, size_t *col);

/**
 * Singular Value Decomposition local implementation
 */
void SVD(matrix* m, matrix *U, matrix *S, matrix *VT);

/**
 * Singular Value Decomposition lapack implementation
 */
void SVDlapack(matrix *m_, matrix *u, matrix *s, matrix *vt);

/**
 * Eigenvectors and Eigenvalues with the QR Method
 */
void EVectEval(matrix *m, dvector *eval, matrix *evect);

/**
 * QR Decomposition with the householder method
 */
void QRDecomposition(matrix *m, matrix *Q, matrix *R);

#endif

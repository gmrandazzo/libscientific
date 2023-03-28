/* tensor.h
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

#ifndef TENSOR_H
#define TENSOR_H

#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "vector.h"

/*  Giuseppe Marco Randazzo
 *  Sun Dec 25 20:54:11 CET 2011
 */


typedef struct{
  matrix **m;
  size_t order;
}tensor;

/*
 * Description:
 * initialize an empty tensor
 */
void initTensor(tensor **t);

/*
 * Description:
 * Create a new empty tensor with a predefined order.
 */
void NewTensor(tensor** t, size_t order_);

/*
 * Description:
 * Create into an empty tensor a matrix of "row x col" at position "order"
 */
void NewTensorMatrix(tensor *t, size_t order, size_t row, size_t col);

/*
 * Description:
 * Append to a tensor a matrix at the end of it.
 */
void AddTensorMatrix(tensor *t, size_t row, size_t col);

/*
 * Description:
 * Delete a tensor and free up the memory.
 */
void DelTensor(tensor**t);

/*
 * Description:
 * Print to video the tensor.
 */
void PrintTensor(tensor *t);

/*
 * Description:
 * Set a value into a tensor 
 */
void setTensorValue(tensor *t, size_t order, size_t row, size_t col, double val);

/*
 * Description:
 * Get a value from a tensor
 */
double getTensorValue(tensor *t, size_t order, size_t row, size_t col);

/*
 * Description:
 * Append a matrix to a tensor at the end of it.
 */
void TensorAppendMatrix(tensor *tdst, matrix *msrc);

/*
 * Description:
 * Append a matrix to a tensor in a specific position of it;
 */
void TensorAppendMatrixAt(tensor *tdst, size_t order, matrix *msrc);

/*
 * Description:
 * Append a column vector into the matrix of a tensor at position "order"
 */
void TensorAppendColumn(tensor *t, size_t order, dvector* column);

/*
 * Description:
 * Append a row vector into the matrix of a tensor at position "order"
 */
void TensorAppendRow(tensor *t, size_t order, dvector* row);

/*
 * Description:
 * Set all the values of a tensor to a specific value "val"
 */
void TensorSet(tensor *t, double val);

/*
 * Description:
 * Copy a tensor from a source into a destination
 */
void TensorCopy(tensor *asrc, tensor **adst);


/*
 * Description:
 * Calculate the mean centered tensor
 */
void MeanCenteredTensor(tensor *t, tensor *tc);

/*
 * Description:
 * Calculate column average of each matrix in a tensor
 */
void TensorColAverage(tensor *t, matrix *colaverage);

/*
 * Description:
 * Calculate column standard deviation of each matrix in a tensor
 */
void TensorColSDEV(tensor *t, matrix *colsdev);

/*
 * Description:
 * Transpose a tensor
 */
void TensorTranspose(tensor *t1, tensor *t2);

/*
 * Description:
 * Calculate the transposed tensor x vector product.
 * The result is a matrix.
 */
void TransposedTensorDVectorProduct(tensor *t, dvector *v, matrix *p);

/*
 * Description:
 * Calculate the vector x tensor product.
 * The result is a matrix.
 */
void DvectorTensorDotProduct(tensor *t, dvector *v, matrix *m);

/*
 * Description:
 * Calculate the tensor-matrix product
 * The result is a matrix.
 */
void TensorMatrixDotProduct(tensor *t, matrix *m, dvector *v);

/*
 * Description:
 * Calculate the Kroneker dot product between a vector and a matrix.
 * The result is a tensor.
 */
void KronekerProductVectorMatrix(dvector *v, matrix *m, tensor *t);

#endif

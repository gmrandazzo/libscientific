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

void initTensor(tensor **t);
void NewTensor(tensor** t, size_t order_);
void NewTensorMatrix(tensor *t, size_t order, size_t row, size_t col);
void AddTensorMatrix(tensor *t, size_t row, size_t col);

void DelTensor(tensor**t);

void PrintTensor(tensor *t);
void setTensorValue(tensor *t, size_t order, size_t row, size_t col, double val);
double getTensorValue(tensor *t, size_t order, size_t row, size_t col);
void TensorAppendMatrix(tensor *tdst, matrix *msrc);
void TensorAppendMatrixAt(tensor *tdst, size_t order, matrix *msrc);
void TensorAppendColumn(tensor *t, size_t order, dvector* column);
void TensorAppendRow(tensor *t, size_t order, dvector* row);
void TensorSet(tensor *t, double val);
void TensorCopy(tensor *asrc, tensor **adst);



/* matrix *getRowPlane(tensor *t, size_t row);
matrix  *getColumnPlane(tensor *t, size_t col);
void TensorAppendRowPlane(tensor **t, matrix *row);
void TensorAppendColPlane(tensor **t, matrix *col);
*/

/*  Matrix Operations */
void MeanCenteredTensor(tensor *t, tensor *tc);
void TensorColAverage(tensor *t, matrix *colaverage);
void TensorColSDEV(tensor *t, matrix *colsdev);
void TensorTranspose(tensor *t1, tensor *t2);

void TransposedTensorDVectorProduct(tensor *t, dvector *v, matrix *p);
void DvectorTensorDotProduct(tensor *t, dvector *v, matrix *m);
void TensorMatrixDotProduct(tensor *t, matrix *m, dvector *v);
void KronekerProductVectorMatrix(dvector *v, matrix *m, tensor *t);


/*
void MatrixDVectorDotProduct(matrix *m, dvector *v, dvector *r);
void DVectorMatrixDotProduct(matrix *mx, dvector *v, dvector *p);
void MatrixDotProduct(matrix *m_t, matrix *m, matrix **r);
void MatrixTranspose(matrix *m, matrix *r);
void MatrixInversion(matrix *A, matrix **Y);
void GenIdentityMatrix(matrix **m);
void MatrixColAverage(matrix *mx, dvector **colaverage);
void MatrixColSDEV(matrix *mx, dvector **colsdev);
*/
#endif

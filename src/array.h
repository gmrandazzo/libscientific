#ifndef ARRAY_H
#define ARRAY_H

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
}array;

void initArray(array **t);
void NewArray(array** t, size_t order_);
void NewArrayMatrix(array **t, size_t order, size_t row, size_t col);
void AddArrayMatrix(array **t, size_t row, size_t col);

void DelArray(array**t);

void PrintArray(array *t);
void setArrayValue(array *t, size_t order, size_t row, size_t col, double val);
double getArrayValue(array *t, size_t order, size_t row, size_t col);
void ArrayAppendMatrix(array **tdst, matrix *msrc);
void ArrayAppendMatrixAt(array **tdst, size_t order, matrix *msrc);
void ArrayAppendColumn(array **t, size_t order, dvector* column);
void ArrayAppendRow(array **t, size_t order, dvector* row);
void ArraySet(array *t, double val);
void ArrayCopy(array *asrc, array **adst);



/* matrix *getRowPlane(array *t, size_t row);
matrix  *getColumnPlane(array *t, size_t col);
void ArrayAppendRowPlane(array **t, matrix *row);
void ArrayAppendColPlane(array **t, matrix *col);
*/

/*  Matrix Operations */
void MeanCenteredArray(array *t, array *tc);
void ArrayColAverage(array *t, matrix **colaverage);
void ArrayColSDEV(array *t, matrix **colsdev);
void ArrayTranspose(array *t1, array *t2);

void TransposedArrayDVectorProduct(array *t, dvector *v, matrix *p);
void DvectorArrayDotProduct(array *t, dvector *v, matrix *m);
void ArrayMatrixDotProduct(array *t, matrix *m, dvector *v);
void KronekerProductVectorMatrix(dvector *v, matrix *m, array *t);


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

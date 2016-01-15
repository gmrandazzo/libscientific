#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "array.h"
#include "matrix.h"
#include "numeric.h"
#include "memwrapper.h"

/*  Giuseppe Marco Randazzo 
 *  Sun Dec 25 20:54:11 CET 2011
 * 
 * Multidimensional Matrices
 */

void initArray(array **t)
{
  (*t) = xmalloc(sizeof(array));
  (*t)->order = 0;
  (*t)->m = NULL;
}

void NewArray(array** t, size_t order_)
{
  size_t i;
  (*t) = xmalloc(sizeof(array));
  (*t)->order = order_;
  (*t)->m = xmalloc(sizeof(matrix*)*order_);
  
  for(i = 0; i < order_; i++){
    (*t)->m[i] = NULL;
  }
}

void NewArrayMatrix(array **t, size_t order, size_t row, size_t col)
{
  if((*t)->order != 0){
    if(order < (*t)->order){
      NewMatrix(&(*t)->m[order], row, col);
    }
    else{
      fprintf(stderr, "Error! Order not in range order. order: %u > ordersize: %u\n", (unsigned int)order, (unsigned int)(*t)->order);
    }
  }
  else{
    fprintf(stderr, "Error! Order Array not Defined! Please create Array with NewArray(array **t, size_t order);\n");
    abort();
  }
}

void AddArrayMatrix(array **t, size_t row, size_t col)
{
   if((*t)->order > 0){
    /*if((*t)->m[(*t)->order-1]->row != row){
      fprintf(stderr, "Error while appending matrix to array! Object size differ: matrix: %u;  array: %u", (unsigned int)row, (unsigned int)(*t)->m[(*t)->order-1]->row);
      abort();
    }
    else{*/
      (*t)->order += 1;
      (*t)->m = xrealloc((*t)->m, sizeof(matrix*)*(*t)->order);
      NewMatrix(&(*t)->m[(*t)->order-1], row, col);
/*     }*/
  }
  else{
    (*t)->order = 1;
    (*t)->m = xmalloc(sizeof(matrix*)*1);
    NewMatrix(&(*t)->m[0], row, col);
  } 
}

void DelArray(array** t)
{
  size_t i;
  for(i = 0; i < (*t)->order; i++)
    DelMatrix(&(*t)->m[i]);
  free((*t)->m);
  free((*t));
}


void PrintArray(array *t)
{
  size_t i, j, k;
  printf("Array - order: %u\n", (unsigned int)t->order);
  for(i = 0; i < t->order; i++){
    printf("Array No: %u of row: %u; col: %u\n", (unsigned int)i+1, (unsigned int)t->m[i]->row, (unsigned int)t->m[i]->col);
    for(j = 0; j < t->m[i]->row; j++){
      for(k = 0; k < t->m[i]->col; k++)
        printf("%8.3f\t", getMatrixValue(t->m[i], j, k));
      printf("\n");
    }
  }
}

void setArrayValue(array *t, size_t i, size_t j, size_t k, double val)
{
  if(i < (*t).order){
    if(j < (*t).m[i]->row){
      if(k < (*t).m[i]->col){
        setMatrixValue((*t).m[i], j, k, val);
      }
      else{
        fprintf(stderr, "setArrayValue Error! Wrong colum index.\n");
        fflush(stderr);
        abort();
      }
    }
    else{
      fprintf(stderr, "setArrayValue Error! Wrong row index.\n");
      fflush(stderr);
      abort();
    }
  }
  else{
    fprintf(stderr, "setArrayValue Error! Wrong order.\n");
    fflush(stderr);
    abort();
  }
}

double getArrayValue(array *t, size_t i, size_t j, size_t k){
    if(i < (*t).order){
    if(j < (*t).m[i]->row){
      if(k < (*t).m[i]->col){
        return getMatrixValue((*t).m[i], j, k);
      }
      else{
        fprintf(stderr, "getArrayValue Error! Wrong colum index.\n");
        fflush(stderr);
        return NAN;
      }
    }
    else{
      fprintf(stderr, "getArrayValue Error! Wrong row index.\n");
      fflush(stderr);
      return NAN;
    }
  }
  else{
    fprintf(stderr, "getArrayValue Error! Wrong order.\n");
    fflush(stderr);
    return NAN;
  }
}

void ArrayAppendMatrix(array **tdst, matrix *msrc)
{
  if((*tdst)->order > 0){
    if((*tdst)->m[(*tdst)->order-1]->row != msrc->row){
      fprintf(stderr, "Error while appending matrix to array! Object size differ: matrix: %u;  array: %u", (unsigned int)msrc->row, (unsigned int)(*tdst)->m[(*tdst)->order-1]->row);
      abort();
    }
    else{
      (*tdst)->order += 1;
      (*tdst)->m = xrealloc((*tdst)->m, sizeof(matrix*)*(*tdst)->order);
      NewMatrix(&(*tdst)->m[(*tdst)->order-1], msrc->row, msrc->col);
      MatrixCopy(msrc, &(*tdst)->m[(*tdst)->order-1]);
    }
  }
  else{
    (*tdst)->order = 1;
    (*tdst)->m = xmalloc(sizeof(matrix*)*1);
    NewMatrix(&(*tdst)->m[0], msrc->row, msrc->col);
    MatrixCopy(msrc, &(*tdst)->m[0]);
  }
}

void ArrayAppendMatrixAt(array **tdst, size_t order, matrix *msrc)
{
  if(order < (*tdst)->order){
    fprintf(stderr, "Module not developed. Work in progress...\n");
    fflush(stderr);
    abort();
    /*array *tmp;
    initArray(&tmp);
    
    
    DelArray(&tmp);*/
  }
  else{
    ArrayAppendMatrix(tdst, msrc);
  }
}

/*
 * t is the array
 * n is the order number
 * column is the colum vector to append
 */
void ArrayAppendColumn(array **t, size_t n, dvector* column)
{
  if(n < (*t)->order){
    MatrixAppendCol(&(*t)->m[n], column);
    /*
    if(column->size == (*t)->m[n]->row){
      MatrixAppendCol(&(*t)->m[n], column);
    }
    else{
      fprintf(stderr, "Error! The objects number differ %u != %u\n", (unsigned int)column->size, (unsigned int)(*t)->m[n]->row);
      fflush(stderr);
      abort();
    }
    */
  }
  else{
    if(n > (*t)->order){
      fprintf(stderr, "Error! Order number too high. %u > %u\n", (unsigned int)n, (unsigned int)(*t)->order );
    }
    else{
      fprintf(stderr, "Error! Order 0!\n");
    }
    fflush(stderr);
    abort();
  }
}

/*
 * t is the array
 * n is the order number
 * column is the colum vector to append
 */
void ArrayAppendRow(array **t, size_t n, dvector* row)
{
  if(n < (*t)->order){
    if(row->size != (*t)->m[n]->col){
      MatrixAppendRow(&(*t)->m[n], row);
    }
    else{
      fprintf(stderr, "Error! The column number differ %u != %u\n", (unsigned int)row->size, (unsigned int)(*t)->m[n]->col);
      fflush(stderr);
      abort();
    }
  }
  else{
    fprintf(stderr, "Error! Order number too high. %u > %u\n", (unsigned int)n, (unsigned int)(*t)->order);
    fflush(stderr);
    abort();
  }
}

void ArraySet(array* t, double val)
{
  size_t i;
  for(i = 0; i < t->order; i++)
    MatrixSet(t->m[i], val);
}

void ArrayCopy(array* asrc, array** adst)
{
  size_t i, j, k;
  if((*adst)->m == NULL){
    (*adst)->order = asrc->order;
    
    (*adst)->m = xmalloc(sizeof(matrix*)*asrc->order);
  
    for(k = 0; k < asrc->order; k++){
      NewMatrix(&((*adst)->m[k]), asrc->m[k]->row, asrc->m[k]->col);
    }
  }
  else{
    if(asrc->order != (*adst)->order){
      /* resize  the order */
      (*adst)->m = xrealloc((*adst)->m, sizeof(array*)*asrc->order);
    }

    /*chek and resize the matrix for each order if is necessary */
    for(k = 0; k < asrc->order; k++){
      if(asrc->m[k]->row != (*adst)->m[k]->row || asrc->m[k]->col != (*adst)->m[k]->col){
        
        (*adst)->m[k]->row = asrc->m[k]->row;
        (*adst)->m[k]->col = asrc->m[k]->col;

        (*adst)->m[k]->data = xrealloc((*adst)->m[k]->data, sizeof(double*)*asrc->m[k]->row);
        
        for(i = 0; i < asrc->m[k]->row; i++){
          (*adst)->m[k]->data[i] = xrealloc((*adst)->m[k]->data[i], sizeof(double)*asrc->m[k]->col);
        }
      }
    }
    
  }
  
  /*copy the data...*/
  for(k = 0; k < asrc->order; k++){
    for(i = 0; i < asrc->m[k]->row; i++){
      for(j = 0; j < asrc->m[k]->col; j++){
        setArrayValue((*adst), k, i, j, getArrayValue(asrc, k, i, j));
      }
    }
  }

}


/*
 * t is the input data array
 * tc is the output array meancentered. This must be created with NewArray()
 */
void MeanCenteredArray(array *t, array *tc)
{
  size_t i;
  for(i = 0; i < t->order; i++){
    MeanCenteredMatrix(t->m[i], tc->m[i]);
  }
}

void ArrayColAverage(array *t, matrix **colaverage)
{
  size_t i;
  dvector *v;
  for(i = 0; i < t->order; i++){
    initDVector(&v);
    MatrixColAverage(t->m[i], &v);
    MatrixAppendCol(colaverage, v);
    DelDVector(&v);
  }
}

void ArrayColSDEV(array *t, matrix **colsdev)
{
  size_t i;
  dvector *v;
  for(i = 0; i < t->order; i++){
    initDVector(&v);
    MatrixColSDEV(t->m[i], &v);
    MatrixAppendCol(colsdev, v);
    DelDVector(&v);
  }
}

void ArrayTranspose(array* t1, array* t2)
{
  size_t i;
  for(i = 0; i < t1->order; i++){
    MatrixTranspose(t1->m[i], t2->m[i]);
  }
}


/*
 * k = order
 * i = row
 * j = column
 * 
 * E'(j,i,k) = E(i,j,k)
 * 
 * P(k,j) = Sum_i[ E'(j,i,k)*t(i) ]
 * 
 */
void TransposedArrayDVectorProduct(array *t, dvector *v, matrix *p)
{
  size_t i, j, k;
  double res;
  for(k = 0; k < t->order; k++){
    for(i = 0; i < t->m[k]->row; i++){
      for(j = 0; j < t->m[k]->col; j++){
        res = getArrayValue(t, k, i, j)*getDVectorValue(v, j);
        if(_isnan_(res) || _isinf_(res)){
          setMatrixValue(p, k, i, getMatrixValue(p, k, i) + 0);
        }
        else{
          setMatrixValue(p, k, i, getMatrixValue(p, k, i) + res);
        }
      }
    }
  }
}


/* the output "m" is size of: 
 *   colum = t->order;
 *   row = v->size
 */
void DvectorArrayDotProduct(array* t, dvector* v, matrix* m)
{
  /*We do not need tests because getMatrixValue and setMatrixValue and getArrayValue handle some errors*/
  size_t i, j, k;
  double res;
  /* m = v't
  *
  * k = order of array
  * j = column size
  * i = row size
  * 
  * m[j][k] =   Σ v[i] * t[k][i][j]
  * 
  */
  for(k = 0; k < t->order; k++){
    for(i = 0; i < t->m[k]->row; i++){
      for(j = 0; j < t->m[k]->col; j++){
        if(v->size == t->m[k]->row && t->m[k]->col == m->row && t->order == m->col){
          res = getDVectorValue(v, i) * getArrayValue(t, k, i, j);
          if(_isnan_(res) || _isinf_(res)){
            setMatrixValue(m, j, k, (getMatrixValue(m, j, k) + 0));
          }
          else{
            setMatrixValue(m, j, k, (getMatrixValue(m, j, k) + res));
          }
        }
        else{
          fprintf(stderr, "Error while computing DvectorArrayDotProduct.\n");
          fflush(stderr);
          abort();
        }
      }
    }
  }
}

void ArrayMatrixDotProduct(array *t, matrix *m, dvector *v)
{
  size_t i, j, k;
  double res;
  for(k = 0; k < t->order; k++){
    if(m->col == t->order && m->row == t->m[k]->col){
      for(i = 0; i < t->m[k]->row; i++){
        for(j = 0; j < t->m[k]->col; j++){
          res = (getArrayValue(t, k, i, j)*getMatrixValue(m, j, k));
          if(_isnan_(res) || _isinf_(res)){
            setDVectorValue(v, i, (getDVectorValue(v, i) + 0));
          }
          else{
            setDVectorValue(v, i, (getDVectorValue(v, i) + res));
          }
        }
      }
    }
    else{
      fprintf(stderr, "Error while computing ArrayMatrixDotProduct.");
      fflush(stderr);
      abort();
    }
  }
}


/*
 * i = row
 * j = column
 * k = order
 * 
 * t(i) = Sum_j[ Sum_k[ E(i,j,k) * P(k,j)] ]
 * 
 */
void ArrayMatrixDotProduct2(array *t, matrix *m, dvector *v)
{
  size_t i, j, k;
  for(i = 0; i < v->size; i++){
    for(k = 0; k < t->order; k++){
      for(j = 0; j < t->m[k]->col; j++){
        setDVectorValue(v, i, getDVectorValue(v, i) + (getArrayValue(t, k, i, j)*getMatrixValue(m, k, j)));
      } 
    }
  }
}

void KronekerProductVectorMatrix(dvector* v, matrix* m, array* t)
{
  size_t i, j, k;
  double res;
  for(i = 0; i < v->size; i++){
    for(j = 0; j < m->row; j++){
      for(k = 0; k < m->col; k++){
        res = getDVectorValue(v, i)*getMatrixValue(m, j, k);
        if(_isnan_(res) || _isinf_(res)){
          setArrayValue(t, k, i, j, 0.f);
        }
        else{
          setArrayValue(t, k, i, j, res);
        }
      }
    }
  }
}

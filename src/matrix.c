/* matrix.c
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

#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "numeric.h"
#include "algebra.h"

#include "memwrapper.h"

#define MAXTHREADS 256

void initMatrix(matrix **m)
{
  (*m) = xmalloc(sizeof(matrix));
  (*m)->row = 0;
  (*m)->col = 0;
  (*m)->data = NULL;
}

void NewMatrix(matrix **m, size_t row_ , size_t col_)
{
  size_t i, j;
  (*m) = xmalloc(sizeof(matrix));
  (*m)->row = row_;
  (*m)->col = col_;
  (*m)->data = xmalloc(sizeof(double*)*row_);
  for(i = 0; i < row_; i++){
    (*m)->data[i] = xmalloc(sizeof(double)*col_);
    for(j = 0; j < col_; j++)
      (*m)->data[i][j] = +0.f;
  }
}

void ResizeMatrix(matrix *m, size_t row_, size_t col_)
{
  size_t i, j;
  if(m != NULL){
    if(m->row == row_ && m->col == col_){
      MatrixSet(m, +0.f);
    }
    else{
      if(m->col > 0 && m->row > 0){
        for(i = 0; i < m->row; i++){
          xfree(m->data[i]);
        }
        xfree(m->data);
      }

      m->data = xmalloc(sizeof(double*)*row_);
      for(i = 0; i < row_; i++){
        m->data[i] = xmalloc(sizeof(double)*col_);
        for(j = 0; j < col_; j++)
          m->data[i][j] = +0.f;
      }

      m->row = row_;
      m->col = col_;
    }
  }
  else{
    NewMatrix(&m, row_, col_);
  }
}

void DelMatrix(matrix **m)
{
  if(m != NULL){
    size_t i;
    for(i = 0; i < (*m)->row; i++)
      xfree((*m)->data[i]);
    xfree((*m)->data);
    xfree((*m));
  }
}

void MatrixCheck(matrix *m)
{
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(_isnan_(m->data[i][j]) || !isfinite(m->data[i][j])){
        m->data[i][j] = MISSING;
      }
    }
  }
}

void FindNan(matrix *m)
{
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(_isnan_(m->data[i][j])){
        printf("Nan at %d %d\n", (int)i, (int)j);
      }
    }
  }
}

void PrintMatrix(matrix *m)
{
  size_t i, j;
  printf("Matrix of row: %u; col: %u\n", (unsigned int)m->row, (unsigned int)m->col);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++)
      printf("%8.3f\t", m->data[i][j]);
    printf("\n");
  }
}

/*if a value is in matrix return 1 else 0*/
int ValInMatrix(matrix* m, double val)
{
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      /*if((size_t)getMatrixValue(m, i, j) == (size_t)id) */
      if(FLOAT_EQ(m->data[i][j], val, 1*10e-8))
        return 1;
      else
        continue;
    }
  }
  return 0;
}

void MatrixSet(matrix *m, double val)
{
  size_t i, j;
  if(m->row == m->col){
    for(i = 0; i < m->row; i++){
      m->data[i][i] = val;
      for(j = i+1; j < m->col; j++){
        m->data[i][j]= m->data[j][i] = val;
      }
    }
  }
  else{
    for(i = 0; i < m->row; i++){
      for(j = 0; j < m->col; j++){
        m->data[i][j] = val;
      }
    }
  }
}

void MatrixInitRandomInt(matrix *m, int low, int high)
{
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      m->data[i][j] = (double)randInt(low, high);
    }
  }
}

void MatrixInitRandomFloat(matrix *m, double low, double high)
{
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      m->data[i][j] = randDouble(low, high);
    }
  }
}

void MatrixCopy(matrix *msrc, matrix **mdst)
{
  size_t i, j;
  if((*mdst)->data == NULL){
    (*mdst)->row = msrc->row;
    (*mdst)->col = msrc->col;
    (*mdst)->data = xmalloc(sizeof(double*)*msrc->row);
    for(i = 0; i < msrc->row; i++){
      (*mdst)->data[i] = xmalloc(sizeof(double)*msrc->col);
      for(j = 0; j < msrc->col; j++)
        (*mdst)->data[i][j] = +0.f;
    }
  }
  else{
    if((*mdst)->row != msrc->row || (*mdst)->col != msrc->col){
      DelMatrix(mdst);
      initMatrix(mdst);
      (*mdst)->row = msrc->row;
      (*mdst)->col = msrc->col;

      (*mdst)->data = xmalloc(sizeof(double*)*msrc->row);
      for(i = 0; i < (*mdst)->row; i++){
        (*mdst)->data[i] = xmalloc(sizeof(double)*msrc->col);
        for(j = 0; j < msrc->col; j++){
          (*mdst)->data[i][j] = +0.f;
        }
      }
    }
  }

  for(i = 0; i < msrc->row; i++){
    for(j = 0; j < msrc->col; j++){
      (*mdst)->data[i][j] = msrc->data[i][j];
    }
  }
}

void setMatrixValue(matrix *m, size_t row, size_t col, double val)
{
  if(row < m->row && col < m->col){
    if(_isnan_(val) || _isinf_(val)){
      (*m).data[row][col] = MISSING;
    }
    else{
      (*m).data[row][col] = val;
    }
  }
  else{
    fprintf(stdout,"setMatrixValue Error: row to set: %u row max: %u; column to set %u; column max: %u out of range.\n", (unsigned int)row, (unsigned int)m->row-1, (unsigned int)col, (unsigned int)m->col-1);
    fflush(stdout);
//     abort();
  }
}

double getMatrixValue(matrix *m, size_t row, size_t col)
{
  if(row < (*m).row && col < (*m).col){
    return (*m).data[row][col];
  }
  else{
    fprintf(stdout,"getMatrixValue Error: row to get: %u row max: %u; column to get %u; column max: %u out of range.\n", (unsigned int)row, (unsigned int)m->row-1, (unsigned int)col, (unsigned int)m->col-1);
    fflush(stdout);
//     abort();
    return NAN;
  }
}

dvector *getMatrixRow(matrix *m, size_t row)
{
  if(row < (*m).row){
    dvector *v;
    size_t j;
    NewDVector(&v, m->col);
    for(j = 0; j < m->col; j++){
      v->data[j] = m->data[row][j];
    }

    /*
    initDVector(&v);
    memcpy(v->data, (*m).data[row], (*m).col);
    v->size = (*m).col;
    */
    return v;
  }
  else{
    fprintf(stdout,"getRow Error: row %u out of range.\n", (unsigned int)row);
    fflush(stdout);
    return NULL;
  }

}

dvector *getMatrixColumn(matrix *m, size_t col)
{
  if(col < (*m).col){
    dvector *v;
    size_t i;
    NewDVector(&v, m->row);
    for(i = 0; i < m->row; i++)
      v->data[i] = m->data[i][col];
    return v;
  }
  else{
    fprintf(stdout,"getColumn Error: column %u out of range.\n", (unsigned int)col);
    fflush(stdout);
    return NULL;
  }
}

void MatrixAppendRow(matrix* m, dvector *row)
{
  size_t i, j;
  size_t rowsize =  m->row + 1;
  size_t colsize;


  if(m->col != 0){
    if(row->size > m->col)
      colsize = row->size;
    else /*if (row->size <= m->col)*/
      colsize = m->col;
  }
  else{
    colsize = row->size;
  }

  /*adding a new row*/
  m->data = xrealloc(m->data, sizeof(double*)*rowsize);

  if(colsize > m->col){
    /*resize the column*/
    for(i = 0; i < m->row; i++){
      m->data[i] = xrealloc(m->data[i], sizeof(double*)*colsize);
      /*initialize the new column value*/
      for(j = m->col; j < colsize; j++){
        m->data[i][j] = +0.f;
      }
    }
    /*allocate the last row added*/
    m->data[rowsize-1] = xmalloc(sizeof(double)*colsize);

    /*copy the row value to the new row matrix*/
    for(i = 0; i < row->size; i++){
      m->data[rowsize-1][i] = row->data[i];
    }
  }
  else /*if(colsize <= m->col)*/{
    /*allocate the last row added*/
    m->data[rowsize-1] = xmalloc(sizeof(double)*colsize);
    for(i = 0; i < m->col; i++){
      if(i < row->size){
        m->data[rowsize -1][i] = row->data[i];
      }
      else{
        m->data[rowsize -1][i] = +0.f;
      }
    }
  }

  if(row->size > m->col)
    m->col = row->size;

  m->row += 1;

}

void MatrixAppendCol(matrix* m, dvector *col)
{
  size_t i, j;
  size_t lastcol;
  size_t colsize;
  size_t rowsize;

  if(m->col != 0){
    colsize =  m->col + 1;
  }
  else{
    /* redefine anyway the row size because the column have 0 size */
    colsize = 1;
  }

  if(m->row != 0){
    if(m->row < col->size){
      rowsize = col->size;
    }
    else{
      rowsize = m->row;
    }
  }
  else{
    rowsize = col->size;
  }

  if(m->row < rowsize){
    m->data = xrealloc(m->data, sizeof(double*)*rowsize);
  }

  for(i = 0; i < rowsize; i++){
    if(i < m->row)
      m->data[i] = xrealloc(m->data[i], sizeof(double)*colsize);
    else
      m->data[i] = xmalloc(sizeof(double)*colsize);
  }

  lastcol = m->col;

  if(rowsize < m->row){
    for(i = 0; i < m->row; i++ ){
      if(i < rowsize)
        m->data[i][lastcol] = col->data[i];
    else
      m->data[i][lastcol] = +0.f;
    }
  }
  else{
    if(rowsize > m->row){
      for(i = 0; i < rowsize; i++ ){
        m->data[i][lastcol] = col->data[i];
      }

      /*Fill with 0 value the new rows except the last column */
      for(i = m->row; i < rowsize; i++)
        for(j = 0; j < colsize-1; j++)
          m->data[i][j] = +0.f;
    }
    else{
      for(i = 0; i < rowsize; i++){
        m->data[i][lastcol] = col->data[i];
      }
    }
  }

  m->col = colsize; /* updating the column matrix size */
  m->row = rowsize; /* updating the row matrix size */
}


void MatrixAppendUIRow(matrix* m, uivector *row)
{
  size_t i, j;
  size_t rowsize =  m->row + 1;
  size_t colsize;


  if(m->col != 0){
    if(row->size > m->col)
      colsize = row->size;
    else /*if (row->size <= m->col)*/
      colsize = m->col;
  }
  else{
    colsize = row->size;
  }

  /*adding a new row*/
  m->data = xrealloc(m->data, sizeof(double*)*rowsize);

  if(colsize > m->col){
    /*resize the column*/
    for(i = 0; i < m->row; i++){
      m->data[i] = xrealloc(m->data[i], sizeof(double*)*colsize);
      /*initialize the new column value*/
      for(j = m->col; j < colsize; j++){
        m->data[i][j] = +0.f;
      }
    }
    /*allocate the last row added*/
    m->data[rowsize-1] = xmalloc(sizeof(double)*colsize);

    /*copy the row value to the new row matrix*/
    for(i = 0; i < row->size; i++){
      m->data[rowsize-1][i] = row->data[i];
    }
  }
  else /*if(colsize <= m->col)*/{
    /*allocate the last row added*/
    m->data[rowsize-1] = xmalloc(sizeof(double)*colsize);
    for(i = 0; i < m->col; i++){
      if(i < row->size){
        m->data[rowsize -1][i] = row->data[i];
      }
      else{
        m->data[rowsize -1][i] = +0.f;
      }
    }
  }

  if(row->size > m->col)
    m->col = row->size;

  m->row += 1;

}

void MatrixAppendUICol(matrix* m, uivector *col)
{
  size_t i, j;
  size_t lastcol;
  size_t colsize;
  size_t rowsize;

  if(m->col != 0){
    colsize =  m->col + 1;
  }
  else{
    /* redefine anyway the row size because the column have 0 size */
    colsize = 1;
  }

  if(m->row != 0){
    if(m->row < col->size){
      rowsize = col->size;
    }
    else{
      rowsize = m->row;
    }
  }
  else{
    rowsize = col->size;
  }

  if(m->row < rowsize){
    m->data = xrealloc(m->data, sizeof(double*)*rowsize);
  }

  for(i = 0; i < rowsize; i++){
    if(i < m->row)
      m->data[i] = xrealloc(m->data[i], sizeof(double)*colsize);
    else
      m->data[i] = xmalloc(sizeof(double)*colsize);
  }

  lastcol = m->col;

  if(rowsize < m->row){
    for(i = 0; i < m->row; i++ ){
      if(i < rowsize)
        m->data[i][lastcol] = col->data[i];
    else
      m->data[i][lastcol] = +0.f;
    }
  }
  else{
    if(rowsize > m->row){
      for(i = 0; i < rowsize; i++ ){
        m->data[i][lastcol] = col->data[i];
      }

      /*Fill with 0 value the new rows except the last column */
      for(i = m->row; i < rowsize; i++)
        for(j = 0; j < colsize-1; j++)
          m->data[i][j] = +0.f;
    }
    else{
      for(i = 0; i < rowsize; i++){
        m->data[i][lastcol] = col->data[i];
      }
    }
  }

  m->col = colsize; /* updating the column matrix size */
  m->row = rowsize; /* updating the row matrix size */
}


/* Description:
 * Delete a specific row in a matrix
 */
void MatrixDeleteRowAt(matrix *m, size_t row)
{
  size_t i, j, k;
  matrix *c;
  NewMatrix(&c, m->row, m->col);
  MatrixCopy(m, &c);
  ResizeMatrix(m, c->row-1, c->col);
  k = 0;
  for(i = 0; i < c->row; i++){
    if(i == row){
      continue;
    }
    else{
      for(j = 0; j < c->col; j++){
        m->data[k][j] = c->data[i][j];
      }
      k++;
    }
  }
  DelMatrix(&c);
}

/* Description:
 * Delete a specific column in a matrix
 */
void MatrixDeleteColAt(matrix *m, size_t col)
{
  size_t i, j, k;
  matrix *c;
  NewMatrix(&c, m->row, m->col);
  MatrixCopy(m, &c);
  ResizeMatrix(m, c->row, c->col-1);
  k = 0;
  for(j = 0; j < c->col; j++){
    if(j == col){
      continue;
    }
    else{
      for(i = 0; i < c->row; i++){
        m->data[i][k] = c->data[i][j];
      }
      k++;
    }
  }
  DelMatrix(&c);
}

/*
 * p[i] =   Σ m[i][j] * v[j]
 */
void MatrixDVectorDotProduct(matrix *m, dvector *v, dvector *p)
{
  /* (m*vect) where t is a column vector not transposed
     the size of the "vect" vector must be equal to the number of matrix row*/
  size_t i, j;
  double res;
  if(m->col == v->size){
    for(i = 0; i < m->row; i++){
      for(j = 0; j < m->col; j++){
        if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1) ||
           FLOAT_EQ(v->data[j], MISSING, 1e-1)){
          continue;
        }
        else{
          res = m->data[i][j] * v->data[j];
          if(_isnan_(res) || _isinf_(res)){
            continue;
          }
          else{
            p->data[i] += res;
          }
        }
      }
    }
  }
  else{
    fprintf(stdout,"MatrixDVectorDotProduct Error while calculating product (X*v)!!\n The column vector size must be equal to the matrix column size.\n");
    fflush(stdout);
    abort();
  }
}

/* MultiThreadMatrixDVectorDotProduct for large data! */
typedef struct{
  size_t from, to; /*used for slicing*/
  matrix *m; /* shared between threads */
  dvector *v; /* shared between threads */
  dvector *res; /* value modified in thread */
} tharg;

void *MatrixDVectorDotProductWorker(void *arg_){
  tharg *arg;
  double res;
  arg = (tharg*) arg_;
  size_t i, j;

  for(i = arg->from; i < arg->to; i++){
    arg->res->data[i] = 0.f;
    for(j = 0; j < arg->m->col; j++){
      if(FLOAT_EQ(arg->m->data[i][j], MISSING, 1e-1) ||
         FLOAT_EQ(arg->v->data[j], MISSING, 1e-1)){
        continue;
      }
      else{
        res = arg->m->data[i][j] * arg->v->data[j];
        if(_isnan_(res) || _isinf_(res)){
          continue;
        }
        else{
          arg->res->data[i] += res;
        }
      }
    }
  }
  return NULL;
}

void MT_MatrixDVectorDotProduct(matrix *m, dvector *v, dvector *p)
{
  if(m->col == v->size){
    /* (m*vect) where t is a column vector not transposed
      the size of the "vect" vector must be equal to the number of matrix row*/
    size_t th, nthreads;
    GetNProcessor(&nthreads, NULL);
    if(nthreads == 1){
      /* Redirect to the single thread version */
      MatrixDVectorDotProduct(m, v, p);
    }
    else{
      pthread_t *threads = xmalloc(sizeof(pthread_t)*nthreads);
      tharg *arg = xmalloc(sizeof(tharg)*nthreads);
      /* initialize threads arguments.. */
      size_t step = (size_t)ceil((double)m->row/(double)nthreads);
      size_t from = 0, to = step;
      for(th = 0; th < nthreads; th++){
        arg[th].m = m; /*SHARE THIS MATRIX. N.B.: THIS MATRIX MUST BE NOT MODIFIED DURING THE CALCULATION */
        arg[th].v = v; /*SHARE THIS VECTOR. N.B.: THIS VECTOR MUST BE NOT MODIFIED DURING THE CALCULATION */
        arg[th].res = p;
        arg[th].from = from;
        arg[th].to = to;
        pthread_create(&threads[th], NULL, MatrixDVectorDotProductWorker, (void*) &arg[th]);

        from = to;
        if(from+step > m->row){
          to = m->row;
        }
        else{
          to+=step;
        }
      }

      for(th = 0; th < nthreads; th++){
        pthread_join(threads[th], NULL);
      }

      xfree(threads);
      xfree(arg);
    }
  }
  else{
    fprintf(stdout,"MatrixDVectorDotProduct Error while calculating product (X*v)!!\n The column vector size must be equal to the matrix column size.\n");
    fflush(stdout);
    abort();
  }
}

/*
 * p[j] =   Σ v[i] * m[i][j]
 */
void DVectorMatrixDotProduct(matrix *m, dvector *v, dvector *p)
{
  size_t i, j;
  double res;
  /* (vect'*m) where vect is colunm vector transposed
     the size of the "vect" vector must be equal to the number of matrix column */
  if(m->row == v->size){
    for(j = 0; j < m->col; j++){
      for(i = 0; i < m->row; i++){
        if(FLOAT_EQ(v->data[i], MISSING, 1e-1) ||
           FLOAT_EQ(m->data[i][j], MISSING, 1e-1)){
          continue;
        }
        else{
          res = v->data[i] * m->data[i][j];
          if(_isnan_(res) || _isinf_(res)){
            continue;
          }
          else{
            p->data[j] += res;
          }
        }
      }
    }
  }
  else{
    fprintf(stdout,"DVectorMatrixDotProduct Error while calculating product of a (v'*X)!!\n The transposed column vector size must be equal to the matrix row size.\n");
    fflush(stdout);
    abort();
  }
}


void *DVectorMatrixDotProductWorker(void *arg_)
{
  tharg *arg;
  double res;
  arg = (tharg*) arg_;
  size_t i, j;

  for(j = arg->from; j < arg->to; j++){
    for(i = 0; i < arg->m->row; i++){
      if(FLOAT_EQ(arg->v->data[i], MISSING, 1e-1) ||
         FLOAT_EQ(arg->m->data[i][j], MISSING, 1e-1)){
           continue;
      }
      else{
        res = arg->v->data[i] * arg->m->data[i][j];
        if(_isnan_(res) || _isinf_(res)){
          continue;
        }
        else{
          arg->res->data[j] += res;
        }
      }
    }
  }
  return NULL;
}

 /*
 * p[j] =   Σ v[i] * m[i][j]
 */
void MT_DVectorMatrixDotProduct(matrix *m, dvector *v, dvector *p)
{
  /* (vect'*m) where vect is colunm vector transposed
     the size of the "vect" vector must be equal to the number of matrix column */
  if(m->row == v->size){
    size_t th, nthreads;
    GetNProcessor(&nthreads, NULL);
    if(nthreads == 1){
      DVectorMatrixDotProduct(m, v, p);
    }
    else{
      pthread_t *threads = xmalloc(sizeof(pthread_t)*nthreads);
      tharg *arg = xmalloc(sizeof(tharg)*nthreads);

      /* initialize threads arguments.. */
      size_t step = (size_t)ceil((double)m->col/(double)nthreads);
      size_t from = 0, to = step;
      for(th = 0; th < nthreads; th++){
        arg[th].m = m; /*SHARE THIS MATRIX. N.B.: THIS MATRIX MUST BE NOT MODIFIED DURING THE CALCULATION */
        arg[th].v = v; /*SHARE THIS VECTOR. N.B.: THIS VECTOR MUST BE NOT MODIFIED DURING THE CALCULATION */
        arg[th].res = p;
        arg[th].from = from;
        arg[th].to = to;
        pthread_create(&threads[th], NULL, DVectorMatrixDotProductWorker, (void*) &arg[th]);

        from = to;
        if(from+step > m->col){
          to = m->col;
        }
        else{
          to+=step;
        }
      }

      for(th = 0; th < nthreads; th++){
        pthread_join(threads[th], NULL);
      }

      xfree(threads);
      xfree(arg);
    }
  }
  else{
    fprintf(stdout,"DVectorMatrixDotProduct Error while calculating product of a (v'*X)!!\n The transposed column vector size must be equal to the matrix row size.\n");
    fflush(stdout);
    abort();
  }
}

void DVectorTrasposedDVectorDotProduct(dvector *v1, dvector *v2, matrix *m)
{
  size_t i, j;

  if(m->row != v1->size && m->col != v2->size)
    ResizeMatrix(m, v1->size, v2->size);

  for(i = 0; i < v1->size; i++){
    for(j = 0; j < v2->size; j++){
      if(FLOAT_EQ(v1->data[i], MISSING, 1e-1) ||
         FLOAT_EQ(v2->data[i], MISSING, 1e-1)){
        m->data[i][j] = MISSING;
      }
      else{
        m->data[i][j] = v1->data[i]*v2->data[j];
      }
    }
  }
}

/*
v/m = (inverse (m') * v')'
*/
void DVectorTransposedMatrixDivision(dvector *v, matrix *m, dvector *r)
{
  if(m->col == v->size){
    matrix *m_t, *m_inv;
    NewMatrix(&m_t, m->col, m->row);
    MatrixTranspose(m, m_t);
    initMatrix(&m_inv);
    MatrixInversion(m_t, m_inv);
    DVectorResize(r, m_inv->col);
    MatrixDVectorDotProduct(m_inv, v, r);
    DelMatrix(&m_inv);
    DelMatrix(&m_t);
  }
  else{
    fprintf(stderr, "Unable to compute vector / matrix. vector size %d != matrix column size %d\n", (int)v->size, (int)m->row);
    abort();
  }
}
/*
 * R = M'M
 * opp
 *
 * The product of an m x n matrix A and an n x p matrix B is an m x p matrix C where
 *
 * c[i][j] = Sum a[i][k]*b[k][j]
 */

/*
 * m_t is the transposed matrix of m. The result is a square matrix named r.
 * void MatrixDotProduct(matrix *m_t, matrix *m, matrix **r)
{
  int i, j, k;
//  (*r).row =(*r).col = m.col; square matrix
  for(i = 0; i < (*r)->row; i++){
    for(j = 0; j < m_t->row; j++){
      for(k = 0; k < m_t->col; k++){
        (*r)->data[j][i] += m_t->data[j][k] * m->data[k][i];
      }
    }
  }
}
*/

void MatrixDotProduct(matrix *a, matrix *b, matrix *r)
{
  if(a->col == b->row){
    size_t i, j, k;
    double res;
    for(i = 0; i < a->row; i++){ /* m */
      for(j = 0; j < b->col; j++){ /* p */
        for(k = 0; k < a->col; k++){ /* n */
          if(FLOAT_EQ(a->data[i][k], MISSING, 1e-1) ||
             FLOAT_EQ(b->data[k][j], MISSING, 1e-1)){
            continue;
          }
          else{
            res = a->data[i][k] * b->data[k][j];
            if(_isnan_(res) || _isinf_(res)){
              r->data[i][j] +=  +0.f;
            }
            else{
              r->data[i][j] +=  res;
            }
          }
        }
      }
    }
  }
  else{
    fprintf(stdout,"MatrixDotProduct Error!!\n The product of an m x n matrix A and an n x p matrix B is an m x p matrix C. col(a): %u != row(b) %u\n", (unsigned int)a->col, (unsigned int)b->row);
    fflush(stdout);
    abort();
  }
}


/*
 * Array product between two vector
 * X = t ⊗ p'
 *
 * X_ij = t_i*p_j
 *
 */
void RowColOuterProduct(dvector *a, dvector *b, matrix *m)
{
  size_t i, j;
  double res;
  for(i = 0; i < a->size; i++){
    for(j = 0; j < b->size; j++){
      if(FLOAT_EQ(a->data[i], MISSING, 1e-1) ||
         FLOAT_EQ(b->data[j], MISSING, 1e-1)){
        m->data[i][j] = MISSING;
      }
      else{
        res = a->data[i] * b->data[j];
        if(_isnan_(res) || _isinf_(res)){
          m->data[i][j] = +0.f;
        }
        else{
          m->data[i][j] = res;
        }
      }
    }
  }
}

void MatrixTranspose(matrix *m, matrix *r){
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      r->data[j][i] = m->data[i][j];
    }
  }
}

/* LU decomoposition of a general matrix */
extern void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
/* Generate inverse of a matrix given its LU decomposition */
extern void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

/*
 * Matrix Inversion according to the LU decomposition
 */
void MatrixLUInversion(matrix *m, matrix *m_inv)
{
  int N = m->row;
  double *M = xmalloc(sizeof(double)*m->row*m->col);
  int *IPIV = xmalloc(sizeof(int)*N);
  int LWORK = N*N;
  double *WORK = xmalloc(sizeof(double)*LWORK);
  int INFO;
  int i, j, k;
  k = 0;
  for(i = 0; i < m->col; i++){
    for(j = 0; j < m->row; j++){
      M[k] = m->data[i][j];
      k++;
    }
  }
  dgetrf_(&N,&N,M,&N,IPIV,&INFO);
  dgetri_(&N,M,&N,IPIV,WORK,&LWORK,&INFO);

  xfree(IPIV);
  xfree(WORK);

  ResizeMatrix(m_inv, m->row, m->col);
  k = 0;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      m_inv->data[i][j] = M[k];
      k++;
    }
  }
  xfree(M);
}

/*Gauss-Jordan
 *
 * 1) If the first row have the first element 0, then excange this with an other row that have the first element != 0. If all the row have the first element 0 go to the step 3.
 * 2) For each row (*AT)i with the first element != 0, except the first row considered, multiply the first row for a coefficient c that must
 *
 * Spiegazione in italiano:
 * Per ogni riga Ai con primo elemento non nullo, eccetto la prima (i > 1), moltiplica la prima riga per un coefficiente scelto in maniera tale che
 * la somma tra la prima riga e Ai abbia il primo elemento nullo (quindi coefficiente = − Ai1 / A11). Sostituisci Ai con la somma appena ricavata.
 *
 *
 */

void MatrixInversion(matrix *m, matrix *m_inv)
{
  if(m->row == m->col){
    size_t i, j, k;
    double ratio, a;

    matrix *AI;
    NewMatrix(&AI, m->row, m->col*2);

    /*copy the m1 value to AI*/
    for(i = 0; i < m->row; i++){
      for(j = 0; j < m->col; j++){
        AI->data[i][j] = m->data[i][j];
      }
    }

    /*build the identity matrix*/
    for(i = 0; i < m->row; i++){
      for(j = m->col; j < 2*m->col; j++){
        if(i==(j-m->row)){
          AI->data[i][j] = 1.f;
        }
        else{
          AI->data[i][j] = +0.f;
        }
      }
    }


    for(i = 0; i < m->row; i++){
      for(j = 0; j < m->col; j++){
        if(i!=j){
          ratio = AI->data[j][i] / AI->data[i][i];
          for(k = 0; k < 2*m->col; k++){
            AI->data[j][k] -= AI->data[i][k]*ratio;
          }
        }
      }
    }

    for(i = 0; i < m->row; i++){
      a = AI->data[i][i];
      for(j = 0; j < 2*m->col; j++){
        AI->data[i][j] = AI->data[i][j]/a;
      }
    }

    if(m_inv->row != m->row || m_inv->col != m->col){
      ResizeMatrix(m_inv, m->row, m->col);
    }

    for(i = 0; i < m->row; i++){
      for(j = 0; j < m->col; j++){
        m_inv->data[i][j] = AI->data[i][m->col+j];
      }
    }

    DelMatrix(&AI);
  }
  else{
    fprintf(stdout,"Matrix Inversion Error!\n The matrix to invert must be squared!\n");
    fflush(stdout);
    abort();
  }
}

/*
 * Let A=UΣVT be an SVD of A, i.e., U and V are orthogonal and
 * Σ=diag(σ1,σ2,…,σr,0,…,0) is real diagonal with σk>0 for all k=1,…,r.
 *
 * Then
 *
 * A+=VΣ+UT,
 *
 * with
 *
 * Σ+=diag(σ_1^−1,σ_2^-1,…,σ_r^-1,0,…,0).
 *
 */
void MatrixPseudoinversion(matrix *m, matrix *m_inv)
{
  matrix *U, *S, *V_T;
  initMatrix(&U);
  initMatrix(&S);
  initMatrix(&V_T);
  SVD(m, U, S, V_T);

  /*
  puts("U");
  PrintMatrix(U);
  puts("S");
  PrintMatrix(S);
  puts("VT");
  PrintMatrix(V_T);
  */

  matrix *Sinv;
  NewMatrix(&Sinv, S->col, S->row);
  MatrixInversion(S, Sinv);
  DelMatrix(&S);

  matrix *USinv;
  NewMatrix(&USinv, U->row, Sinv->col);
  MatrixDotProduct(U, Sinv, USinv);
  DelMatrix(&U);
  DelMatrix(&Sinv);
  ResizeMatrix(m_inv, m->row, m->col);
  MatrixDotProduct(USinv, V_T, m_inv);
  DelMatrix(&V_T);
  DelMatrix(&USinv);
}


void MatrixMoorePenrosePseudoinverse(matrix *m, matrix *inv)
{
  /*A+ = (A'A)-1 A'*/
  matrix *m_T, *m_t_m, *i_m_t_m;

  /*A'*/
  NewMatrix(&m_T, m->col, m->row);
  MatrixTranspose(m, m_T);

  /*A'A*/
  NewMatrix(&m_t_m, m_T->row, m->col);
  MatrixDotProduct(m_T, m, m_t_m);

  /*(A'A)-1*/
  NewMatrix(&i_m_t_m, m_t_m->row, m_t_m->col);
  //MatrixLUInversion(m_t_m, i_m_t_m);
  MatrixPseudoinversion(m_t_m, i_m_t_m);
  DelMatrix(&m_t_m);

  /*(A'A)-1 A*/
  ResizeMatrix(inv, i_m_t_m->row, m_T->col);
  MatrixDotProduct(i_m_t_m, m_T, inv);

  DelMatrix(&i_m_t_m);
  DelMatrix(&m_T);
}

void GenIdentityMatrix(matrix *m)
{
  size_t i;
  if(m->row == m->col){
    for(i=0; i < m->row; i++){
      m->data[i][i]=1;
    }
  }
}


double MatrixTrace(matrix *m)
{
  size_t i;
  double tr = 0.f;
  if(m->row == m->col){
    for(i=0; i < m->row; i++){
      tr += m->data[i][i];
    }
  }
  return tr;
}

void MeanCenteredMatrix(matrix *m, matrix *mc)
{
  size_t i, j, n;
  double average, res;
  for(j = 0; j < m->col; j++){
    if(m->row > 1){
      /*Calculate the average */
      average = +0.f;
      n = 0;
      for(i = 0; i < m->row; i++ ){
        if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1)){
          continue;
        }
        else{
          average += m->data[i][j];
          n++;
        }
      }

      average /= (double)n;

      for(i = 0; i < m->row; i++){
        if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1)){
          continue;
        }
        else{
          res = m->data[i][j] - average;
          if(_isnan_(res) || _isinf_(res)){
            (*mc).data[i][j] = +0.f;
          }
          else{
            (*mc).data[i][j] = res;
          }
        }
      }
    }
    else{
      for(i = 0; i < m->row; i++){
        (*mc).data[i][j] = m->data[i][j];
      }
    }
  }
}

/*
 * The equation for the Pearson product moment correlation coefficient, r, is:
 * RSQ = Sum (x_i -x_med)*(y_i - y_med) / sqrt( Sum (x_i - x_med)^2 Sum (y_i - y_med)^2 )
 */

void PearsonCorrelMatrix(matrix* msrc, matrix* mdst)
{
  size_t i, j, k;
  double n, a, b, xres, yres;
  dvector *mean;
  ResizeMatrix(mdst, msrc->col, msrc->col);
  initDVector(&mean);
  MatrixColAverage(msrc, mean);
  for(k = 0; k < msrc->col; k++){
    mdst->data[k][k] = 1.f;
    for(j = k+1; j < msrc->col; j++){
      n = a = b = +0.f;
      for(i = 0; i < msrc->row; i++){
        xres = msrc->data[i][k] - mean->data[k];
        yres = msrc->data[i][j] - mean->data[j];
        n += xres * yres;
        a += square(xres);
        b += square(yres);
      }
      if((int)floor(a*b) == 0){
        mdst->data[k][j] = mdst->data[j][k] = +0.f;
      }
      else{
        mdst->data[k][j] = mdst->data[j][k] = square(n / sqrt(a*b));
      }
    }
  }
  DelDVector(&mean);
}

/*
 * The equation for the Spearman product moment correlation coefficient, r, is:
 * rho = 1 - (Sum (rank_x_i - rank_y_i)^2 / (n(n^2-1)))
 * where n are the observation
 */


void SpearmanCorrelMatrix(matrix* msrc, matrix* mdst)
{
  size_t i, j, k, l;
  matrix *rankm;
  dvector *vtosort;
  double n;
  ResizeMatrix(mdst, msrc->col, msrc->col);
  NewMatrix(&rankm, msrc->row, 4);
  NewDVector(&vtosort, msrc->row);
  for(k = 0; k < msrc->col; k++){
    mdst->data[k][k] = 1.f;
    for(i = 0; i < msrc->row; i++){
      rankm->data[i][0] = vtosort->data[i] = msrc->data[i][k];
    }

    DVectorSort(vtosort);

    for(i = 0; i < vtosort->size; i++){
      for(j = 0; j < rankm->row; j++){
        if(FLOAT_EQ(vtosort->data[i], rankm->data[j][0], EPSILON)){
          rankm->data[j][2] = i+1;
          break;
        }
        else{
          continue;
        }
      }
    }

    for(j = k+1; j < msrc->col; j++){
      for(i = 0; i < msrc->row; i++){
        rankm->data[i][1] = vtosort->data[i] = msrc->data[i][j];
      }

      /* rank second column */
      DVectorSort(vtosort);
      for(i = 0; i < vtosort->size; i++){
        for(l = 0; l < rankm->row; l++){
          if(FLOAT_EQ(vtosort->data[i], rankm->data[l][1], EPSILON)){
            rankm->data[l][3] = i+1;
            break;
          }
          else{
            continue;
          }
        }
      }

      /*calculate d^2*/
      n = +0.f;
      for(i = 0; i < rankm->row; i++){
        n += square(rankm->data[i][2] - rankm->data[i][3]);
      }
      mdst->data[k][j] = mdst->data[j][k] = 1 - ((6*n) / (rankm->row*((square(rankm->row)-1))));
    }
  }
  DelMatrix(&rankm);
  DelDVector(&vtosort);
}

void MatrixColAverage(matrix *m, dvector *colaverage)
{
  size_t i, j, n;

  double average;
  for(j = 0; j < m->col; j++){
    /*Calculate the average */
    average = +0.f;
    n = 0;
    for(i = 0; i < m->row; i++){
      if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1))
        continue;
      else{
        average += m->data[i][j];
        n++;
      }
    }

    if(FLOAT_EQ(average, 0.f, 1e-6))
        average = 0.f;
    else
        average /= (double)n;
    DVectorAppend(colaverage, average);
  }
}

void MatrixRowAverage(matrix *m, dvector *rowaverage)
{
  size_t i, j, n;

  double average;
  for(i = 0; i < m->row; i++){
    /*Calculate the average */
    average = +0.f;
    n = 0;
    for(j = 0; j < m->col; j++){
      if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1))
        continue;
      else{
        average += m->data[i][j];
        n++;
      }
    }
    average /= (double)n;
    DVectorAppend(rowaverage, average);
  }
}

void MatrixColSDEV(matrix* m, dvector *colsdev)
{
  size_t i, j, n;
  double var, average;
  for(j = 0; j < m->col; j++){
    average=+0.f;
    n = 0;
    for(i = 0; i < m->row; i++){
      if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1))
        continue;
      else{
        average += m->data[i][j];
        n++;
      }
    }

    /* average of the column j;*/
    average /= (double)n;

    var = +0.f;
    n = 0;
    for(i = 0; i< m->row; i++){
      if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1))
        continue;
      else{
        var += square(m->data[i][j] - average);
        n++;
      }
    }
    /* sample variance: is used whe the average of data is not known so you need to extimate the data average */
    var = var/(n-1);

    /* standard deviation calculation */
    DVectorAppend(colsdev, sqrt(var));
  }
}

void MatrixColRMS(matrix* m, dvector *colrms)
{
  size_t i, j, n;
  double a;

  for(j = 0; j < m->col; j++){
    a = +0.f;
    n = 0;
    for(i = 0; i < m->row; i++){
      if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1))
        continue;
      else{
        a += square(m->data[i][j]);
        n++;
      }
    }

    /* average of the column j;*/
    a /= (double)n;
    /* standard deviation calculation */
    DVectorAppend(colrms, sqrt(a));
  }
}

void MatrixColVar(matrix* m, dvector *colvar)
{
  size_t i, j, n;
  double var, average;
  for(j = 0; j < m->col; j++){
    average = +0.f;
    n = 0;
    for(i = 0; i < m->row; i++){
      if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1))
        continue;
      else{
        average += m->data[i][j];
        n++;
      }
    }

    /* average of the column j;*/
    average /= (double)n;

    var = +0.f;
    n = 0;
    for(i = 0; i<m->row; i++){
      if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1))
        continue;
      else{
        double a = (m->data[i][j] - average);
        var += a*a;
        n++;
      }
    }
    /* sample variance: is used whe the average of data is not known so you need to extimate the data average */
    var = var/(double)(n-1);

    /* standard deviation calculation */
    DVectorAppend(colvar, var);
  }
}

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
 *  - Coefficient of variation Population (CV)
 *  - Coefficient of variation Sample (CV)
 *  - N. zeros
 *  - N. missing values
 */
void MatrixColDescStat(matrix *m, matrix *ds)
{
  int i, j, n;
  size_t n_zeros;
  size_t n_missing;
  double avg;
  double median;
  double armonic;
  double var;
  double min;
  double max;
  dvector *v;

  /* n = m->row if no MISSING value */
  ResizeMatrix(ds, m->col, 13);
  NewDVector(&v, m->row);

  for(j = 0; j < m->col; j++){
    avg = 0.f;
    var = 0.f;
    median = 0.f;
    armonic = 0.f;

    /* get column min and max,
     * cacluate the average, the median and the armonic average
     */
    min = max = m->data[0][j];
    n_zeros = 0;
    n_missing = 0;
    n = 0;
    for(i = 0; i < m->row; i++){
      if(FLOAT_EQ(m->data[i][j], 0.f, 1e-6))
        n_zeros++;

      if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1)){
        n_missing++;
      }
      else{
        avg += m->data[i][j];
        armonic += 1.f/m->data[i][j];
        v->data[i] = m->data[i][j];
        if(m->data[i][j] > max)
          max = m->data[i][j];

        if(m->data[i][j] < min)
          min = m->data[i][j];
        n++;
      }
    }
    avg /= (double)n;

    DVectorMedian(v, &median);
    for(i = 0; i < m->row; i++){
      if(FLOAT_EQ(m->data[i][j], MISSING, 1e-1)){
        continue;
      }
      else{
        var += square(m->data[i][j] - avg);
      }
    }

    ds->data[j][0] = avg;
    ds->data[j][1] = median;
    ds->data[j][2] = (double)(m->row)/armonic;
    ds->data[j][3] = var/(double)n;
    ds->data[j][4] = var/(double)(n-1);
    ds->data[j][5] = sqrt(var/(double)n);
    ds->data[j][6] = sqrt(var/(double)(n-1));
    ds->data[j][7] = ds->data[j][5]/avg * 100;
    ds->data[j][8] = ds->data[j][6]/avg * 100;
    ds->data[j][9] = min;
    ds->data[j][10] = max;
    ds->data[j][11] = n_zeros;
    ds->data[j][12] = n_missing;
  }
  DelDVector(&v);
}

/* calculation of the covariance matrix */
void MatrixCovariance(matrix *m, matrix *cm)
{
  size_t i, j, k;
  double sum;
  dvector *colaverage;

  ResizeMatrix(cm, m->col, m->col);

  initDVector(&colaverage);
  MatrixColAverage(m, colaverage);


  for(i = 0; i < m->col; i++){
    for(j = 0; j < m->col; j++){
      sum = +0.f;
      for(k =0; k < m->row; k++){
        sum += (m->data[k][i] - colaverage->data[i]) * (m->data[k][j] - colaverage->data[j]);
      }
      cm->data[i][j] = sum/(m->row-1);
    }
  }

  DelDVector(&colaverage);
}

/* Transform a matrix into a logaritmic matrix */
void Matrix2LogMatrix(matrix *m_in, matrix *m_out)
{
  size_t i, j;
  ResizeMatrix(m_out, m_in->row, m_in->col);
  for(i = 0; i < m_in->row; i++){
    for(j = 0; j < m_in->col; j++){
      m_out->data[i][j] = log10(m_in->data[i][j]+1);
    }
  }
}

/* Transform a matrix into a SQUARE matrix */
void Matrix2SquareMatrix(matrix *m_in, matrix *m_out)
{
  size_t i, j;
  ResizeMatrix(m_out, m_in->row, m_in->col);
  for(i = 0; i < m_in->row; i++)
    for(j = 0; j < m_in->col; j++)
      m_out->data[i][j] = square(m_in->data[i][j]);
}

/* Transform a matrix into a SQRT matrix */
void Matrix2SQRTMatrix(matrix *m_in, matrix *m_out)
{
  size_t i, j;
  ResizeMatrix(m_out, m_in->row, m_in->col);
  for(i = 0; i < m_in->row; i++)
    for(j = 0; j < m_in->col; j++)
      m_out->data[i][j] = sqrt(m_in->data[i][j]);
}

/* Transform a matrix into ABS matrix */
void Matrix2ABSMatrix(matrix *m_in, matrix *m_out)
{
  size_t i, j;
  ResizeMatrix(m_out, m_in->row, m_in->col);
  for(i = 0; i < m_in->row; i++)
    for(j = 0; j < m_in->col; j++)
      m_out->data[i][j] = fabs(m_in->data[i][j]);
}


/* Develop an interaction factors matrix
 * Es. Use in DOE
 * WARNING: EXPERIMENTAL!!
 */
void Matrix2IntFactorsMatrix(matrix *m_in, size_t factors, matrix *m_out)
{
  size_t i, j, k, l, c;
  size_t nifc = 0;
  for(i = 1; i < m_in->col; i++)
    nifc += i;
  ResizeMatrix(m_out, m_in->row, (size_t)nifc+(2*m_in->col));
  for(i = 0; i < m_in->row; i++){
    puts("#######");
    for(k = 0, c = 0; k < factors; k++){
      for(j = 0; j < m_in->col; j++){
        double res = m_in->data[i][j];
        for(l = 0; l < k; l++){
          res *= m_in->data[i][j+l];
          printf("%d * %d \n", (int)j, (int)l);
        }
        m_out->data[i][c] = res;
        c++;
      }
    }
  }
}


/* Transform a matrix into a row centered scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixRowCenterScaling(matrix *m_in, matrix *m_out)
{
  size_t i, j;
  double rowsum;
  ResizeMatrix(m_out, m_in->row, m_in->col);
  for(i = 0; i < m_in->row; i++){
    rowsum = 0.f;
    for(j = 0; j < m_in->col; j++){
      rowsum += m_in->data[i][j];
    }

    for(j = 0; j < m_in->col; j++){
      m_out->data[i][j] = m_in->data[i][j]/rowsum;
    }
  }
}

/* Transform a matrix into a SVN row scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixSVNScaling(matrix *m_in, matrix *m_out)
{
  size_t i, j;
  double rowaverage, rowstdev;
  ResizeMatrix(m_out, m_in->row, m_in->col);
  for(i = 0; i < m_in->row; i++){
    rowaverage = 0.f;
    for(j = 0; j < m_in->col; j++){
      rowaverage += m_in->data[i][j];
    }
    rowaverage /= (double)m_in->col;

    rowstdev = 0.f;
    for(j = 0; j < m_in->col; j++){
      rowstdev += square(m_in->data[i][j]-rowaverage);
    }

    rowstdev /= (double)(m_in->col-1);
    rowstdev = sqrt(rowstdev);

    for(j = 0; j < m_in->col; j++){
      m_out->data[i][j] = (m_in->data[i][j]-rowaverage)/rowstdev;
    }
  }
}


/*
 * ||X|| = the square root of the sum of the squares of all the elements in the matrix
 */
double Matrixnorm(matrix *m)
{
  size_t i, j;
  double norm = +0.f, v;
  for(j = 0; j < m->col; j++){
    for(i = 0; i < m->row; i++){
      v = m->data[i][j];
      if(_isnan_(v) || !isfinite(v)){
        continue;
      }
      else{
        norm += square(v);
      }
    }
  }
  return sqrt(norm);
}

double Matrix1norm(matrix *m)
{
  size_t i, j;
  double norm = +0.f;
  for(j = 0; j < m->col; j++){
    for(i = 0; i < m->row; i++){
      norm += fabs(m->data[i][j]);
    }
  }
  return norm;
}

double MatrixDeterminant(matrix *m)
{
  if(m->row == m->col){
    size_t i, j, k, l;
    double d = 0;
    matrix *sub_m;

    if (m->row < 1){
      return MISSING;
    }
    else if(m->row == 1){
      d = m->data[0][0];
    }
    else if (m->row == 2){
      d = m->data[0][0] * m->data[1][1] - m->data[1][0] * m->data[0][1];
    }
    else{
      d = 0.;
      for(k = 0; k < m->row; k++){
        NewMatrix(&sub_m, m->row-1, m->row-1);
        for(i = 1; i < m->row; i++){
          l = 0;
          for(j = 0; j < m->row; j++){
            if(j == k)
              continue;
            else{
              sub_m->data[i-1][l] = m->data[i][j];
              l++;
            }
          }
        }
        d += pow(-1.0, 1.0+k+1.0) * m->data[0][k] * MatrixDeterminant(sub_m);
        DelMatrix(&sub_m);
      }
    }
    return d;
  }
  else{
    /* Error! */
    return MISSING;
  }
}

/*
double MatrixDeterminant(matrix *m)
{
  if(m->row == m->col){
    return Determinant(m->data, m->row);
  }
  else{
    return MISSING;
  }
}
*/

void MatrixNorm(matrix *m, matrix *nm)
{
  if((*nm).row == m->row && m->col == (*nm).col){
    size_t i, j;
    double mod, res;

    mod = Matrixnorm(m);

    for(i = 0; i < m->row; i++){
      for(j = 0; j < m->col; j++){
        res = m->data[i][j]/mod;
        if(_isnan_(res) || _isinf_(res)){
          nm->data[i][j] = +0.f;
        }
        else{
          nm->data[i][j] = res;
        }
      }
    }
  }
  else{
    fprintf(stdout,"MatrixNorm Error!\n");
    fflush(stdout);
    abort();
  }
}

void MatrixColumnMinMax(matrix* m, size_t col, double* min, double* max)
{
  if(m->row > 0 && col < m->col ){
    size_t i;
    double a;
    (*min) = (*max) = m->data[0][col];
    for(i = 1; i < m->row; i++){
      a = m->data[i][col];
      if(FLOAT_EQ(a, MISSING, 1e-1)){
        continue;
      }
      else{
        if(a < (*min)){
          (*min) = a;
        }

        if(a > (*max)){
          (*max) = a;
        }
      }
    }
  }
  else{
    fprintf(stdout,"Get Column Max Min Error!\n");
    fflush(stdout);
    (*min) = (*max) = MISSING;
  }
}

void MatrixSort(matrix* m, size_t col_n)
{
  size_t i, j, k;
  double temp;
  for(i = 0; i < m->row; i++){
    for(j = i+1; j < m->row; j++){
      if(m->data[i][col_n] > m->data[j][col_n]){
        for(k = 0; k < m->col; k++){
          temp = m->data[i][k];
          m->data[i][k] = m->data[j][k];
          m->data[j][k] = temp;
        }
      }
      else{
        continue;
      }
    }
  }
}

void MatrixReverseSort(matrix* m, size_t col_n)
{
  size_t i, j, k;
  double temp;
  for(i = 0; i < m->row; i++){
    for(j = i+1; j < m->row; j++){
      if(m->data[i][col_n] < m->data[j][col_n]){
        for(k = 0; k < m->col; k++){
          temp = m->data[i][k];
          m->data[i][k] = m->data[j][k];
          m->data[j][k] = temp;
        }
      }
      else{
        continue;
      }
    }
  }
}

void MatrixGetMaxValueIndex(matrix* m, size_t* row, size_t* col)
{
  size_t i, j;
  double tmp_value, best_value;

  if(col != NULL)
    (*col) = 0;

  if(row != NULL)
    (*row) = 0;

  best_value = m->data[0][0];

  for(j = 0; j < m->col; j++){
    for(i = 1; i < m->row; i++){
      tmp_value = m->data[i][j];
      if(tmp_value > best_value || FLOAT_EQ(tmp_value, best_value, EPSILON)){
        best_value = tmp_value;
        if(col != NULL)
          (*col) = j;

        if(row != NULL)
          (*row) = i;
      }
    }
  }
}


void MatrixGetMinValueIndex(matrix* m, size_t* row, size_t* col)
{
  size_t i, j;
  double tmp_value, best_value;

  if(col != NULL)
    (*col) = 0;

  if(row != NULL)
    (*row) = 0;

  best_value = m->data[0][0];

  for(j = 0; j < m->col; j++){
    for(i = 1; i < m->row; i++){
      tmp_value = m->data[i][j];
      if(tmp_value < best_value || FLOAT_EQ(tmp_value, best_value, EPSILON)){
        best_value = tmp_value;
        if(col != NULL)
          (*col) = j;

        if(row != NULL)
          (*row) = i;
      }
    }
  }
}


int cmpfunc(const void *a, const void *b )
{
  const double *a_ = *(const double **)a;
  const double *b_ = *(const double **)b;
  return b_[0] - a_[0];
}

void SVD(matrix* m, matrix *U, matrix *S, matrix *VT)
{
  size_t i;
  matrix *w1, *w2, *m_t, *v/*, *to_sort*/;
  dvector *eval1, *eval2;
  NewMatrix(&w1, m->row, m->row); // A A^T
  NewMatrix(&w2, m->col, m->col); // A^T A
  NewMatrix(&m_t, m->col, m->row);

  MatrixTranspose(m, m_t);

  MatrixDotProduct(m, m_t, w1);
  MatrixDotProduct(m_t, m, w2);

  initDVector(&eval1);
  initDVector(&eval2);

  initMatrix(&v);
  EVectEval(w1, eval1, v);
  EVectEval(w2, eval2, U);

  ResizeMatrix(VT, v->col, v->row);
  MatrixTranspose(v, VT);
  ResizeMatrix(S, m->row, m->col);

  /*NewMatrix(&to_sort, (*S)->row, 2);*/

  for(i = 0; i < S->col; i++){
    if(FLOAT_EQ(eval1->data[i], 0.f, 1e-6) || eval1->data[i] < 0)
      S->data[i][i] = 0.f;
    else{
      S->data[i][i] = sqrt(eval1->data[i]);
    }
  }

  DelMatrix(&v);
  DelMatrix(&m_t);
  DelDVector(&eval1);
  DelDVector(&eval2);
  DelMatrix(&w2);
  DelMatrix(&w1);
}

/* DGESDD prototype */
extern void dgesdd_( char* jobz, int* m, int* n, double* a,
                int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
                double* work, int* lwork, int* iwork, int* info );

void conv2matrix(int m, int n, double* a, int lda, matrix *m_)
{
  size_t i, j;
  ResizeMatrix(m_, m, n);
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++)
      m_->data[i][j] = a[i+j*lda];
  }
}

void SVDlapack(matrix *m_, matrix *u, matrix *s, matrix *vt)
{
  int i, j, k;
  int m = m_->row, n = m_->col, lda = m_->row, ldu = m_->row, ldvt = m_->col, info, lwork;

  double wkopt;
  double* work = NULL;
  /* Local arrays */
  double *s_, *u_, *vt_, *a;
  s_ = xmalloc(sizeof(double)*m_->row);
  u_ = xmalloc(sizeof(double)*m_->row*m_->row);
  vt_ = xmalloc(sizeof(double)*m_->col*m_->col);
  a = xmalloc(sizeof(double)*m_->row*m_->col);
  k = 0;
  for(j = 0; j < m_->col; j++){
    for(i = 0; i < m_->row; i++){
      a[k] = m_->data[i][j];
      k+=1;
    }
  }

  /* dgesdd_ implementation */
  int iwork[8*n];
  lwork = -1;
  /* U * SIGMA * V^H */
  dgesdd_( "Singular vectors", &m, &n, a, &lda, s_, u_, &ldu, vt_, &ldvt, &wkopt, &lwork, iwork, &info );
  lwork = (int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );
  /* Compute SVD */
  dgesdd_( "Singular vectors", &m, &n, a, &lda, s_, u_, &ldu, vt_, &ldvt, work, &lwork, iwork, &info );

  if(info > 0) {
    printf("SVD convergence failed!\n" );
    return;
  }

  /* s are the eigenvectors singular values diagonal matrix*/
  ResizeMatrix(s, m, n);
  for(i = 0; i < m_->col; i++){
    s->data[i][i] = s_[i];
  }
  //conv2matrix(1, n, s_, 1, s);
  /* u is left singular vectors */
  conv2matrix(m, n, u_, ldu, u);
  /*vt is the right singular vectors */
  conv2matrix( n, n, vt_, ldvt, vt);
  /* Free workspace */
  xfree(work);
  xfree(a);
  xfree(s_);
  xfree(vt_);
  xfree(u_);
}

void print_eigenvalues( char* desc, int n, double* wr, double* wi ) {
  int j;
  printf( "\n %s\n", desc );
  for( j = 0; j < n; j++ ) {
    if( wi[j] == (double)0.0 ) {
      printf( " %6.2f", wr[j] );
    } else {
      printf( " (%6.2f,%6.2f)", wr[j], wi[j] );
    }
  }
  printf( "\n" );
}

/* Auxiliary routine: printing eigenvectors */
void print_eigenvectors( char* desc, int n, double* wi, double* v, int ldv ) {
  int i, j;
  printf( "\n %s\n", desc );
  for( i = 0; i < n; i++ ) {
    j = 0;
    while( j < n ) {
      if( wi[j] == (double)0.0 ) {
        printf( " %6.2f", v[i+j*ldv] );
        j++;
      } else {
        printf( " (%6.2f,%6.2f)", v[i+j*ldv], v[i+(j+1)*ldv] );
        printf( " (%6.2f,%6.2f)", v[i+j*ldv], -v[i+(j+1)*ldv] );
        j += 2;
      }
    }
    printf( "\n" );
  }
}

extern void dgeev_(char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info);

void EVectEval(matrix *m, dvector *eval, matrix *evect)
{
  /* Locals */
  size_t i, j, k, N;

  N = m->row;
  int LDA = N;
  int LDVL = N;
  int LDVR = N;

  int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
  double wkopt;
  double* work;

  /* Local arrays */
  double *wr = xmalloc(sizeof(double)* N);
  double *wi = xmalloc(sizeof(double)* N);
  double *vl = xmalloc(sizeof(double)* N*LDVL);
  double *vr = xmalloc(sizeof(double)* N*LDVR);

  double *a = xmalloc(sizeof(double)*N*LDA);

  k = 0;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      a[k] = m->data[i][j];
      k++;
    }
  }

  /* Query and allocate the optimal workspace */
  lwork = -1;
  dgeev_("N", "V", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info);
  lwork = (int)wkopt;
  work = (double*)xmalloc( lwork*sizeof(double) );
  /* Solve eigenproblem */

  dgeev_("N", "V", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  /* Check for convergence */
  if(info > 0){
    printf( "The algorithm failed to compute eigenvalues.\n" );
  }

//  print_eigenvalues( "Eigenvalues", n, wr, wi );
//    Print left eigenvectors
//  print_eigenvectors( "Left eigenvectors", n, wi, vl, ldvl );
//    Print right eigenvectors
//  print_eigenvectors( "Right eigenvectors", n, wi, vr, ldvr );

  DVectorResize(eval, N);
  /*Eigenvalues
  * get only the real part
  */
  for(i = 0; i < N; i++){
    eval->data[i] = wr[i];
  }

  /*Eigenvectors Right
  *get only the real part!
  */
  ResizeMatrix(evect, N, N);
  for(i = 0; i < N; i++ ) {
    for(j = 0; j < N; j++){
       evect->data[i][j] = vr[i+j*ldvr];
    }
  }

  /* Free workspace */
  xfree(wr);
  xfree(wi);
  xfree(vl);
  xfree(vr);
  xfree((void*)work);
  xfree(a);
}

void QRMatrixVectNorm(matrix *m, size_t col, dvector *nv)
{
  size_t i;
  double s = +0.f;

  if(nv->size != m->row)
    DVectorResize(nv, m->row);

  for(i = 0; i < m->row; i++){
    s += square(m->data[i][col]);
  }
  s = sqrt(s);

  for(i = 0; i < m->row; i++){
    nv->data[i] = m->data[i][col]/s;
  }
}

void QRDecomposition(matrix *m, matrix *Q, matrix *R)
{
  size_t i, j, k;
  double D, p;
  dvector *d, *v;
  matrix *a, *a1, *P, *PQ;
  initMatrix(&a);
  MatrixCopy(m, &a);
  NewMatrix(&a1, a->row, a->col);

  NewDVector(&d, a->row);
  NewDVector(&v, a->row);


  NewMatrix(&P, v->size, v->size);
  NewMatrix(&PQ, a->row, a->row);

  for(k = 0; k < m->col && k < m->row - 1; k++){
    QRMatrixVectNorm(a, k, d);

    D = +0.f;
    for(i = k; i < d->size; i++){
      D += d->data[i]*d->data[i];
    }
    D = sqrt(D);

    if(d->data[k] > 0)
      D = -D;

    for(i = 0; i < k; i++)
      v->data[i] = +0.f;

    v->data[k] = sqrt((1/2.)* (1 - (d->data[k]/D)));
    p = -1*D*v->data[k];

    for(i = k+1; i < v->size; i++){
      v->data[i] = d->data[i]/(2*p);
    }

    for(i = 0; i < v->size; i++){
      for(j = 0; j < v->size; j++){
        P->data[i][j] = -2 * v->data[i] * v->data[j];
      }
    }

    for(i = 0; i < v->size; i++){
      P->data[i][i] += 1;
    }

    MatrixSet(a1, +0.f);
    MatrixDotProduct(P, a, a1);

    if(k == 0){
      MatrixCopy(P, &Q);
    }
    else{
      MatrixDotProduct(P, Q, PQ);
      MatrixCopy(PQ, &Q);
      MatrixSet(PQ, +0.f);
    }

    MatrixCopy(a1, &a);
    MatrixSet(P, +0.f);
  }

  for(i = 0; i < Q->row; i++){
    if(FLOAT_EQ(Q->data[i][i], 0, 1e-6))
      Q->data[i][i] = 1.f;
  }

  ResizeMatrix(R, m->row, m->col);
  MatrixDotProduct(Q, m, R);

  MatrixCopy(Q, &a1);

  MatrixTranspose(a1, Q);

  DelMatrix(&P);
  DelMatrix(&PQ);
  DelMatrix(&a1);
  DelMatrix(&a);
  DelDVector(&d);
  DelDVector(&v);
}

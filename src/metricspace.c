/* metricspace.c
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

#include "memwrapper.h"
#include "metricspace.h"
#include "numeric.h"
#include <math.h>
#include <pthread.h>

enum cmethod {EUCLIDEAN = 0,
              SQUARE_EUCLIDEAN,
              MANHATTAN,
              COSINE}; 

typedef struct{
  matrix *m1, *m2;
  matrix *distances;
  dvector *condensed_distances;
  size_t r_from, r_to;
  enum cmethod method;
} dst_th_arg;

void *CalcWorker(void *arg_)
{
  size_t i;
  size_t j;
  size_t k;
  dst_th_arg *arg = (dst_th_arg*) arg_;
  if(arg->method == EUCLIDEAN ||
     arg->method == SQUARE_EUCLIDEAN){
    double dist;
    for(i = arg->r_from; i < arg->r_to; i++){
      for(k = 0; k < arg->m2->row; k++){
        dist = 0.f;
        for(j = 0; j < arg->m2->col; j++){
          dist += square(arg->m1->data[i][j] - arg->m2->data[k][j]);
        }

        if(arg->method == SQUARE_EUCLIDEAN){
          arg->distances->data[k][i] = dist;
        }
        else{
          arg->distances->data[k][i] = sqrt(dist);
        }
      }
    }
  }
  else if(arg->method == MANHATTAN){
    double dist;
    for(i = arg->r_from; i < arg->r_to; i++){
      for(k = 0; k < arg->m2->row; k++){
        dist = 0.f;
        for(j = 0; j < arg->m2->col; j++){
          dist += fabs(arg->m1->data[i][j] - arg->m2->data[k][j]);
        }
        arg->distances->data[k][i] = dist;
      }
    }
  }
  else if(arg->method == COSINE){
    double n;
    double d_a;
    double d_b;
    for(i = arg->r_from; i < arg->r_to; i++){
      for(k = 0; k < arg->m2->row; k++){
        n = 0.f; d_a = 0.f; d_b = 0.f;
        for(j = 0; j < arg->m2->col; j++){
          n += arg->m1->data[i][j] * arg->m2->data[k][j];
          d_a += square(arg->m1->data[i][j]);
          d_b += square(arg->m2->data[k][j]);
        }
        arg->distances->data[k][i] = n/(sqrt(d_a)*sqrt(d_b));
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute any distance. Unknown distance method selected!\n");
    fflush(stderr);
    abort();
  }
  return 0;
}

void EuclideanDistance(matrix* m1, matrix* m2, matrix *distances, size_t nthreads)
{
  if(m1->col == m2->col){
    pthread_t *threads;
    dst_th_arg *args;
    size_t th;
    /* each column is a distance that correspond to m1->row */
    ResizeMatrix(distances, m2->row, m1->row);
    threads = xmalloc(sizeof(pthread_t)*nthreads);
    args = xmalloc(sizeof(dst_th_arg)*nthreads);

    size_t step = (size_t)ceil((double)m1->row/(double)nthreads);
    size_t from = 0;
    size_t to = step;

    for(th = 0; th < nthreads; th++){
      args[th].m1 = m1;
      args[th].m2 = m2;
      args[th].distances = distances;
      args[th].r_from = from;
      args[th].r_to = to;
      args[th].method = EUCLIDEAN;

      from = to;
      if(from+step > m1->row){
        to = m1->row;
      }
      else{
        to+=step;
      }
      pthread_create(&threads[th], NULL, CalcWorker, (void*) &args[th]);
    }

    for(th = 0; th < nthreads; th++){
      pthread_join(threads[th], NULL);
    }

    xfree(threads);
    xfree(args);
  }
  else{
    fprintf(stderr, "Unable to compute Euclidean Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void EuclideanDistance_ST(matrix *m1, matrix *m2, matrix *distances)
{
  size_t i;
  size_t j;
  size_t k;
  double dist;
  ResizeMatrix(distances, m2->row, m1->row);
  for(i = 0; i < m1->row; i++){
    for(k = 0; k < m2->row; k++){
      dist = 0.f;
      for(j = 0; j < m1->col; j++){
        dist += square(m1->data[i][j] - m2->data[k][j]);
      }
      distances->data[k][i] = sqrt(dist);
    }
  }
}

void SquaredEuclideanDistance(matrix *m1, matrix *m2, matrix *distances, size_t nthreads)
{
  if(m1->col == m2->col){
    pthread_t *threads;
    dst_th_arg *args;
    size_t th;
    /* each column is a distance that correspond to m1->row */
    ResizeMatrix(distances, m2->row, m1->row);
    threads = xmalloc(sizeof(pthread_t)*nthreads);
    args = xmalloc(sizeof(dst_th_arg)*nthreads);

    size_t step = (size_t)ceil((double)m1->row/(double)nthreads);
    size_t from = 0;
    size_t to = step;
    for(th = 0; th < nthreads; th++){
      args[th].m1 = m1;
      args[th].m2 = m2;
      args[th].distances = distances;
      args[th].r_from = from;
      args[th].r_to = to;
      args[th].method = SQUARE_EUCLIDEAN;

      pthread_create(&threads[th], NULL, CalcWorker, (void*) &args[th]);
      from = to;
      if(from+step > m1->row){
        to = m1->row;
      }
      else{
        to+=step;
      }
    }

    /* Wait till threads are complete before main continues. Unless we  */
    /* wait we run the risk of executing an exit which will terminate   */
    /* the process and all threads before the threads have completed.   */
    for(th = 0; th < nthreads; th++){
      pthread_join(threads[th], NULL);
    }

    xfree(threads);
    xfree(args);
  }
  else{
    fprintf(stderr, "Unable to compute Euclidean Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void SquaredEuclideanDistance_ST(matrix *m1, matrix *m2, matrix *distances){
  size_t i, j, k;
  double dist;

  ResizeMatrix(distances, m2->row, m1->row);

  for(i = 0; i < m1->row; i++){
    for(k = 0; k < m2->row; k++){
      dist = 0.f;
      for(j = 0; j < m1->col; j++){
        dist += square(m1->data[i][j] - m2->data[k][j]);
      }
      distances->data[k][i] = dist;
    }
  }
}

void ManhattanDistance(matrix *m1, matrix *m2, matrix *distances, size_t nthreads)
{
  if(m1->col == m2->col){
    pthread_t *threads;
    dst_th_arg *args;
    size_t th;
    /* each column is a distance that correspond to m1->row */
    ResizeMatrix(distances, m2->row, m1->row);
    threads = xmalloc(sizeof(pthread_t)*nthreads);
    args = xmalloc(sizeof(dst_th_arg)*nthreads);

    size_t step = (size_t)ceil((double)m1->row/(double)nthreads);
    size_t from = 0;
    size_t to = step;
    for(th = 0; th < nthreads; th++){
      args[th].m1 = m1;
      args[th].m2 = m2;
      args[th].distances = distances;
      args[th].r_from = from;
      args[th].r_to = to;
      args[th].method = MANHATTAN;

      pthread_create(&threads[th], NULL, CalcWorker, (void*) &args[th]);
      from = to;
      if(from+step > m1->row){
        to = m1->row;
      }
      else{
        to+=step;
      }
    }

    /* Wait till threads are complete before main continues. Unless we  */
    /* wait we run the risk of executing an exit which will terminate   */
    /* the process and all threads before the threads have completed.   */
    for(th = 0; th < nthreads; th++){
      pthread_join(threads[th], NULL);
    }

    xfree(threads);
    xfree(args);
  }
  else{
    fprintf(stderr, "Unable to compute Manhattan Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void ManhattanDistance_ST(matrix *m1, matrix *m2, matrix *distances)
{
  size_t i, j, k;
  double dist;

  ResizeMatrix(distances, m2->row, m1->row);

  for(i = 0; i < m1->row; i++){
    for(k = 0; k < m2->row; k++){
      dist = 0.f;
      for(j = 0; j < m1->col; j++){
        dist += fabs(m1->data[i][j] - m2->data[k][j]);
      }
      distances->data[k][i] = dist;
    }
  }
}

void CosineDistance(matrix *m1, matrix *m2, matrix *distances, size_t nthreads)
{
  if(m1->col == m2->col){
    pthread_t *threads;
    dst_th_arg *args;
    size_t th;
    /* each column is a distance that correspond to m1->row */
    ResizeMatrix(distances, m2->row, m1->row);
    threads = xmalloc(sizeof(pthread_t)*nthreads);
    args = xmalloc(sizeof(dst_th_arg)*nthreads);

    size_t step = (size_t)ceil((double)m1->row/(double)nthreads);
    size_t from = 0;
    size_t to = step;
    for(th = 0; th < nthreads; th++){
      args[th].m1 = m1;
      args[th].m2 = m2;
      args[th].distances = distances;
      args[th].r_from = from;
      args[th].r_to = to;
      args[th].method = COSINE;

      pthread_create(&threads[th], NULL, CalcWorker, (void*) &args[th]);
      from = to;
      if(from+step > m1->row){
        to = m1->row;
      }
      else{
        to+=step;
      }
    }

    /* Wait till threads are complete before main continues. Unless we  */
    /* wait we run the risk of executing an exit which will terminate   */
    /* the process and all threads before the threads have completed.   */
    for(th = 0; th < nthreads; th++){
      pthread_join(threads[th], NULL);
    }

    xfree(threads);
    xfree(args);
  }
  else{
    fprintf(stderr, "Unable to compute Cosine Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void CosineDistance_ST(matrix *m1, matrix *m2, matrix *distances)
{
  size_t i;
  size_t j;
  size_t k;
  double n, d_a, d_b;
  ResizeMatrix(distances, m2->row, m1->row);
  for(i = 0; i < m1->row; i++){
    for(k = 0; k < m2->row; k++){
      n = 0.f; d_a = 0.f; d_b = 0.f;
      for(j = 0; j < m1->col; j++){
        n += m1->data[i][j] * m2->data[k][j];
        d_a += square(m1->data[i][j]);
        d_b += square(m2->data[k][j]);
      }
      distances->data[k][i] = n/(sqrt(d_a)*sqrt(d_b));
    }
  }
}

double MatrixMatrixDistance(matrix *m1, matrix *m2)
{
  double dist = MISSING;
  size_t i;
  size_t j;
  size_t k;
  size_t n;
  if(m1->col == m2->col){
    dist = 0.f;
    n = 0;
    for(j = 0; j < m1->col; j++){
      for(i = 0; i < m1->row; i++){
        for(k = 0; k < m2->row; k++){
          dist += square(m1->data[i][j] - m2->data[k][j]);
          n++;
        }
      }
    }
    dist /= (double)n;
    dist = sqrt(dist);
  }
  else{
    fprintf(stderr, "Unable to compute MatrixMahalanobisDistance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
  return dist;
}

void CovarianceDistanceMap(matrix* mi, matrix *mo)
{
  size_t i, j;
  dvector *x;
  dvector *p;
  dvector *colavg;
  matrix *covmx;
  matrix *inv_cov;
  
  initDVector(&colavg);
  MatrixColAverage(mi, colavg);
  
  initMatrix(&covmx);
  MatrixCovariance(mi, covmx);

  NewMatrix(&inv_cov, mi->col, mi->col);
  MatrixInversion(covmx, inv_cov);
  DelMatrix(&covmx);

  ResizeMatrix(mo, mi->row, mi->col);

  NewDVector(&x, mi->row);
  NewDVector(&p, mi->row);

  for(j = 0; j < mi->col; j++){
    for(i = 0; i < mi->row; i++){
      x->data[i] = mi->data[i][j]-colavg->data[j];
    }

    MT_DVectorMatrixDotProduct(inv_cov, x, p);

    for(i = 0; i < mi->col; i++){
      mo->data[i][j] = p->data[i];
    }

    DVectorSet(p, 0.f);
  }

  DelDVector(&p);
  DelDVector(&x);
  DelMatrix(&inv_cov);
  DelDVector(&colavg);
}

size_t square_to_condensed_index(size_t i, size_t j, size_t n)
{
  if(i == j)
    return n+1;
  else{
    size_t ii;
    size_t jj;
    if(i < j){
      ii = j;
      jj = i;
    }
    else{
      ii = i;
      jj = j;
    }
    return n*jj - jj*(jj+1)/2 + ii - 1 - jj;
  }
}

typedef struct{
  matrix *m;
  dvector *distances;
  size_t r_from, r_to;
  enum cmethod method;
} cdst_th_arg;


void *CalcCondensedWorker(void *arg_)
{
  size_t i;
  size_t j;
  size_t k;
  cdst_th_arg *arg = (cdst_th_arg*) arg_;
  if(arg->method == EUCLIDEAN ||
     arg->method == SQUARE_EUCLIDEAN){
    double dist;
    for(i = arg->r_from; i < arg->r_to; i++){
      for(k = i+1; k < arg->m->row; k++){
        dist = 0.f;
        for(j = 0; j < arg->m->col; j++){
          dist += square(arg->m->data[i][j] - arg->m->data[k][j]);
        }
        size_t indx = square_to_condensed_index(i, k, arg->m->row);
        if(arg->method == EUCLIDEAN){
          arg->distances->data[indx] = sqrt(dist);
        }
        else{
          arg->distances->data[indx] = dist;
        }
      }
    }
  }
  else if(arg->method == MANHATTAN){
    double dist;
    for(i = arg->r_from; i < arg->r_to; i++){
      for(k = i+1; k < arg->m->row; k++){
        dist = 0.f;
        for(j = 0; j < arg->m->col; j++){
          dist += fabs(arg->m->data[i][j] - arg->m->data[k][j]);
        }
        size_t indx = square_to_condensed_index(i, k, arg->m->row);
        arg->distances->data[indx] = dist;
      }
    }
  }
  else if(arg->method == COSINE){
    double n;
    double d_a;
    double d_b;
    for(i = arg->r_from; i < arg->r_to; i++){
      for(k = i+1; k < arg->m->row; k++){
        n = 0.f; d_a = 0.f; d_b = 0.f;
        for(j = 0; j < arg->m->col; j++){
          n += arg->m->data[i][j] * arg->m->data[k][j];
          d_a += square(arg->m->data[i][j]);
          d_b += square(arg->m->data[k][j]);
        }
        size_t indx = square_to_condensed_index(i, k, arg->m->row);
        arg->distances->data[indx] = n/(sqrt(d_a)*sqrt(d_b));
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute any distance. Unknown distance method selected!\n");
    fflush(stderr);
    abort();
  }
  return 0;
}

void EuclideanDistanceCondensed(matrix* m, dvector *distances, size_t nthreads)
{
  pthread_t *threads;
  cdst_th_arg *args;
  size_t th;

  /*
   * m->row*m->row is the normal square matrix
   * with -m->row you remove the diagonal
   * dividing the result by 2 you get the half part of
   * the square matrix!
   */
  DVectorResize(distances, ((m->row*m->row)-m->row)/2);

  threads = xmalloc(sizeof(pthread_t)*nthreads);
  args = xmalloc(sizeof(cdst_th_arg)*nthreads);

  size_t step = (size_t)ceil((double)m->row/(double)nthreads);
  size_t from = 0;
  size_t to = step;


  for(th = 0; th < nthreads; th++){
    args[th].m = m;
    args[th].distances = distances;
    args[th].r_from = from;
    args[th].r_to = to;
    args[th].method = EUCLIDEAN;

    from = to;
    if(from+step > m->row){
      to = m->row;
    }
    else{
      to+=step;
    }
    pthread_create(&threads[th], NULL, CalcCondensedWorker, (void*) &args[th]);
  }

  for(th = 0; th < nthreads; th++){
    pthread_join(threads[th], NULL);
  }

  xfree(threads);
  xfree(args);
}

/* Description: calculate the square euclidean distance in a condensed way */
void SquaredEuclideanDistanceCondensed(matrix *m, dvector *distances, size_t nthreads)
{
  pthread_t *threads;
  cdst_th_arg *args;
  size_t th;

  /*
   * m->row*m->row is the normal square matrix
   * with -m->row you remove the diagonal
   * dividing the result by 2 you get the half part of
   * the square matrix!
   */
  DVectorResize(distances, ((m->row*m->row)-m->row)/2);

  threads = xmalloc(sizeof(pthread_t)*nthreads);
  args = xmalloc(sizeof(cdst_th_arg)*nthreads);

  size_t step = (size_t)ceil((double)m->row/(double)nthreads);
  size_t from = 0;
  size_t to = step;


  for(th = 0; th < nthreads; th++){
    args[th].m = m;
    args[th].distances = distances;
    args[th].r_from = from;
    args[th].r_to = to;
    args[th].method = SQUARE_EUCLIDEAN;

    from = to;
    if(from+step > m->row){
      to = m->row;
    }
    else{
      to+=step;
    }
    pthread_create(&threads[th], NULL, CalcCondensedWorker, (void*) &args[th]);
  }

  for(th = 0; th < nthreads; th++){
    pthread_join(threads[th], NULL);
  }

  xfree(threads);
  xfree(args);
}

/* Description: calculate the manhattan distance in a condensed way */
void ManhattanDistanceCondensed(matrix *m, dvector *distances, size_t nthreads)
{
  pthread_t *threads;
  cdst_th_arg *args;
  size_t th;

  /*
   * m->row*m->row is the normal square matrix
   * with -m->row you remove the diagonal
   * dividing the result by 2 you get the half part of
   * the square matrix!
   */
  DVectorResize(distances, ((m->row*m->row)-m->row)/2);

  threads = xmalloc(sizeof(pthread_t)*nthreads);
  args = xmalloc(sizeof(cdst_th_arg)*nthreads);

  size_t step = (size_t)ceil((double)m->row/(double)nthreads);
  size_t from = 0;
  size_t to = step;


  for(th = 0; th < nthreads; th++){
    args[th].m = m;
    args[th].distances = distances;
    args[th].r_from = from;
    args[th].r_to = to;
    args[th].method = MANHATTAN;

    from = to;
    if(from+step > m->row){
      to = m->row;
    }
    else{
      to+=step;
    }
    pthread_create(&threads[th], NULL, CalcCondensedWorker, (void*) &args[th]);
  }

  for(th = 0; th < nthreads; th++){
    pthread_join(threads[th], NULL);
  }

  xfree(threads);
  xfree(args);
}

/* Description: calculate the cosine distance in a condensed way */
void CosineDistanceCondensed(matrix *m, dvector *distances, size_t nthreads)
{
  pthread_t *threads;
  cdst_th_arg *args;
  size_t th;

  /*
   * m->row*m->row is the normal square matrix
   * with -m->row you remove the diagonal
   * dividing the result by 2 you get the half part of
   * the square matrix!
   */
  DVectorResize(distances, ((m->row*m->row)-m->row)/2);

  threads = xmalloc(sizeof(pthread_t)*nthreads);
  args = xmalloc(sizeof(cdst_th_arg)*nthreads);

  size_t step = (size_t)ceil((double)m->row/(double)nthreads);
  size_t from = 0;
  size_t to = step;


  for(th = 0; th < nthreads; th++){
    args[th].m = m;
    args[th].distances = distances;
    args[th].r_from = from;
    args[th].r_to = to;
    args[th].method = COSINE;

    from = to;
    if(from+step > m->row){
      to = m->row;
    }
    else{
      to+=step;
    }
    pthread_create(&threads[th], NULL, CalcCondensedWorker, (void*) &args[th]);
  }

  for(th = 0; th < nthreads; th++){
    pthread_join(threads[th], NULL);
  }

  xfree(threads);
  xfree(args);
}

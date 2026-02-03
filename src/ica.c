/* ica.c
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
#include <math.h>

#include "memwrapper.h"
#include "numeric.h"
#include "ica.h"
#include "pca.h"
#include "vector.h"
#include "matrix.h"
#include "preprocessing.h"
#include "scientificinfo.h"


void NewICAModel(ICAMODEL **m)
{
  (*m) = xmalloc(sizeof(ICAMODEL));
  initMatrix(&((*m)->S));
  initMatrix(&((*m)->W));
  initMatrix(&((*m)->whitening_matrix));
  initDVector(&((*m)->colaverage));
  initDVector(&((*m)->colscaling));
}

void DelICAModel(ICAMODEL** m)
{
  if (*m) {
    DelDVector(&((*m)->colscaling));
    DelDVector(&((*m)->colaverage));
    DelMatrix(&((*m)->whitening_matrix));
    DelMatrix(&((*m)->W));
    DelMatrix(&((*m)->S));
    xfree((*m));
    *m = NULL;
  }
}

void g(dvector *x, double alpha)
{
  size_t i;
  for(i = 0; i < x->size; i++){
    x->data[i] =  tanh(x->data[i] * alpha);
  }
}

void g_der(dvector *x, double alpha)
{
  size_t i;
  double v;
  for(i = 0; i < x->size; i++){
    v = tanh(x->data[i] * alpha);
    x->data[i] = alpha * (1.f-(v*v));
  }
}

/*
 * new_w = (X_whitened * g(w^T*X_whitened)^T)/M - (g_der(w^T*X_whitened)*1_vector_size_M*w)/M
 */
void newW(dvector *w, matrix *X, dvector *new_w, double alpha)
{
  size_t i;
  dvector *wX;
  dvector *wX_copy;
  dvector *XwX;
  dvector *one;
  dvector *g_der_one_w;
  /* Calculation of X * g(w^T*X)^T)/M */
  NewDVector(&wX, X->row);
  NewDVector(&wX_copy, X->row);
  NewDVector(&XwX, X->col);
  NewDVector(&one, X->row);
  DVectorSet(one, 1.f);

  /* w^T * X_whitened in Python is X * w in C because X is samples x features */
  MatrixDVectorDotProduct(X, w, wX);
  DVectorCopy(wX, wX_copy);/* copy of wX to calculate the derivative */
  g(wX, alpha);
  
  /* X * g(w^T*X)^T in Python is X^T * g(X*w) in C */
  DVectorMatrixDotProduct(X, wX, XwX);
  for(i = 0; i < XwX->size; i++){
    XwX->data[i]/=(double)X->row;
  }

  /* Calculation of (g_der(w^T*X)*one*w)/M */
  g_der(wX_copy, alpha);
  double g_der_one = DVectorDVectorDotProd(wX_copy, one);
  
  NewDVector(&g_der_one_w, w->size);
  DVectorCopy(w, g_der_one_w);
  
  for(i = 0; i < w->size; i++){
    g_der_one_w->data[i] *= (g_der_one / (double)X->row);
  }
  
  /* Calculation of new_w = (X * g(w^T*X)^T)/M - (g_der(w^T*X)*one*w)/M 
   * by using the previous calculations.
   */
  for(i = 0; i < new_w->size; i++){
    new_w->data[i] = XwX->data[i] - g_der_one_w->data[i];
  }

  DelDVector(&one);
  DelDVector(&g_der_one_w);
  DelDVector(&XwX);
  DelDVector(&wX_copy);
  DelDVector(&wX);
}

/*
 * FAST ica Algorithm:
 * - https://en.wikipedia.org/wiki/FastICA#Prewhitening_the_data
 * - doi:10.1016/j.trac.2013.03.013
 * 
 * Algorithm
 * 
 * inputs: 
 *     - mx M x N matrix
 *     - n_signals = number of source signals to extracts
 *     - alpha = parameter for contrast function g (default 1)
 *     - thresh = convergence threshold (default 1e-8)
 *     - max_iter = maximum number of iterations (default 5000)
 * outputs:
 *     - S new signal space size  M x n_signals
 *     - W unmixing matrix size M x n_signals
 * 
 */
void ICA(matrix *mx,
         size_t scaling,
         size_t n_signals,
         ICAMODEL *model)
{
  ICA_ext(mx, scaling, n_signals, 1.0, 1e-8, 5000, model);
}

void ICA_ext(matrix *mx,
         size_t scaling,
         size_t n_signals,
         double alpha,
         double thresh,
         size_t max_iter,
         ICAMODEL *model)
{
  /* Input validation */
  if (!mx || !model) {
    fprintf(stderr, "Error: Invalid input parameters to ICA\n");
    return;
  }
  
  if (mx->row == 0 || mx->col == 0 || n_signals == 0) {
    fprintf(stderr, "Error: Invalid matrix dimensions or number of signals\n");
    return;
  }
  
  if (n_signals > mx->col) {
    fprintf(stderr, "Error: Number of signals cannot exceed matrix columns\n");
    return;
  }

  size_t j;
  size_t k;
  size_t ic;

  matrix *mx_whitened; /* data matrix of autoscaled / mean centred object */
  dvector *w;
  dvector *w_new;
  dvector *wj;
  dvector *s;

  double mod;

  /* check if matrix have nan or infinite and convert them to MISSING */
  MatrixCheck(mx);

  NewMatrix(&mx_whitened, mx->row, mx->col);
  matrix *a;
  NewMatrix(&a, mx->row, mx->col);
  /* Center and scale the input matrix */
  MatrixPreprocess(mx,
                   scaling,
                   model->colaverage,
                   model->colscaling,
                   a);

  MatrixWhitening(a,
                  model->whitening_matrix,
                  mx_whitened);
  DelMatrix(&a);

  for(ic = 0; ic < n_signals; ic++){
    NewDVector(&w, mx_whitened->col);
    NewDVector(&w_new, mx_whitened->col);
    for(j = 0; j < w->size; j++){
      w->data[j] = randDouble(-1, 1);
    }
    /* Initial normalization */
    DVectNorm(w, w);
    
    size_t iter_count = 0;
    double lim = 1.0;
    
    while(lim > thresh && iter_count < max_iter){
      newW(w, mx_whitened, w_new, alpha);
      /*
       * Decorrelate weights (Gram-Schmidt)
       * w_new = w_new - Sum (w_new * wj) * wj
       */
      if(ic >= 1){
        NewDVector(&s, w_new->size);
        DVectorSet(s, 0.0);
        for(j = 0; j < model->W->col; j++){
          wj = getMatrixColumn(model->W, j);
          mod = DVectorDVectorDotProd(w_new, wj);
          for(k = 0; k < s->size; k++){
            s->data[k] += mod*wj->data[k];
          }
          DelDVector(&wj); 
        }

        for(k = 0; k < s->size; k++){
          w_new->data[k] -= s->data[k];
        }
        DelDVector(&s);
      }

      /* Normalization */
      DVectNorm(w_new, w_new);

      /* Calculate limit condition: lim = np.abs(np.abs((wNew * w).sum()) - 1) */
      double dot = DVectorDVectorDotProd(w_new, w);
      lim = fabs(fabs(dot) - 1.0);

      /* Update weights */
      DVectorCopy(w_new, w);
      
      iter_count++;
    }

    /* Append the calculated w result to model->W */
    MatrixAppendCol(model->W, w);

    DelDVector(&w_new);
    DelDVector(&w);
  }
  
  /* S = X_whitened * W */
  ResizeMatrix(model->S, mx_whitened->row, n_signals);
  MatrixDotProduct(mx_whitened, model->W, model->S);
  DelMatrix(&mx_whitened);
}


void ICASignalPredictor(matrix *mx,
                        ICAMODEL *model,
                        matrix *p_signals)
{
  /* Input validation */
  if (!mx || !model || !p_signals) {
    fprintf(stderr, "Error: Invalid input parameters to ICASignalPredictor\n");
    return;
  }
  
  if (mx->col != model->colscaling->size) {
    fprintf(stderr, "Error: Input matrix columns don't match model dimensions\n");
    return;
  }

  matrix *mx_whitened;
  initMatrix(&mx_whitened);
  MatrixPreprocess(mx,
                   -1,
                   model->colaverage,
                   model->colscaling, mx_whitened);

  MatrixWhitening(mx_whitened,
                  model->whitening_matrix,
                  mx_whitened);

  ResizeMatrix(p_signals, mx_whitened->row, model->W->col);
  MatrixDotProduct(mx_whitened, model->W, p_signals);
  DelMatrix(&mx_whitened);
}

void PrintICA(ICAMODEL *m)
{
  puts("New Signal Sources");
  PrintMatrix(m->S);
  printf("\n");

  puts("\nICA Weights");
  PrintMatrix(m->W);

  puts("\nICA whitening matrix");
  PrintMatrix(m->whitening_matrix);
}
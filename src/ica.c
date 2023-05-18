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
#include "vector.h"
#include "matrix.h"
#include "preprocessing.h"
#include "scientificinfo.h"


void NewICAModel(ICAMODEL **m)
{
  (*m) = xmalloc(sizeof(ICAMODEL));
  initMatrix(&((*m)->S));
  initMatrix(&((*m)->W));
  initDVector(&((*m)->colaverage));
  initDVector(&((*m)->colscaling));
}

void DelICAModel(ICAMODEL** m)
{
  DelDVector(&((*m)->colscaling));
  DelDVector(&((*m)->colaverage));
  DelMatrix(&((*m)->W));
  DelMatrix(&((*m)->S));
  xfree((*m));
}

void g(dvector *x)
{
  size_t i;
  for(i = 0; i < x->size; i++){
    x->data[i] =  tanh(x->data[i]);
  }
}

void g_der(dvector *x)
{
  size_t i;
  double v;
  for(i = 0; i < x->size; i++){
    v = tanh(x->data[i]);
    x->data[i] = 1.f-(v*v);
  }
}

void newW(dvector *w, matrix *X, dvector *new_w)
{
  size_t i;
  size_t j;
  dvector *a;
  NewDVector(&a, X->row);
  DVectorMatrixDotProduct(X, w, a);
  dvector *g_a;
  dvector *g_der_a;
  initDVector(&g_a);
  initDVector(&g_der_a);
  DVectorCopy(a, g_a);
  g(g_a);
  DVectorCopy(a, g_der_a);
  g_der(g_der_a);
  double mean_g_der = 0.f;
  for(i = 0; i < g_der_a->size; i++)
    mean_g_der += g_der_a->data[i];
  mean_g_der /= (double)g_der_a->size;

  dvector *b;
  NewDVector(&b, X->col);
  
  for(i = 0; i < X->row; i++){
    for(j = 0; j < X->col; j++){
      b->data[j] += X->data[i][j] * g_a->data[i];
    }
  }

  /* 
   * w_new = (X * g(np.dot(w.T, X))).mean(axis=1) - g_der(np.dot(w.T, X)).mean() * w
   *
   * (X * g(np.dot(w.T, X))).mean(axis=1)  = b
   * 
   * g_der(np.dot(w.T, X)).mean()  = mean_g_der
   * 
   */

  for(j = 0; j < X->col; j++){
    b->data[j] /= (double)X->row;
    b->data[j] -= mean_g_der;
    new_w->data[j] = b->data[j]*w->data[j];
  }
  
  double sum_sqrt = 0.f;
  for(j = 0; j < X->col; j++){
    sum_sqrt += b->data[j]*b->data[j];
  }

  sum_sqrt = sqrt(sum_sqrt);
   for(j = 0; j < X->col; j++){
    b->data[j] /= sum_sqrt;
  }

  DelDVector(&g_a);
  DelDVector(&g_der_a);
  DelDVector(&a);
}

/*
 * Algorithm:
 * - doi:10.1016/j.trac.2013.03.013
 * - https://towardsdatascience.com/independent-component-analysis-ica-in-python-a0ef0db0955e
 * 
 */
void ICA(matrix *mx,
         size_t scaling,
         size_t n_signals,
         ICAMODEL *model,
         ssignal *s)
{
 /* 1. Calculate the PCA extracting
  * 2. Rotation of loadings
  * Central Limit Theorem:
  * S = W*X
  * W is the demelange matrix
  * X = loadings from PCA
  */

  size_t i;
  size_t j;
  size_t ic;

  matrix *E; /* data matrix of autoscaled / mean centred object */
  dvector *w;
  dvector *w_new;

  NewMatrix(&E, mx->row, mx->col);

  /* check if matrix have nan or infinite and convert them to MISSING */
  MatrixCheck(mx);

  /* Center and scale the input matrix */
  MatrixPreprocess(mx,
                   scaling,
                   model->colaverage,
                   model->colscaling,
                   E);

  MatrixWhitening(E, E);
  
  ResizeMatrix(model->W, n_signals, E->col);
  ResizeMatrix(model->S, E->row, n_signals);

  for(ic = 0; ic < n_signals; ic++){
    NewDVector(&w, E->col);
    NewDVector(&w_new, E->col);
    for(j = 0; j < w->size; j++){
      w->data[j] = randDouble(-1, 1);
    }

    for(i = 0; i < 1000; i++){
      newW(w, E, w_new);
      if(ic >= 1){ /* remove the other source s */
        /*w_new -= np.dot(np.dot(w_new, W[:i].T), W[:i])*/
      }
    }    
    DelDVector(&w);
    DelDVector(&w_new);
  }
  
  /* S = np.dot(W, X) */
  
}


void PrintICA(ICAMODEL *m)
{
  puts("New Signal Sources");
  PrintMatrix(m->S);
  printf("\n");

  puts("\nICA Weights");
  PrintMatrix(m->W);
}

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
  DelDVector(&((*m)->colscaling));
  DelDVector(&((*m)->colaverage));
  DelMatrix(&((*m)->whitening_matrix));
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
#include<unistd.h>

/*
 * new_w = (X_whitened * g(w^T*X_whitened)^T)/M - (g_der(w^T*X_whitened)*1_vector_size_M*w)/M
 */
void newW(dvector *w, matrix *X, dvector *new_w)
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

  /*printf("%zu x %zu x %zu x 1 = %zu x 1 \n", X->row, X->col, w->size, wX->size);*/
  MatrixDVectorDotProduct(X, w, wX);
  DVectorCopy(wX, wX_copy);/* copy of wX to calculate the derivative */
  g(wX);
  /*printf("%zu x %zu x %zu x 1 = %zu x 1 \n", X->row, X->col, wX->size, XwX->size);*/
  DVectorMatrixDotProduct(X, wX, XwX);
  for(i = 0; i < XwX->size; i++){
    XwX->data[i]/=(double)X->row;
  }

  /* Calculation of (g_der(w^T*X)*one*w)/M */
  g_der(wX_copy);
  double g_der_one = DVectorDVectorDotProd(wX_copy, one);
  
  initDVector(&g_der_one_w);
  DVectorCopy(w, g_der_one_w);
  
  for(i = 0; i < w->size; i++){
    g_der_one_w->data[i] *= g_der_one;
    g_der_one_w->data[i] /= (double)X->row;
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
 * outputs:
 *     - S new signal space size  M x n_signals
 *     - W unmixing matrix size M x n_signals
 * 
 * center matrix to zero  mx = X
 * whitening matrix X -> X_whitened (store the whitening matrix for future predictions)
 * 
 * for p in 1 to C:
 *     w = random vector of lenght N
 *     while w changes
 *        new_w = (X_whitened * g(w^T*X_whitened)^T)/M - (g_der(w^T*X_whitened)*1*w)/M
 *        new_w = new_w - Sum (new_w - all_prev_w) * new_w
 *        new_w = new_w/||new_w||
 *        check if new_w ~ w. If yes stop else continue
 *     W.append(new_w)
 * 
 *  S = X_whitened * W 
 * 
 * 
 * 
 * 
 */
void ICA(matrix *mx,
         size_t scaling,
         size_t n_signals,
         ICAMODEL *model)
{
 /* 1. Calculate the PCA extracting
  * 2. Rotation of loadings
  * Central Limit Theorem:
  * S = W*X
  * W is the demelange matrix
  * X = loadings from PCA
  */

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
  puts("Centered matrix");
  PrintMatrix(a);

  MatrixWhitening(a,
                  model->whitening_matrix,
                  mx_whitened);
  DelMatrix(&a);

  puts("Whitened Matrix");
  printf("%zu %zu\n", model->whitening_matrix->row, model->whitening_matrix->col);
  PrintMatrix(mx_whitened);
  sleep(2);

  mx_whitened->data[0][0] = 15581708.89772898; mx_whitened->data[0][1] = 3738282.0789715727; mx_whitened->data[0][2] = 616737.7235975976;
  mx_whitened->data[1][0] = -33024973.944985524; mx_whitened->data[1][1] = -61100142.84298632; mx_whitened->data[1][2] = -68499847.63121258;
  mx_whitened->data[2][0] = -273674700.82157916; mx_whitened->data[2][1] = -303367702.828234; mx_whitened->data[2][2] = -311193815.5725274;
  mx_whitened->data[3][0] = -310021478.96733695; mx_whitened->data[3][1] = -353761425.7790786; mx_whitened->data[3][2] = -365289860.2342865;
  mx_whitened->data[4][0] = 32355196.352794625; mx_whitened->data[4][1] = 15106939.908953723; mx_whitened->data[4][2] = 10560859.556477588;
  mx_whitened->data[5][0] = 12257796.148508463; mx_whitened->data[5][1] = -2169869.944647597; mx_whitened->data[5][2] = -5972535.927271978;
  mx_whitened->data[6][0] = -242395274.77423197; mx_whitened->data[6][1] = -261796514.41807708; mx_whitened->data[6][2] = -266910052.90352032;
  mx_whitened->data[7][0] = -116787811.31648967; mx_whitened->data[7][1] = -160695662.42759553; mx_whitened->data[7][2] = -172268349.50678033;
  mx_whitened->data[8][0] = 213605536.8350418; mx_whitened->data[8][1] = 200453330.40735683; mx_whitened->data[8][2] = 196986834.20210207;
  mx_whitened->data[9][0] = 155-24815055.1608388781708; mx_whitened->data[9][1] = -57687855.3787338; mx_whitened->data[9][2] = -66352060.85964788;

  for(ic = 0; ic < n_signals; ic++){
    NewDVector(&w, mx_whitened->col);
    NewDVector(&w_new, mx_whitened->col);
    for(j = 0; j < w->size; j++){
      w->data[j] = randDouble(-1, 1);
    }
    
    do{ /* Loop untill convergence! */
      newW(w, mx_whitened, w_new);
      /*
       * Remove the other sources
       * new_w = new_w - Sum (new_w - all_prev_w) * new_w
       */
      if(ic >= 1){
        NewDVector(&s, w_new->size);
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

        /*
         * new_w normalization 
         * new_w = new_w/||new_w||
         */
        mod = DVectorDVectorDotProd(w_new, w_new);
        for(k = 0; k < s->size; k++){
          w_new->data[k] /= mod;
        }
        DelDVector(&s);
      }
      
      printf("Loop %f\n", calcConvergence(w_new, w));
      for(j = 0; j < w_new->size; j++)
        w->data[j] = w_new->data[j];
    } while(calcConvergence(w_new, w) > ICACONVERGENCE);

    /* Append the calculated w_new result to model->W
     * Every column represent an independent component
     */
    MatrixAppendCol(model->W, w_new);

    DelDVector(&w_new);
    DelDVector(&w);
  }
  
  /* S = np.dot(W, X) */
  ResizeMatrix(model->S, mx_whitened->row, n_signals);
  MatrixDotProduct(mx_whitened, model->W, model->S);
  DelMatrix(&mx_whitened);
}


void ICASignalPredictor(matrix *mx,
                        ICAMODEL *model,
                        matrix *p_signals)
{
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

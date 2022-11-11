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
#include "scientificinfo.h"


void NewICAModel(ICAMODEL **m)
{
  (*m) = xmalloc(sizeof(ICAMODEL));
  initMatrix(&((*m)->scores));
  initMatrix(&((*m)->loadings));
  initMatrix(&((*m)->dmodx));
  initDVector(&((*m)->varexp));
  initDVector(&((*m)->colaverage));
  initDVector(&((*m)->colscaling));
}

void DelICAModel(ICAMODEL** m)
{
  DelDVector(&((*m)->colscaling));
  DelDVector(&((*m)->colaverage));
  DelDVector(&((*m)->varexp));
  DelMatrix(&((*m)->dmodx));
  DelMatrix(&((*m)->loadings));
  DelMatrix(&((*m)->scores));
  xfree((*m));
}


void whitening(matrix *X, matrix *X_whiten)
{
  size_t i;
  matrix *Xcov, *D, *Dinv, *evect, *evectT, *eX, *DinveX;
  dvector *eval;
  NewMatrix(&Xcov, X->row, X->row);
  //d, E = np.linalg.eigh(cov)
  MatrixCovariance(X, Xcov);
  initDVector(&eval);
  initMatrix(&evect);
  EVectEval(Xcov, eval, evect);
  NewMatrix(&D, eval->size, eval->size);
  for(i = 0; i < eval->size; i++)
    D->data[i][i] = eval->data[i];

  NewMatrix(&Dinv, D->row, D->col);
  MatrixInversion(D, Dinv);

  NewMatrix(&evectT, evect->col, evect->row);
  MatrixTranspose(evect, evectT);

  NewMatrix(&eX, evectT->row, X->col);
  MatrixDotProduct(evectT, X, eX);

  NewMatrix(&DinveX, Dinv->row, eX->col);
  MatrixDotProduct(Dinv, eX, DinveX);

  DelMatrix(&eX);
  ResizeMatrix(X_whiten, evect->row, DinveX->col);
  MatrixDotProduct(evect, DinveX, X_whiten);

  DelMatrix(&DinveX);
  DelMatrix(&evect);
  DelMatrix(&Dinv);
  DelMatrix(&evectT);
  DelMatrix(&D);
  DelDVector(&eval);
  DelMatrix(&Xcov);
}

void newW(dvector *w, matrix *X, dvector *new_w)
{

}

/*
 * Algorithm:
 * - doi:10.1016/j.trac.2013.03.013
 * - https://towardsdatascience.com/independent-component-analysis-ica-in-python-a0ef0db0955e
 * 
 */
void ICA(matrix *mx, size_t scaling, size_t n_signals, ICAMODEL *model, ssignal *s)
{
  // 1. Calculate the PCA extracting
  // 2. Rotation of loadings
  // Central Limit Theorem:
  // S = W*X
  // W is the demelange matrix
  // X = loadings from PCA

  size_t i, j;
  double min, max;
  matrix *E; /* data matrix of autoscaled / mean centred object */
  NewMatrix(&E, mx->row, mx->col);

  /* check if matrix have nan or infinite and convert them to MISSING */
  MatrixCheck(mx);

  /*CENTERING */
  MatrixColAverage(mx, model->colaverage);
  for(j = 0; j < mx->col; j++){
    for(i = 0; i < mx->row; i++){
      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1)){
        continue;
      }
      else{
        E->data[i][j] = mx->data[i][j] - model->colaverage->data[j];
      }
    }
  }

  if(scaling > 0){
    if(scaling == 1){
      MatrixColSDEV(mx, model->colscaling);
    }
    else if(scaling == 2){
      MatrixColRMS(mx, model->colscaling);
    }
    else if(scaling == 3){ /* PARETO Autoscaling */
      MatrixColSDEV(mx, model->colscaling);
      for(i = 0; i < model->colscaling->size; i++){
        model->colscaling->data[i] = sqrt(model->colscaling->data[i]);
      }
    }
    else if(scaling == 4){ /* Range Scaling */
      for(i = 0; i < mx->col; i++){
        MatrixColumnMinMax(mx, i, &min, &max);
        DVectorAppend(model->colscaling, (max - min));
      }
    }
    else if(scaling == 5){ /* Level Scaling  */
      DVectorCopy(model->colaverage, model->colscaling);
    }
    else{
      for(i = 0; i < model->colaverage->size; i++){
        DVectorAppend(model->colscaling, 1.0);
      }
    }

    for(j = 0; j < E->col; j++){
      if(FLOAT_EQ(getDVectorValue(model->colscaling, j), 0, EPSILON)){
        for(i = 0; i< E->row; i++){
          E->data[i][j] = 0.f;
        }
      }
      else{
        for(i = 0; i < E->row; i++){
          if(FLOAT_EQ(E->data[i][j], MISSING, 1e-1)){
            continue;
          }
          else{
            E->data[i][j] /= model->colscaling->data[j];
          }
        }
      }
    }
  }

}


void PrintICA(ICAMODEL *m)
{
  size_t i, j;
  printf("Variance Explained\n");
  for(i = 0; i < m->varexp->size; i++){
    printf("PC%d: %.4f\n", (int)i+1, m->varexp->data[i]);
  }

  puts("Scores");
  for(j = 0; j < m->scores->col; j++){
    printf("   PC %d\t", (int)j);
  }
  printf("\n");

  for(i = 0; i < m->scores->row; i++){
    for(j = 0; j < m->scores->col; j++){
      printf("%8.4f\t", m->scores->data[i][j]);
    }
    printf("\n");
  }

  puts("\nLoadings");
  for(j = 0; j < m->loadings->col; j++){
    printf("   PC %d\t", (int)j);
  }
  printf("\n");

  for(i = 0; i < m->loadings->row; i++){
    for(j = 0; j < m->loadings->col; j++){
      printf("%8.4f\t", m->loadings->data[i][j]);
    }
    printf("\n");
  }

  puts("\nDmodX");
  for(j = 0; j < m->loadings->col; j++){
    printf("   PC %d\t", (int)j);
  }
  printf("\n");

  for(i = 0; i < m->dmodx->row; i++){
    for(j = 0; j < m->dmodx->col; j++){
      printf("%8.4f\t", m->dmodx->data[i][j]);
    }
    printf("\n");
  }
}

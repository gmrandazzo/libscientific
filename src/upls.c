/* upls.c
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
#include <stdbool.h>
#include <math.h>

#include "memwrapper.h"
#include "upls.h"
#include "upca.h"
#include "pca.h"
#include "preprocessing.h"
#include "matrix.h"
#include "vector.h"
#include "numeric.h"

void NewUPLSModel(UPLSMODEL** m)
{
  (*m) = xmalloc(sizeof(UPLSMODEL));
  initMatrix(&(*m)->xscores);
  initTensor(&(*m)->xloadings);
  initTensor(&(*m)->xweights);
  initDVector(&(*m)->xvarexp);
  initDVectorList(&(*m)->xcolaverage);
  initDVectorList(&(*m)->xcolscaling);

  initMatrix(&(*m)->yscores);
  initTensor(&(*m)->yloadings);
  initDVectorList(&(*m)->ycolaverage);
  initDVectorList(&(*m)->ycolscaling);

  initDVector(&(*m)->b);

  initDVector(&(*m)->r2x_model);
  initDVector(&(*m)->r2x_validation);
  initTensor(&(*m)->r2y_model);
  initTensor(&(*m)->r2y_validation);
  initTensor(&(*m)->q2y);
  initTensor(&(*m)->sdep);
  initTensor(&(*m)->sdec);
  initTensor(&(*m)->recalculated_y);
  initTensor(&(*m)->predicted_y);
  initTensor(&(*m)->recalc_residuals);
  initTensor(&(*m)->pred_residuals);

  initTensor(&(*m)->q2y_yscrambling);
  initTensor(&(*m)->sdep_yscrambling);
}

void DelUPLSModel(UPLSMODEL** m)
{
  DelMatrix(&(*m)->xscores);
  DelTensor(&(*m)->xloadings);
  DelTensor(&(*m)->xweights);
  DelDVector(&(*m)->xvarexp);
  DelDVectorList(&(*m)->xcolaverage);
  DelDVectorList(&(*m)->xcolscaling);

  DelMatrix(&(*m)->yscores);
  DelTensor(&(*m)->yloadings);
  DelDVectorList(&(*m)->ycolaverage);
  DelDVectorList(&(*m)->ycolscaling);

  DelDVector(&(*m)->b);

  DelDVector(&(*m)->r2x_model);
  DelDVector(&(*m)->r2x_validation);
  DelTensor(&(*m)->r2y_model);
  DelTensor(&(*m)->r2y_validation);
  DelTensor(&(*m)->q2y);
  DelTensor(&(*m)->sdep);
  DelTensor(&(*m)->sdec);
  DelTensor(&(*m)->recalculated_y);
  DelTensor(&(*m)->predicted_y);

  DelTensor(&(*m)->recalc_residuals);
  DelTensor(&(*m)->pred_residuals);

  DelTensor(&(*m)->q2y_yscrambling);
  DelTensor(&(*m)->sdep_yscrambling);

  xfree((*m));
}

int CheckPredTensors(tensor *X_, tensor *Y_)
{
  size_t i, j;
  for(i = 0; i < X_->order; i++){
    for(j = 0; j < Y_->order; j++)
      if(Y_->m[j]->col != X_->order || Y_->m[j]->row != X_->m[i]->col)
        return -1;
      else
        continue;
  }
  return 0;
}

int CheckTensors(tensor *X_, tensor *Y_)
{
  size_t i, j;
  for(i = 0; i < X_->order; i++){
    for(j = i+1; j < X_->order; j++){
      if(X_->m[i]->row == X_->m[j]->row
        && X_->m[i]->col == X_->m[j]->col)
        continue;
      else
        return -1;
    }

    for(j = 0; j < Y_->order; j++){
      if(X_->m[i]->row == Y_->m[j]->row)
        continue;
      else
        return -1;
    }
  }
  return 0;
}


/*
 * PCA Nipals algorithm for Multi-Way Principal Component Analysis
 * Journal of Chemometrics, Vol 1, 41-56 (1987)
 *
 * Svante Woold, Paul Geladi and Kim Esbensen
 *
 * X_ = input Tensor for object descriptors
 * Y_ = input Tensor for Y response value
 * X = tensor MeanCentered and Autoscaled
 * Y = tensor MeanCentered and Autoscaled
 * P = Loading Matrix for X
 * W = Weights Matrix
 * Q = Loading Matrix for Y
 *
 *
 * step 1: Preprocessing data the X_ and Y_ are scaled to unit variance
 *
 * for each component:
 *
 * step 2: u-start = column in Y with largest variance.
 *
 * step 3: W = u'X/u'u
 *
 * step 4: W = W / ||W|| normalize W to 1
 *
 * step 5: NIPALS shunt in order to decompose W into say sr'
 *
 *    step 5.1: s-start is the largest column in W
 *    for h = 0 to h = 2 do
 *
 *    step 5.1: r' = s'W/s's
 *
 *    step 5.2: r = r / ||r||   normalize r to 1
 *
 *    step 5.3: s = Wr
 *
 *    step 5.4: if h < 2 then bac to step 5.1 else stop and set W = sr'
 *
 *
 * step 6:  t = X..W / ||W||^2
 *
 * step 7: d = (t_new - t_old)'(t_new - t_old)/(N t_new't_new)
 *
 *  if d > 1*10^-10 return to step 3 else go to step 8
 *
 * step 8: Q = t'Y/t't
 *
 * step 9: u = Y..Q/||Q||^2
 *
 * step 10: XLoadings => P = t'X/t't
 *
 * step 11: b = t'u/t't
 *
 * step 12: X = X - t (x) P     ;  Y = Y - t (x)  Q * b   N.B.: (x) = Scalar Product
 *
 *  return to step 2 until all the component are calculated
 *
 *
 */

void UPLS(tensor* X_,
          tensor* Y_,
          size_t npc,
          size_t xautoscaling,
          size_t yautoscaling,
          UPLSMODEL* m,
          ssignal *s)
{

  if(CheckTensors(X_, Y_) == 0 && npc >= 1){
    size_t pc, i, j, k, order, column;
    double a, ssx;
    tensor *X;
    dvector *t_old; /* row vector with size X->m->row */
    dvector *t_new;
    dvector *xeval; /* t't */
    matrix *P;
    matrix *W;

    tensor *Y;
    dvector *colvar;
    dvector *u;
    dvector *yeval;
    matrix *Q;

    /* nipals shunt */
    dvector *s;
    dvector *r;

    NewTensor(&X, X_->order);
    for(k = 0; k < X_->order; k++){
      NewTensorMatrix(X, k, X_->m[k]->row, X_->m[k]->col);
    }

    TensorPreprocess(X_, xautoscaling, m->xcolaverage, m->xcolscaling, X);

    if(npc > X->m[0]->col)
      npc = X->m[0]->col;


    ssx = 0.f;
    for(k = 0; k < X->order; k++)
      for(i = 0; i < X->m[k]->row; i++)
        for(j = 0; j < X->m[k]->col; j++)
          ssx += square(getTensorValue(X, k, i, j));

    #ifdef DEBUG
    puts("Processed X Tensor");
    PrintTensor(X);
    #endif

    NewDVector(&t_old, X->m[0]->row); /* the object number is the same for all the matrix in tensor */
    NewDVector(&t_new, X->m[0]->row); /* the object number is the same for all the matrix in tensor */
    NewMatrix(&W, X->m[0]->col, X->order);
    NewMatrix(&P, X->m[0]->col, X->order);
    NewDVector(&xeval, npc);


    NewTensor(&Y, Y_->order);
    for(k = 0; k < Y_->order; k++){
      NewTensorMatrix(Y, k, Y_->m[k]->row, Y_->m[k]->col);
    }

    TensorPreprocess(Y_, yautoscaling, m->ycolaverage, m->ycolscaling, Y);
    /*END Disable */

    #ifdef DEBUG
    puts("Processed Y Tensor");
    PrintTensor(Y);
    #endif

    NewDVector(&u, Y->m[0]->row);
    NewMatrix(&Q, Y->m[0]->col, Y->order);
    NewDVector(&yeval, npc);


    for(pc = 0; pc < npc; pc++){

      /* step 2: chosing u-start with the best variance in Y*/
      order = column = 0;
      for(i = 0; i < Y->order; i++){
        initDVector(&colvar);
        MatrixColVar(Y->m[i], colvar);

        #ifdef DEBUG
        puts("Column Variance for chosing the best t-start column vector");
        PrintDVector(colvar);
        #endif

        /* check the column in tensor with the largest column variance */
        k = 0;
        for(j = 1; j < colvar->size; j++){
          if(getDVectorValue(colvar, j) > getDVectorValue(colvar, k)){
            k = j;
          }
          else{
            continue;
          }
        }

        if(i > 0){
          if(getDVectorValue(colvar, k) > a){
            a = getDVectorValue(colvar, k);
            order = i;
            column = k;
          }
          else{
            DelDVector(&colvar);
            continue;
          }
        }
        else{
          a = getDVectorValue(colvar, k);
          order = i;
          column = k;
        }

        DelDVector(&colvar);
      }


      /* copy the best column vector into u */
      for(i = 0; i < Y->m[order]->row; i++){
        setDVectorValue(u, i, getTensorValue(Y, order, i, column));
      }

      #ifdef DEBUG
      puts("The best u-start vector");
      PrintDVector(u);
      #endif

      while(1){
        /*step 3: W = u'X/u'u */
        MatrixSet(W, 0.f);
        DvectorTensorDotProduct(X, u, W);

        a = DVectorDVectorDotProd(u, u);

        for(i = 0; i < W->row; i++){
          for(j = 0; j < W->col; j++){
            setMatrixValue(W, i, j, getMatrixValue(W, i, j) / a );
          }
        }


        /*step 4: W = W / ||W|| */
        MatrixNorm(W, W);

        #ifdef DEBUG
        puts("Weight Matrix");
        PrintMatrix(W);
        #endif

        /* step 5:
        * step 5: NIPALS shunt in order to decompose W into say sr'
        *
        *    step 5.1: s-start is the largest column in W
        *    for h = 0 to h = 2 do
        *
        *    step 5.1: r' = s'W/s's
        *
        *    step 5.2: r = r / ||r||   normalize r to 1
        *
        *    step 5.3: s = Wr
        *
        *    step 5.4: if h < 2 then bac to step 5.1 else stop and set W = sr'
        */

        #ifdef DEBUG
        puts("Nipals SHUNT");
        #endif
        NewDVector(&r, W->col);
        NewDVector(&s, W->row);

        /*s-start is the larget column in W*/
        initDVector(&colvar);
        MatrixColAverage(W, colvar);

        #ifdef DEBUG
        puts("Column Variance for select objects...");
        PrintDVector(colvar);
        #endif

        j = 0;
        for(i = 1; i < colvar->size; i++){
          if(getDVectorValue(colvar, i) > getDVectorValue(colvar, j)){
            j = i;
          }
          else{
            continue;
          }
        }
        DelDVector(&colvar);

        /*copy the best column vector from P to u*/
        for(i = 0; i < W->row; i++){
          setDVectorValue(s, i, getMatrixValue(W, i, j));
        }

        #ifdef DEBUG
        puts("S-start");
        PrintDVector(s);
        puts("........");
        #endif

        for(i = 0; i < 2; i++){ /* for the simplest case */
          /*  r' = s'W/s's  */
          DVectorSet(r, 0.f);
          DVectorMatrixDotProduct(W, s, r);
          a = DVectorDVectorDotProd(s, s);

          for(j = 0; j < r->size; j++){
            setDVectorValue(r, j, getDVectorValue(r, j) / a);
          }

          /*||r|| = 1; */
          DVectNorm(r, r);

          #ifdef DEBUG
          puts("r vector");
          PrintDVector(r);
          #endif

          /*s = Wr */
          DVectorSet(s, 0.f); /* reset the s vector*/
          MatrixDVectorDotProduct(W, r, s);

          #ifdef DEBUG
          puts("s vector");
          PrintDVector(s);
          #endif
        }

        MatrixSet(W, 0.f);
        RowColOuterProduct(s, r, W);

        DelDVector(&s);
        DelDVector(&r);

        #ifdef DEBUG
        puts("Recomputed W Loadings Matrix from sr'");
        PrintMatrix(W);
        #endif

        /* step 6:  t = X..W/||W||^2*/
        DVectorSet(t_new, 0.f);
        TensorMatrixDotProduct(X, W, t_new);

        a = Matrixnorm(W);
        a = square(a);

        for(i = 0; i < t_new->size; i++)
          setDVectorValue(t_new, i, getDVectorValue(t_new, i)/a);


        /* step 7: d = (t_new - t_old)'(t_new - t_old)/(N t_new't_new)
        *
        *  if d > 1*10^-10 return to step 3 else go to step 8
        */


        #ifdef DEBUG
        puts("new t score vector");
        PrintDVector(t_new);
        #endif

        /*Check Distance */
        if(calcConvergence(t_new, t_old)){
          /* storing t score result */
          MatrixAppendCol(m->xscores, t_new);
          TensorAppendMatrix(m->xweights, W);
          setDVectorValue(xeval, pc, DVectorDVectorDotProd(t_new, t_new));
          break;
        }
        else if(_isnan_(DVectorDVectorDotProd(t_new, t_new))){
          fprintf(stderr, "UPLS Error! The Solver Engine was Unable to Converge! Please Check your data.\n");
          fflush(stderr);
          return;
          /*abort();*/
        }
        else{
          DVectorCopy(t_new, t_old);
        }
      }

      /* step 8: Q = t'Y/t't */
      MatrixSet(Q, 0.f);
      DvectorTensorDotProduct(Y, t_new, Q);

      a = DVectorDVectorDotProd(t_new, t_new);

      for(i = 0; i < Q->row; i++)
        for(j = 0; j < Q->col; j++)
          setMatrixValue(Q, i, j, getMatrixValue(Q, i, j) / a);

      /*storing Q Loadings matrix */
      TensorAppendMatrix(m->yloadings, Q);

      /* step 9: u = Y..Q/||Q||^2 */
      DVectorSet(u, 0.f);
      TensorMatrixDotProduct(Y, Q, u);

      a = Matrixnorm(Q);
      a = square(a);

      for(i = 0; i < u->size; i++)
        setDVectorValue(u, i, getDVectorValue(u, i)/a);


      /* storing u score result */
      MatrixAppendCol(m->yscores, u);


      /* step 10: XLoadings => P = t'X/t't */
      DvectorTensorDotProduct(X, t_new, P);

      a = DVectorDVectorDotProd(t_new, t_new);

      for(i = 0; i < P->row; i++)
        for(j = 0; j < P->col; j++)
          setMatrixValue(P, i, j, getMatrixValue(P, i, j) / a);

      /*storing P Loadings matrix */
      TensorAppendMatrix(m->xloadings, P);

      /* step 11: b = t'u/t't */
      DVectorAppend((m->b), (DVectorDVectorDotProd(u, t_new)) / a);

      /* step 12: X = X - t (x) P     ;  Y = Y - t (x)  Q * b   N.B.: (x) = Scalar Product */
      /* remove computed component from X */
      for(k = 0; k < X->order; k++){
        for(i = 0; i < X->m[k]->row; i++){
          for(j = 0; j < X->m[k]->col; j++){
            setTensorValue(X, k, i, j, getTensorValue(X, k, i, j) - (getDVectorValue(t_new, i)*getMatrixValue(P, j, k)));
          }
        }
      }

      /* remove computed component from Y */
      for(k = 0; k < Y->order; k++){
        for(i = 0; i < Y->m[k]->row; i++){
          for(j = 0; j < Y->m[k]->col; j++){
            setTensorValue(Y, k, i, j, getTensorValue(Y, k, i, j) - (getDVectorValue(t_new, i)*getMatrixValue(Q, j, k)*getDVectorValue(m->b, pc)) );
          }
        }
      }

      /*Adding EigenValue for estimate the principal component variance explained */
      setDVectorValue(xeval, pc, DVectorDVectorDotProd(t_new, t_new));

      /*Reset all*/
      MatrixSet(W, 0.f);
      MatrixSet(P, 0.f);
      DVectorSet(t_new, 0.f);
      DVectorSet(t_old, 0.f);

      MatrixSet(Q, 0.f);
      DVectorSet(u, 0.f);

    }


    calcVarExpressed(ssx, xeval, m->xvarexp);

    /*Recalculated y*/
    for(k = 0; k  < Y_->order; k++){
      AddTensorMatrix(m->recalculated_y, 0, 0);
    }

    for(i = 1; i <= npc; i++){
      tensor *recalculated_y;
      initTensor(&recalculated_y);
      UPLSYPredictor(m->xscores, m, i, recalculated_y);

      for(k = 0; k  < recalculated_y->order; k++){
        for(j = 0; j < recalculated_y->m[k]->col; j++){
          dvector *v = getMatrixColumn(recalculated_y->m[k], j);
          TensorAppendColumn(m->recalculated_y, k, v);
          DelDVector(&v);
        }
      }
      DelTensor(&recalculated_y);
    }

    /*Calculate recalculated residuals */
    for(k = 0; k < Y_->order; k++){
      AddTensorMatrix(m->recalc_residuals, m->recalculated_y->m[k]->row, m->recalculated_y->m[k]->col);
      for(i = 0; i < m->recalculated_y->m[k]->row; i++){
        for(j = 0; j < m->recalculated_y->m[k]->col; j++){
          setTensorValue(m->recalc_residuals, k, i, j, getTensorValue(m->recalculated_y, k, i, j) - getTensorValue(Y_, k, i, (size_t)floor(j/npc)));
        }
      }
    }

    DelDVector(&yeval);
    DelDVector(&u);
    DelMatrix(&Q);
    DelTensor(&Y);

    DelMatrix(&W);
    DelDVector(&t_old);
    DelDVector(&t_new);
    DelDVector(&xeval);
    DelMatrix(&P);
    DelTensor(&X);
  }
  else{
    fprintf(stderr, "Unable To run Multi Way PLS Module. Please Check your input data.\n");
    fflush(stderr);
  }
}


void UPLSScorePredictor(tensor *X_,  UPLSMODEL *m, size_t npc, matrix *ptscores)
{
  if(CheckPredTensors(X_, m->xloadings) == 0){
    size_t pc, i, j , k;
    tensor *X;
    dvector *t;

    double a;

    NewTensor(&X, X_->order);
    for(k = 0; k < X_->order; k++){
      NewTensorMatrix(X, k, X_->m[k]->row, X_->m[k]->col);
    }

    NewDVector(&t, X->m[0]->row);

    /*preprocessing data using preprocess*/
    for(k = 0; k < X_->order; k++){
      for(j = 0; j < X_->m[k]->col; j++){
        /* Mean Centering */
        if(m->xcolaverage->size > 0){
          for(i = 0; i < X_->m[k]->row; i++){
            X->m[k]->data[i][j] = X_->m[k]->data[i][j] - m->xcolaverage->d[k]->data[j];
          }
        }
        else{
          for(i = 0; i < X_->m[k]->row; i++){
            X->m[k]->data[i][j] = X_->m[k]->data[i][j];
            setTensorValue(X, k, i, j, getTensorValue(X_, k, i, j));
          }
        }

        if(m->xcolscaling->size > 0){
          for(i = 0; i < X->m[k]->row; i++){
            X->m[k]->data[i][j] = X_->m[k]->data[i][j]/m->xcolscaling->d[k]->data[j];
          }
        }
      }
    }

    if(m->xloadings->order < npc)
      npc = m->xweights->order;

    for(pc = 0; pc < npc; pc++){
      DVectorSet(t, 0.f);
      TensorMatrixDotProduct(X, m->xweights->m[pc], t);

      a = Matrixnorm(m->xweights->m[pc]);
      a = square(a);

      for(i = 0; i < t->size; i++)
        setDVectorValue(t, i, getDVectorValue(t, i)/a);

      MatrixAppendCol(ptscores, t);

      for(k = 0; k < X->order; k++){
        for(i = 0; i < X->m[k]->row; i++){
          for(j = 0; j < X->m[k]->col; j++){
            setTensorValue(X, k, i, j, getTensorValue(X, k, i, j) - (getDVectorValue(t, i)*getMatrixValue(m->xloadings->m[pc], j, k)));
          }
        }
      }
    }

    DelDVector(&t);
    DelTensor(&X);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute MultiWay PLS Prediction! Please check the data to predict.\n");
    fflush(stderr);
    /*abort();*/
  }
}

/* Calculate predicted Y-values as Y = t1b1Q1 + t2b2Q2 + ...*/
void UPLSYPredictor(matrix *tscores, UPLSMODEL *m, size_t npc, tensor *py)
{
  size_t pc, i, j , k;
  /*
  yloadings is organized as:
  - order = number of principal component
  - row = number of y dependent value
  - col = number of k layer that are equal to the order number of the initial matrix mx

  for each principal component
    for each principal component
      for each object to predict
        for each y of the layer
          y = t1b1Q1 + .....

  Y Loadings
  Tensor - order: 2
  Tensor No: 1 of row: 1; col: 2
    1.821           1.339
  Tensor No: 2 of row: 1; col: 2
    3.536           1.768

  t1      1.050
  t2      -0.150



    IL GIUSTO RISULTATO
    y11 = 1.821 * 1.050 + 3.536 * -0.150
    y12 = 1.339 * 1.05 + 1.768 *-0.150

  */
  if(npc > tscores->col)
    npc = tscores->col;

  for(k = 0; k < m->yloadings->m[0]->col; k++){
    AddTensorMatrix(py, tscores->row, m->yloadings->m[0]->row );
  }

  for(pc = 0; pc < npc; pc++){
    for(k = 0; k < py->order; k++){
      for(i = 0; i <  tscores->row; i++){
        for(j = 0; j < py->m[k]->col; j++){
          /*printf("[%u][%u][%u]  %f %f %f\n", (unsigned int)k, (unsigned int)i, (unsigned int)j, getDVectorValue(m->b, pc), getMatrixValue(tscores, i, pc), getTensorValue(m->yloadings, pc, j, k));*/
          setTensorValue(py, k, i, j, getTensorValue(py, k, i, j) +  (getDVectorValue(m->b, pc) * getMatrixValue(tscores, i, pc) * getTensorValue(m->yloadings, pc, j, k)));
        }
      }
    }
  }

  /* Adjusting to real Y value, so no autoscaled and autocentered */
  if(m->ycolaverage->size > 0 && m->ycolscaling->size > 0){

    for(k = 0; k < m->ycolaverage->size; k++){
      for(j = 0; j < m->ycolaverage->size; j++){
        for(i = 0; i < py->m[k]->row; i++){
          py->m[k]->data[i][j] *= m->ycolscaling->d[k]->data[j];
        }

        for(i = 0; i < py->m[k]->row; i++){
          py->m[k]->data[i][j] += m->ycolaverage->d[k]->data[j];
        }
      }
    }
  }
}

void UPLSRSquared(tensor *X_,
                  tensor *Y_,
                  UPLSMODEL *m,
                  size_t npc,
                  dvector *r2x,
                  tensor *r2y,
                  tensor *sdec)
{
  size_t pc, i, j, k;
  tensor *recalcy;
  tensor *recalcx;
  matrix *recalcscores;
  matrix *r2y_;
  matrix *sdec_;
  double n, d;

  if(npc > m->b->size)
    npc = m->b->size;

  NewMatrix(&r2y_, Y_->m[0]->col, Y_->order);
  NewMatrix(&sdec_, Y_->m[0]->col, Y_->order);

  for(pc = 1; pc <= npc; pc++){
    initTensor(&recalcy);
    initTensor(&recalcx);
    initMatrix(&recalcscores);

    UPLSScorePredictor(X_, m, pc, recalcscores);

    UPLSYPredictor(recalcscores, m, pc, recalcy);

    UPCAIndVarPredictor(recalcscores,
                        m->xloadings,
                        m->xcolaverage,
                        m->xcolscaling,
                        pc,
                        recalcx);

    for(k = 0; k < recalcy->order; k++){
      n = d = 0.f;
      for(j = 0; j < recalcy->m[k]->col; j++){
        for(i = 0; i < recalcy->m[k]->row; i++){
          n += square(recalcy->m[k]->data[i][j] - Y_->m[k]->data[i][j]);
          d += square(Y_->m[k]->data[i][j] - m->ycolaverage->d[k]->data[j]);
        }
        setMatrixValue(sdec_, j, k, sqrt(n/Y_->m[k]->row));
        setMatrixValue(r2y_, j, k, 1 - (n/d));
      }
    }

    TensorAppendMatrix(sdec, sdec_);
    TensorAppendMatrix(r2y, r2y_);

    if(r2x != NULL){
      n = d = 0.f;
      for(k = 0; k < recalcx->order; k++){
        for(j = 0; j < recalcx->m[k]->col; j++){
          for(i = 0; i < recalcx->m[k]->row; i++){
            n += square(X_->m[k]->data[i][j] - recalcx->m[k]->data[i][j]);
            d += square(X_->m[k]->data[i][j] - m->xcolaverage->d[k]->data[j]);
          }
        }
      }

      DVectorAppend(r2x, 1 - (n/d));
    }

    MatrixSet(r2y_, 0.f);
    MatrixSet(sdec_, 0.f);

    DelTensor(&recalcy);
    DelTensor(&recalcx);
    DelMatrix(&recalcscores);
  }

  DelMatrix(&sdec_);
  DelMatrix(&r2y_);
}

void UPLSRSquared_SSErr_SStot(tensor *X_,
                              tensor *Y_,
                              UPLSMODEL *m,
                              size_t npc,
                              dvector *xss_err,
                              dvector *xss_tot,
                              tensor *yss_err,
                              tensor *yss_tot,
                              tensor *predicted_y)
{
  size_t pc, i, j, k;
  tensor *recalcy;
  tensor *recalcx;
  matrix *recalcscores;
  matrix *yss_err_;
  matrix *yss_tot_;
  dvector *tmp;

  if(npc > m->b->size){
    npc = m->b->size;
  }

  NewMatrix(&yss_err_, Y_->m[0]->col, Y_->order);
  NewMatrix(&yss_tot_, Y_->m[0]->col, Y_->order);

  if(predicted_y != NULL){
    for(i = 0; i < Y_->order; i++){
      AddTensorMatrix(predicted_y, Y_->m[i]->row, 0);
    }
  }

  for(pc = 1; pc <= npc; pc++){

    initTensor(&recalcy);
    initTensor(&recalcx);
    initMatrix(&recalcscores);

    UPLSScorePredictor(X_, m, pc, recalcscores);

    UPLSYPredictor(recalcscores, m, pc, recalcy);

    UPCAIndVarPredictor(recalcscores, m->xloadings, m->xcolaverage, m->xcolscaling, pc, recalcx);

    for(k = 0; k < recalcy->order; k++){
      double n, d;
      n = d = 0.f;
      for(j = 0; j < recalcy->m[k]->col; j++){
        for(i = 0; i < recalcy->m[k]->row; i++){
          n += square(recalcy->m[k]->data[i][j] -  Y_->m[k]->data[i][j]);
          d += square(Y_->m[k]->data[i][j] - m->ycolaverage->d[k]->data[j]);
        }
        setMatrixValue(yss_err_, j, k, n);
        setMatrixValue(yss_tot_, j, k, d);

        if(predicted_y != NULL){
          tmp = getMatrixColumn(recalcy->m[k], j);
          MatrixAppendCol(predicted_y->m[k], tmp);
          DelDVector(&tmp);
        }
      }
    }

    TensorAppendMatrix(yss_err, yss_err_);
    TensorAppendMatrix(yss_tot, yss_tot_);

    if(xss_err != NULL && xss_tot != NULL){
      double n, d;
      n = d = 0.f;
      for(k = 0; k < recalcx->order; k++){
        for(j = 0; j < recalcx->m[k]->col; j++){
          for(i = 0; i < recalcx->m[k]->row; i++){
            n += square(X_->m[k]->data[i][j] - recalcx->m[k]->data[i][j]);
            d += square(X_->m[k]->data[i][j] - m->xcolaverage->d[k]->data[j]);
          }
        }
      }

      DVectorAppend(xss_err, n);
      DVectorAppend(xss_tot, d);
    }

    MatrixSet(yss_err_, 0.f);
    MatrixSet(yss_tot_, 0.f);

    DelTensor(&recalcy);
    DelTensor(&recalcx);
    DelMatrix(&recalcscores);
  }

  DelMatrix(&yss_err_);
  DelMatrix(&yss_tot_);
}

void UPLSYScrambling(tensor *X_, tensor *Y_,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t npc, size_t block,
                        size_t valtype, size_t rgcv_group, size_t rgcv_iterations,
                        tensor **q2y, tensor **sdep, ssignal *s)
{
  size_t scrambiterations, iterations_, i, j, order, q, k, n, y_, blocksize;
  int id;
  double temp;
  matrix *sorty, *gid;
  tensor *randomY, *sorted_y_id;
  tensor *tmpq2, *tmpsdep;

  srand_(X_->order*X_->m[0]->row*X_->m[0]->col*block);
  NewTensor(&randomY, Y_->order);
  NewTensor(&sorted_y_id, Y_->order);
  for(i = 0; i < Y_->order; i++){
    NewTensorMatrix(randomY, i, Y_->m[i]->row, Y_->m[i]->col);
    NewTensorMatrix(sorted_y_id, i, Y_->m[i]->row, Y_->m[i]->col);
  }

  NewMatrix(&sorty, Y_->m[0]->row, 2);
  for(order = 0; order < Y_->order; order++){
    for(j = 0; j < Y_->m[order]->col; j++){
      for(i = 0; i < Y_->m[order]->row; i++){
        sorty->data[i][0] = Y_->m[order]->data[i][j];
        sorty->data[i][1] = i;
      }
      MatrixSort(sorty, 0);

      for(i = 0; i < Y_->m[order]->row; i++){
        sorted_y_id->m[order]->data[i][j] = sorty->data[i][1];
      }
    }
  }
  DelMatrix(&sorty);

  /*calcualte the block size for the rotate matrix*/
  blocksize = (size_t)ceil(Y_->m[0]->row/(double)block);
  blocksize += (size_t)ceil((float)((blocksize*block) - Y_->m[0]->row)/  block);

  NewMatrix(&gid, block, blocksize);
  MatrixSet(gid, -2);
  /* Crate the boxes to fill -2 means no value to fill, -1 means value to fill*/
  for(i = 0, j = 0, k = 0; i < Y_->m[0]->row; i++){
    if(j < block){
      gid->data[j][k] = -1;
      j++;
    }
    else{
      j = 0;
      k++;
      gid->data[j][k] = -1;
      j++;
    }
  }

  /*get number of iterations*/
  scrambiterations = 0;
  for(i = 0; i < gid->row; i++){
    iterations_ = 0;
    for(j = 0; j < gid->col; j++){
      if((int)gid->data[i][j] == -1){
        iterations_++;
      }
      else{
        continue;
      }
    }

    if(iterations_ > scrambiterations){
      scrambiterations = iterations_;
    }
    else{
      continue;
    }
  }

  for(order = 0; order < sorted_y_id->order; order++){
    for(y_ = 0; y_ < sorted_y_id->m[order]->col; y_++){

      /* START WITH THE ORDERED Y_*/
      k = 0;
      for(i = 0; i < gid->row; i++){
        for(j = 0; j < gid->col; j++){
          if(gid->data[i][j] >= -1){
            gid->data[i][j] = sorted_y_id->m[order]->data[k][y_];
            k++;
          }
          else{
            continue;
          }
        }
      }
      /*
      puts("GID Y_");
      PrintMatrix(gid);
      */

      iterations_ = 0;
      while(iterations_ <  scrambiterations){
        if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
          break;
        }
        else{
          /* Shuffle the Y bloks */
          for(i = 0; i < gid->row; i++){
            for(j = (gid->col-1); j > 0; j--){
              if(gid->data[i][j] > -1){
                /*then shift from this id to the first value..*/
                temp = gid->data[i][j];
                for(k = j; k > 0; k--){
                  gid->data[i][k] = gid->data[i][k-1];
                }
                gid->data[i][0] = temp;
                break;
              }
              else{
                continue;
              }
            }
          }

          /*Fill the shifted y*/
          for(q = 0; q < Y_->order; q++){
            n = 0;
            for(i = 0; i < gid->row; i++){
              for(j = 0; j < gid->col; j++){
                id = gid->data[i][j];
                if(id > -1){
                  for(k = 0; k < Y_->m[q]->col; k++){
                    randomY->m[q]->data[n][k] = Y_->m[q]->data[id][k];
                  }
                  n++;
                }
                else{
                  continue;
                }
              }
            }
          }

         /*
          puts("MX");
          PrintTensor(X_);
          puts("MY");
          PrintTensor(Y_);
          puts("RANDOMY");
          PrintTensor(randomY);
          */
          /* Calculate calculate Q2 for y predicted...*/

          initTensor(&tmpq2);
          initTensor(&tmpsdep);

          if(valtype == 0){
            UPLSLOOCV(X_, randomY, xautoscaling, yautoscaling, npc, NULL, &tmpq2, &tmpsdep, NULL, NULL, s);
          }
          else{
            UPLSRandomGroupsCV(X_, randomY, xautoscaling, yautoscaling, npc, rgcv_group, rgcv_iterations, NULL, &tmpq2, &tmpsdep, NULL, NULL, s);
          }

          if(order == 0 && iterations_ == 0 && y_ == 0){
            TensorCopy(tmpq2, q2y);
            TensorCopy(tmpsdep, sdep);
          }
          else{
            for(q = 0; q < tmpq2->order; q++){
              for(j = 0; j < tmpq2->m[q]->col; j++){
                for(i = 0; i < tmpq2->m[q]->row; i++){
                  (*q2y)->m[q]->data[i][j] += tmpq2->m[q]->data[i][j];
                  (*sdep)->m[q]->data[i][j] += tmpsdep->m[q]->data[i][j];
                }
              }
            }
          }

          DelTensor(&tmpq2);
          DelTensor(&tmpsdep);
          iterations_++;
        }
      }
    }
  }

    /*Finalize the output by dividing for the number of iterations*/
  for(q = 0; q < (*q2y)->order; q++){
    for(j = 0; j < (*q2y)->m[q]->col; j++){
      for(i = 0; i < (*q2y)->m[q]->row; i++){
        (*q2y)->m[q]->data[i][j] /= (double)(Y_->m[0]->col*scrambiterations*Y_->order);
        (*sdep)->m[q]->data[i][j] /= (double)(Y_->m[0]->col*scrambiterations*Y_->order);
      }
    }
  }
  DelMatrix(&gid);
  DelTensor(&randomY);
  DelTensor(&sorted_y_id);
}

/* DEPRECATED */
void UPLSRandomGroupsCV(tensor *X_, tensor *Y_, size_t xautoscaling, size_t yautoscaling, size_t npc, size_t group, size_t iterations, dvector **r2x, tensor **q2y, tensor **sdep, tensor **predicted_y, tensor **pred_residuals, ssignal *s)
{
  if(npc > 0 && X_->m[0]->row == Y_->m[0]->row && group > 0 && iterations > 0){
    size_t iterations_, i, j, k, n, l, g;
    matrix *gid; /* randomization and storing id for each random group into a matrix */

    /* Matrix for compute the PLS models for groups */
    tensor *subX;
    tensor *subY;

    UPLSMODEL *subm;

    /* matrix for the randomg group to predict */
    tensor *predictX;
    tensor *realY;

    dvector *predictxss_err;
    dvector *predictxss_tot;
    tensor *predictyss_err;
    tensor *predictyss_tot;
    tensor *predicty;
    uivector *predictcount;

    dvector *sumxss_err;
    dvector *sumxss_tot;
    tensor *sumyss_err;
    tensor *sumyss_tot;

    initDVector(&sumxss_err);
    initDVector(&sumxss_tot);
    initTensor(&sumyss_err);
    initTensor(&sumyss_tot);

    NewMatrix(&gid, group, (size_t)ceil(X_->m[0]->row/(double)group));

    if(npc > X_->m[0]->col){
      npc = X_->m[0]->col;
    }

    if(predicted_y != NULL){
      if(pred_residuals != NULL){
        for(k = 0; k < Y_->order; k++){
          AddTensorMatrix((*predicted_y), Y_->m[k]->row, Y_->m[k]->col*npc); /* each component have my->col ypsilon */
          AddTensorMatrix((*pred_residuals), Y_->m[k]->row, Y_->m[k]->col*npc); /* each component have my->col ypsilon */
        }
      }
      else{
        for(k = 0; k < Y_->order; k++){
          AddTensorMatrix((*predicted_y), Y_->m[k]->row, Y_->m[k]->col*npc); /* each component have my->col ypsilon */
        }
      }

      NewUIVector(&predictcount, Y_->m[0]->row);
    }

    /*srand(time(NULL));*/
    srand_(group*X_->m[0]->row*iterations);
    iterations_ = 0;
    while(iterations_ <  iterations){
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* printf("iter: %u\n", (unsigned int)iterations_); */
        /* Divide in group  all the Dataset */
        MatrixSet(gid, -1);

        /*step 1*/
        k = 0;
        for(i = 0; i <  gid->row; i++){
          for(j = 0; j <  gid->col; j++){
            do{
              n = randInt(0, X_->m[0]->row);
            } while(ValInMatrix(gid, n) == 1 && k < (X_->m[0]->row));
            if(k < X_->m[0]->row){
              setMatrixValue(gid, i, j, n);
              k++;
            }
            else
              continue;
          }
        }

        /*
        puts("Gid Matrix");
        PrintMatrix(gid);
        */

        /*step 2*/
        for(g = 0; g < gid->row; g++){ /*For aeach group */
          /* Estimate how many objects are inside the sub model without the group "g" */
          n = 0;
          for(k = 0; k < gid->row; k++){
            if(k != g){
              for(j = 0; j < gid->col; j++){
                if((int)getMatrixValue(gid, k, j) != -1)
                  n++;
                else
                  continue;
              }
            }
            else
              continue;
          }

          /*Allocate the submodel*/
          NewTensor(&subX, X_->order);
          NewTensor(&subY, Y_->order);

          for(k = 0; k < X_->order; k++){
            NewTensorMatrix(subX, k, n, X_->m[k]->col);
          }

          for(k = 0; k < Y_->order; k++){
            NewTensorMatrix(subY, k, n, Y_->m[k]->col);
          }


          /* Estimate how many objects are inside the group "g" to predict*/
          n = 0;
          for(j = 0; j < gid->col; j++){
            if((int)getMatrixValue(gid, g, j) != -1)
              n++;
            else
              continue;
          }


          /*Allocate the */
          NewTensor(&predictX, X_->order);
          NewTensor(&realY, Y_->order);

          for(k = 0; k < X_->order; k++){
            NewTensorMatrix(predictX, k, n, X_->m[k]->col);
          }

          for(k = 0; k < Y_->order; k++){
            NewTensorMatrix(realY, k, n, Y_->m[k]->col);
          }

          /* copy the submodel values */

          for(i = 0, l = 0; i < gid->row; i++){
            if(i != g){
              for(j = 0; j < gid->col; j++){
                int a =  (int)getMatrixValue(gid, i, j); /* get the row index */
                if(a != -1){
                  for(k = 0; k < subX->order; k++){
                    for(n = 0; n < X_->m[k]->col; n++){
                      setTensorValue(subX, k, l, n, getTensorValue(X_, k, a, n));
                    }
                  }

                  for(k = 0; k < subY->order; k++){
                    for(n = 0; n < Y_->m[k]->col; n++){
                      setTensorValue(subY, k, l, n, getTensorValue(Y_, k, a, n));
                    }
                  }
                  l++;
                }
                else
                  continue;
              }

            }
            else
              continue;
          }

          /* copy the objects to predict into predictmx*/
          for(j = 0, l = 0; j < gid->col; j++){
            int a = (int)getMatrixValue(gid, g, j);
            if(a != -1){
              for(k = 0; k < predictX->order; k++){
                for(n = 0; n < X_->m[k]->col; n++){
                  setTensorValue(predictX, k, l, n, getTensorValue(X_, k, a, n));
                }
              }
              /* copy the real value Y_ to realY */
              for(k = 0; k < realY->order; k++){
                for(n = 0; n < Y_->m[k]->col; n++){
                  setTensorValue(realY, k, l, n, getTensorValue(Y_, k, a, n));
                }
              }
              l++;
            }
            else
              continue;
          }
          /*
          printf("Excuded the group number %u\n", (unsigned int)g);
          puts("Sub Model\nX:");
          PrintTensor(subX);
          puts("Y:");
          PrintTensor(subY);

          puts("\n\nPredict Group\nX:");
          PrintTensor(predictX);
          puts("RealY:");
          PrintTensor(realY);
          */

          NewUPLSModel(&subm);

          UPLS(subX, subY, npc, xautoscaling, yautoscaling, subm, s);

          initDVector(&predictxss_err);
          initDVector(&predictxss_tot);
          initTensor(&predictyss_err);
          initTensor(&predictyss_tot);

          if(predicted_y != NULL){
            initTensor(&predicty);
            UPLSRSquared_SSErr_SStot(predictX, realY, subm, npc, predictxss_err, predictxss_tot, predictyss_err, predictyss_tot, predicty);
          }
          else{
            UPLSRSquared_SSErr_SStot(predictX, realY, subm, npc, predictxss_err, predictxss_tot, predictyss_err, predictyss_tot, NULL);
          }

          if(iterations_ == 0 && g == 0){
            for(i = 0; i < predictxss_err->size; i++){
              DVectorAppend(sumxss_err, getDVectorValue(predictxss_err, i));
              DVectorAppend(sumxss_tot, getDVectorValue(predictxss_tot, i));
            }

            TensorCopy(predictyss_err, &sumyss_err);
            TensorCopy(predictyss_tot, &sumyss_tot);

            TensorCopy(predictyss_err, sdep);

            l = predictX->m[0]->row; /*the number of object for the current sdep*/
            for(k = 0; k < predictyss_err->order; k++){
              for(j = 0; j < predictyss_err->m[k]->col; j++){
                for(i = 0; i < predictyss_err->m[k]->row; i++){
                  setTensorValue((*sdep), k, i, j, getTensorValue((*sdep), k, i, j)/l);
                }
              }
            }

            if(predicted_y != NULL){
              for(j = 0, l = 0; j < gid->col; j++){
                int a = (int)getMatrixValue(gid, g, j);
                if(a != -1){
                  setUIVectorValue(predictcount, a, getUIVectorValue(predictcount, a)+1);
                  for(k = 0; k < predicty->order; k++){
                    for(n = 0; n < predicty->m[k]->col; n++){
                      setTensorValue((*predicted_y), k, a, n, (getTensorValue((*predicted_y), k, a, n) + getTensorValue(predicty, k, l, n)));
                    }
                  }
                  l++;
                }
                else{
                  continue;
                }
              }
            }
          }
          else{
            for(i = 0; i < predictxss_err->size; i++){
              setDVectorValue(sumxss_err, i, getDVectorValue(sumxss_err, i) + getDVectorValue(predictxss_err, i));
              setDVectorValue(sumxss_tot, i, getDVectorValue(sumxss_tot, i) + getDVectorValue(predictxss_tot, i));
            }

            for(k = 0; k < predictyss_err->order; k++){
              for(j = 0; j < predictyss_err->m[k]->col; j++){
                for(i = 0; i < predictyss_err->m[k]->row; i++){
                  setTensorValue(sumyss_err, k, i, j, getTensorValue(sumyss_err, k, i, j) + getTensorValue(predictyss_err, k, i, j));
                  setTensorValue(sumyss_tot, k, i, j, getTensorValue(sumyss_tot, k, i, j) + getTensorValue(predictyss_tot, k, i, j));
                }
              }
            }

            l = predictX->m[0]->row; /*the number of object for the current sdep*/
            for(k = 0; k < predictyss_err->order; k++){
              for(j = 0; j < predictyss_err->m[k]->col; j++){
                for(i = 0; i < predictyss_err->m[k]->row; i++){
                  setTensorValue((*sdep), k, i, j, getTensorValue((*sdep), k, i, j) + (getTensorValue(predictyss_err, k, i, j)/l));
                }
              }
            }

            if(predicted_y != NULL){
              for(j = 0, l = 0; j < gid->col; j++){
                int a = (int)getMatrixValue(gid, g, j);
                if(a != -1){
                  setUIVectorValue(predictcount, a, getUIVectorValue(predictcount, a)+1);
                  for(k = 0; k < predicty->order; k++){
                    for(n = 0; n < predicty->m[k]->col; n++){
                      setTensorValue((*predicted_y), k, a, n, (getTensorValue((*predicted_y), k, a, n) + getTensorValue(predicty,k, l, n)));
                    }
                  }
                  l++;
                }
                else{
                  continue;
                }
              }
            }
          }

          if(predicted_y != NULL){
            DelTensor(&predicty);
          }

          DelTensor(&predictyss_err);
          DelTensor(&predictyss_tot);
          DelDVector(&predictxss_err);
          DelDVector(&predictxss_tot);

          DelUPLSModel(&subm);

          DelTensor(&realY);
          DelTensor(&predictX);

          DelTensor(&subY);
          DelTensor(&subX);
        }
        iterations_++;
      }
    }
    DelMatrix(&gid);

    /*Finalize the output by dividing for the number of iterations*/
    for(i = 0; i < sumxss_err->size; i++){
      DVectorAppend((*r2x), 1 - (getDVectorValue(sumxss_err, i)/getDVectorValue(sumxss_tot, i)));
    }

    for(k = 0; k < sumyss_err->order; k++){
      TensorAppendMatrix((*q2y), sumyss_err->m[k]);
      for(j = 0; j < sumyss_err->m[k]->col; j++){
        for(i = 0; i < sumyss_err->m[k]->row; i++){
          setTensorValue((*q2y), k, i, j, 1 - (getTensorValue((*q2y), k, i, j)/getTensorValue(sumyss_tot, k, i, j)));
        }
      }
    }

    for(k = 0; k < (*sdep)->order; k++){
      for(j = 0; j < (*sdep)->m[k]->col; j++){
        for(i = 0; i < (*sdep)->m[k]->row; i++){
          /*setTensorValue((*sdep), k, i, j, sqrt(getTensorValue((*sdep), k, i, j)/(group*iterations)));*/
          setTensorValue((*sdep), k, i, j, sqrt(getTensorValue((*sdep), k, i, j)/(iterations)));
        }
      }
    }


    if(predicted_y != NULL){
      if(pred_residuals != NULL){
        for(k = 0; k < (*predicted_y)->order; k++){
          for(j = 0; j < (*predicted_y)->m[k]->col; j++){
            for(i = 0; i < (*predicted_y)->m[k]->row; i++){
              setTensorValue((*predicted_y), k, i, j, getTensorValue((*predicted_y), k, i, j)/getUIVectorValue(predictcount, i));
              setTensorValue((*pred_residuals), k, i, j, getTensorValue((*predicted_y), k, i, j) - getTensorValue(Y_, k, i, (size_t)floor(j/npc)));
            }
          }
        }
      }
      else{
        for(k = 0; k < (*predicted_y)->order; k++){
          for(j = 0; j < (*predicted_y)->m[k]->col; j++){
            for(i = 0; i < (*predicted_y)->m[k]->row; i++){
              setTensorValue((*predicted_y), k, i, j, getTensorValue((*predicted_y), k, i, j)/getUIVectorValue(predictcount, i));
            }
          }
        }
      }

      DelUIVector(&predictcount);
    }

    DelDVector(&sumxss_err);
    DelDVector(&sumxss_tot);
    DelTensor(&sumyss_err);
    DelTensor(&sumyss_tot);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute Cross Validation!!\n");
  }
}

/* Leave One Out Validation
 *
 * 1) remove one object
 * 2) calculate the model
 * 3) predict the removed object and so the r2 and q2
 */
void UPLSLOOCV(tensor* X_, tensor* Y_, size_t xautoscaling, size_t yautoscaling, size_t npc, dvector** r2x, tensor** q2y, tensor** sdep, tensor** predicted_y, tensor **pred_residuals, ssignal *s)
{
  if(npc > 0 && X_->m[0]->row == Y_->m[0]->row){
    size_t i, j, k, l, m, model;
    /* Matrix for compute the PLS models for groups */
    tensor *subX;
    tensor *subY;
    UPLSMODEL *subm;

    tensor *predictX;
    tensor *realY;

    dvector *predictxss_err;
    dvector *predictxss_tot;
    tensor *predictyss_err;
    tensor *predictyss_tot;
    tensor *predicty;


    dvector *sumxss_err;
    dvector *sumxss_tot;
    tensor *sumyss_err;
    tensor *sumyss_tot;

    if(npc > X_->m[0]->col){
      npc = X_->m[0]->col;
    }

    if(predicted_y != NULL){
      if(pred_residuals != NULL){
        for(k = 0; k < Y_->order; k++){
          AddTensorMatrix((*predicted_y), Y_->m[k]->row, Y_->m[k]->col*npc); /* each component have my->col ypsilon */
          AddTensorMatrix((*pred_residuals), Y_->m[k]->row, Y_->m[k]->col*npc); /* each component have my->col ypsilon */
        }
      }
      else{
        for(k = 0; k < Y_->order; k++){
          AddTensorMatrix((*predicted_y), Y_->m[k]->row, Y_->m[k]->col*npc); /* each component have my->col ypsilon */
        }
      }
    }

    NewTensor(&subX, X_->order);
    NewTensor(&subY, Y_->order);
    NewTensor(&predictX, X_->order);
    NewTensor(&realY, Y_->order);

    initDVector(&sumxss_err);
    initDVector(&sumxss_tot);
    initTensor(&sumyss_err);
    initTensor(&sumyss_tot);

    for(k = 0; k < X_->order; k++){
      NewTensorMatrix(subX, k, X_->m[k]->row-1, X_->m[k]->col);
      NewTensorMatrix(predictX, k, 1, X_->m[k]->col);
    }

    for(k = 0; k < Y_->order; k++){
      NewTensorMatrix(subY, k, Y_->m[k]->row-1, Y_->m[k]->col);
      NewTensorMatrix(realY, k, 1, Y_->m[k]->col);
    }

    for(model = 0; model < X_->m[0]->row; model++){
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* copy the data into subX,subY, predictX and realY */
        for(k = 0; k < X_->order; k++){
          m = 0;
          for(j = 0; j < X_->m[k]->row; j++){
            if(j != model){
              for(l = 0; l < X_->m[k]->col; l++){
                setTensorValue(subX, k, m, l, getTensorValue(X_, k, j, l));
              }
              m++;
            }
            else{
              for(l = 0; l < X_->m[k]->col; l++){
                setTensorValue(predictX, k, 0, l, getTensorValue(X_, k, j, l));
              }
            }
          }
        }

        for(k = 0; k < Y_->order; k++){
          m = 0;
          for(j = 0; j < Y_->m[k]->row; j++){
            if(j != model){
              for(l = 0; l < Y_->m[k]->col; l++){
                setTensorValue(subY, k, m, l, getTensorValue(Y_, k, j, l));
              }
              m++;
            }
            else{
              for(l = 0; l < Y_->m[k]->col; l++){
                setTensorValue(realY, k, 0, l, getTensorValue(Y_, k, j, l));
              }
            }
          }
        }

        NewUPLSModel(&subm);

        UPLS(subX, subY, npc, xautoscaling, yautoscaling, subm, s);

        initDVector(&predictxss_err);
        initDVector(&predictxss_tot);
        initTensor(&predictyss_err);
        initTensor(&predictyss_tot);
        initTensor(&predicty);

        UPLSRSquared_SSErr_SStot(predictX, realY, subm, npc, predictxss_err, predictxss_tot, predictyss_err, predictyss_tot, predicty);

        if(model == 0){
          for(i = 0; i < predictxss_err->size; i++){
            DVectorAppend(sumxss_err, getDVectorValue(predictxss_err, i));
            DVectorAppend(sumxss_tot, getDVectorValue(predictxss_tot, i));
          }

          TensorCopy(predictyss_err, &sumyss_err);
          TensorCopy(predictyss_tot, &sumyss_tot);

          TensorCopy(predictyss_err, sdep);

          l = predictX->m[0]->row; /*the number of object for the current sdep*/
          for(k = 0; k < predictyss_err->order; k++){
            for(j = 0; j < predictyss_err->m[k]->col; j++){
              for(i = 0; i < predictyss_err->m[k]->row; i++){
                setTensorValue((*sdep), k, i, j, getTensorValue((*sdep), k, i, j)/l);
              }
            }
          }

          if(predicted_y != NULL){
            if(pred_residuals != NULL){
              for(k = 0; k < predicty->order; k++){
                /*for(i = 0; i < predicty->m[k]->row; i++){*/
                for(j = 0; j < predicty->m[k]->col; j++){
                  setTensorValue((*predicted_y), k, model, j, getTensorValue(predicty, k, 0, j));
                  setTensorValue((*pred_residuals), k, model, j, getTensorValue(predicty, k, 0, j) - getTensorValue(Y_, k, model, (size_t)floor(j/npc)));
                }
                /*}*/
              }
            }
            else{
              for(k = 0; k < predicty->order; k++){
                /*for(i = 0; i < predicty->m[k]->row; i++){*/
                for(j = 0; j < predicty->m[k]->col; j++){
                  setTensorValue((*predicted_y), k, model, j, getTensorValue(predicty, k, 0, j));
                }
                /*}*/
              }
            }
          }
        }
        else{
          for(i = 0; i < predictxss_err->size; i++){
            setDVectorValue(sumxss_err, i, getDVectorValue(sumxss_err, i) + getDVectorValue(predictxss_err, i));
            setDVectorValue(sumxss_tot, i, getDVectorValue(sumxss_tot, i) + getDVectorValue(predictxss_tot, i));
          }

          for(k = 0; k < predictyss_err->order; k++){
            for(j = 0; j < predictyss_err->m[k]->col; j++){
              for(i = 0; i < predictyss_err->m[k]->row; i++){
                setTensorValue(sumyss_err, k, i, j, getTensorValue(sumyss_err, k, i, j) + getTensorValue(predictyss_err, k, i, j));
                setTensorValue(sumyss_tot, k, i, j, getTensorValue(sumyss_tot, k, i, j) + getTensorValue(predictyss_tot, k, i, j));
              }
            }
          }

          l = predictX->m[0]->row; /*the number of object for the current sdep*/
          for(k = 0; k < predictyss_err->order; k++){
            for(j = 0; j < predictyss_err->m[k]->col; j++){
              for(i = 0; i < predictyss_err->m[k]->row; i++){
                setTensorValue((*sdep), k, i, j, getTensorValue((*sdep), k, i, j) + (getTensorValue(predictyss_err, k, i, j)/l));
              }
            }
          }

          if(predicted_y != NULL){
            if(pred_residuals != NULL){
              for(k = 0; k < predicty->order; k++){
                /*for(i = 0; i < predicty->m[k]->row; i++){*/
                for(j = 0; j < predicty->m[k]->col; j++){
                  setTensorValue((*predicted_y), k, model, j, getTensorValue(predicty, k, 0, j));
                  setTensorValue((*pred_residuals), k, model, j, getTensorValue(predicty, k, 0, j) - getTensorValue(Y_, k, model, (size_t)floor(j/npc)));
                }
                /*}*/
              }
            }
            else{
              for(k = 0; k < predicty->order; k++){
                /*for(i = 0; i < predicty->m[k]->row; i++){*/
                for(j = 0; j < predicty->m[k]->col; j++){
                  setTensorValue((*predicted_y), k, model, j, getTensorValue(predicty, k, 0, j));
                }
                /*}*/
              }
            }
          }
        }

        DelTensor(&predicty);
        DelTensor(&predictyss_err);
        DelTensor(&predictyss_tot);
        DelDVector(&predictxss_err);
        DelDVector(&predictxss_tot);
        DelUPLSModel(&subm);
      }
    }
    DelTensor(&realY);
    DelTensor(&predictX);
    DelTensor(&subX);
    DelTensor(&subY);

    /*Finalize the output by dividing for the number of iterations*/
    if(r2x != NULL){
      for(i = 0; i < sumxss_err->size; i++){
        DVectorAppend((*r2x), 1 - (getDVectorValue(sumxss_err, i)/getDVectorValue(sumxss_tot, i)));
      }
    }

    for(k = 0; k < sumyss_err->order; k++){
      TensorAppendMatrix((*q2y), sumyss_err->m[k]);
      for(j = 0; j < sumyss_err->m[k]->col; j++){
        for(i = 0; i < sumyss_err->m[k]->row; i++){
          setTensorValue((*q2y), k, i, j, 1 - (getTensorValue((*q2y), k, i, j)/getTensorValue(sumyss_tot, k, i, j)));
        }
      }
    }

    for(k = 0; k < (*sdep)->order; k++){
      for(j = 0; j < (*sdep)->m[k]->col; j++){
        for(i = 0; i < (*sdep)->m[k]->row; i++){
          setTensorValue((*sdep), k, i, j, sqrt(getTensorValue((*sdep), k, i, j)/X_->m[0]->row));
        }
      }
    }

    DelDVector(&sumxss_err);
    DelDVector(&sumxss_tot);
    DelTensor(&sumyss_err);
    DelTensor(&sumyss_tot);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute UPLS Leave One Out Validation!!\n");
  }
}

/* This function from Q^2, or R^2 return the number of components to use for make predictions.
 * This cutoff components look at the best r^2 or q^2 and pay attention if r^2 or q^2 get low value
 * and after high value..
 *
 * N.B.: 0 means the first principal component, 1 is the second principal component. etc...
 */
size_t UPLSGetPCModelCutOff(tensor *rq2)
{
  size_t i, j, k, best_pc;
  matrix *rq2mean;
  double tmp_value, best_value;

  NewMatrix(&rq2mean, rq2->order, rq2->m[0]->row);

  for(k = 0; k < rq2->order; k++){
    for(i = 0; i < rq2->m[k]->row; i++){
      for(j = 0; j < rq2->m[k]->col; j++){
        setMatrixValue(rq2mean, k, i, getMatrixValue(rq2mean, k, i) + getTensorValue(rq2, k, i, j));
      }
    }
  }

  for(i = 0; i < rq2mean->row; i++){
    for(j = 0; j < rq2mean->col; j++){
      setMatrixValue(rq2mean, i, j, getMatrixValue(rq2mean, i, j)/rq2->m[0]->col);
    }
  }

  best_pc = 0;
  for(j = 0; j < rq2mean->col; j++){
    size_t curr_pc = 0;
    best_value = getMatrixValue(rq2mean, 0, j);

    for(i = 1; i < rq2mean->row; i++){
      tmp_value = getMatrixValue(rq2mean, i, j);
      if(tmp_value > best_value || FLOAT_EQ(tmp_value, best_value, EPSILON)){
        best_value = tmp_value;
        curr_pc = i;
      }
    }

    if(best_pc == 0){
      best_pc = curr_pc;
    }
    else{
      if(curr_pc < best_pc){
        best_pc = curr_pc;
      }
      else{
        continue;
      }
    }
    
  }

  DelMatrix(&rq2mean);
  return best_pc;
}

/*
 * Variable Selection using the Particle Swarm Optimization algorithm
 *
 * Input:
 *   mx and my matrix for make model; optionally px and py for testset validation.
 *
 * Output:
 * - varselected vector: vector of varialble selected
 * - validation vector:  r2 if you have a test set or q2 if the variable selection is made into de model
 *
 *
 * Each particle is a model. Each model have is velocity and his position.
 *
 * Algoritm:
 *
 * 1) Initialize the particle
 *
 * Population = The model with all varible setted
 * P_g_best = Q^2 of the model as is with all variable or R^2 if you have a prediction!
 *
 * 2) Initialize the models particle.
 * for(i = 1 to i = Population Size):
 *   P_velocity = randomvelocity for the model
 *   P_position = randomposition. This means that the variable could be setted on or off in the model.
 *
 *   P_best = Cost of model with the actual position. So calculation of R^2 or Q^2.
 *
 *   if the cost of actual model P_best is <= to the global cost P_g_best:
 *     P_g_best = P_best
 * end
 *
 * While(R^2 or Q^2 is not good enough...):
 *   For Model P in population:
 *     P_new_velocity of the model = UpdateVelocity(P_old_velocity, P_g_best, P_best)
 *     P_new_position of the model = UpdatePosition(P_old_position, P_velocity)
 *
 *     P_actual_model = R^2 or Q^2 of the actual model setted..
 *
 *     if P_actual_model < = P_best:
 *       P_best = P_actual_model
 *       if P_best <= P_g_best:
 *         P_g_best = P_best
 *
 *
 *
 * N.B.: 0 < rand() < 1
 *
 * UpdateVelocity(P_old_velocity, P_g_best, P_best):
 *   P_new_velocity = P_old_velocity + (c1 *rand() * (P_best - P_old_position)) +
 *     (c2 * rand() * (P_g_best + P_old_position))
 *
 * UpdatePosition(P_old_position, P_new_velocity):
 *   P_new_position = P_old_position + P_new_velocity
 *
 * P_new_position is the position to set on or off so is a number between 1 and the number of variables
 *
 */

void UPLSCostPopulation(tensor *ax, tensor *ay,
                      tensor *px, tensor *py,
                      uivector *popvector,
                      size_t xautoscaling, size_t yautoscaling, size_t npc,
                      int validation_type, size_t ngroup, size_t niter,
                      double *r2, double *q2, ssignal *s)
{
  size_t i, j, k, cc, subsize, cutoff;
  double r2_, q2_;
  dvector *r2x;
  tensor *subax, *subpx, *r2y, *sdec, *q2y, *sdep;
  UPLSMODEL *m;

  subsize = 0;
  for(i = 0; i < popvector->size; i++){
    if(getUIVectorValue(popvector, i) == 1){
      subsize++;
    }
    else{
      continue;
    }
  }

  NewTensor(&subax, ax->order);
  for(k = 0; k < ax->order; k++){
    NewTensorMatrix(subax, k, ax->m[k]->row, subsize);
  }


  for(k = 0; k < ax->order; k++){
    for(i = 0; i < ax->m[k]->row; i++){
      cc = 0;
      for(j = 0; j < ax->m[k]->col; j++){
        if(getUIVectorValue(popvector, j) == 1){
          setTensorValue(subax, k, i, cc, getTensorValue(ax, k, i, j));
          cc++;
        }
        else{
          continue;
        }
      }
    }
  }


  if(px != NULL && py != NULL && validation_type == 0){ /* cost = R^2*/
    NewUPLSModel(&m);
    UPLS(subax, ay, npc, xautoscaling, yautoscaling, m, s);

    NewTensor(&subpx, px->order);
    for(k = 0; k < px->order; k++){
      NewTensorMatrix(subpx, k, px->m[k]->row, subsize);
    }

    for(k = 0; k < px->order; k++){
      for(i = 0; i < px->m[k]->row; i++){
        cc = 0;
        for(j = 0; j < px->m[k]->col; j++){
          if(getUIVectorValue(popvector, j) == 1){
            setTensorValue(subpx, k, i, cc, getTensorValue(px, k, i, j));
            cc++;
          }
          else{
            continue;
          }
        }
      }
    }

    initTensor(&r2y);
    initTensor(&sdec);
    UPLSRSquared(subpx, py, m, npc, NULL, r2y, sdec);

    if(r2y->order > 0){
      cutoff = UPLSGetPCModelCutOff(r2y);

      r2_ = 0.f;
      cc = 0;
      for(i = 0; i < r2y->m[cutoff]->row; i++){
        for(j = 0; j < r2y->m[cutoff]->col; j++){
          r2_ += getTensorValue(r2y, cutoff, i, j);
          cc++;
        }
      }
      r2_ /=cc;
    }
    else{
      r2_ = 0.f;
    }

    (*r2) = (*q2) = r2_;

    DelTensor(&sdec);
    DelTensor(&r2y);
    DelTensor(&subpx);
    DelUPLSModel(&m);
  }
  else if(validation_type == 1){ /*LEAVE ONE OUT*/
    initDVector(&r2x);
    initTensor(&q2y);
    initTensor(&sdep);

    UPLSLOOCV(subax, ay, xautoscaling, yautoscaling, npc,
                      &r2x,
                      &q2y,
                      &sdep,
                      NULL, NULL, s);

    NewUPLSModel(&m);
    UPLS(subax, ay, npc, xautoscaling, yautoscaling, m, s);

    UPLSRSquared(subax, ay,  m, npc, m->r2x_model, m->r2y_model, m->sdec);

    if(q2y->order > 0){
      cutoff = UPLSGetPCModelCutOff(q2y);

      q2_ = 0.f;
      cc = 0;
      for(i = 0; i < q2y->m[cutoff]->row; i++){
        for(j = 0; j < q2y->m[cutoff]->col; j++){
          q2_ += getTensorValue(q2y, cutoff, i, j);
          cc++;
        }
      }
      q2_ /=cc;

      r2_ = 0.f;
      cc = 0;
      for(i = 0; i < m->r2y_model->m[cutoff]->row; i++){
        for(j = 0; j < m->r2y_model->m[cutoff]->col; j++){
          r2_ += getTensorValue(m->r2y_model, cutoff, i, j);
          cc++;
        }
      }
      r2_ /=cc;
    }
    else{
      r2_ = q2_ = 0.f;
    }

    (*q2) = q2_;
    (*r2) = r2_;

    DelUPLSModel(&m);
    DelDVector(&r2x);
    DelTensor(&q2y);
    DelTensor(&sdep);
  }
  else{ /*CROSS VALIDATION*/
    initDVector(&r2x);
    initTensor(&q2y);
    initTensor(&sdep);
    UPLSRandomGroupsCV(subax, ay, xautoscaling, yautoscaling, npc, ngroup, niter,
                  &r2x,
                  &q2y,
                  &sdep,
                  NULL, NULL, s);

    NewUPLSModel(&m);
    UPLS(subax, ay, npc, xautoscaling, yautoscaling, m, s);
    UPLSRSquared(subax, ay,  m, npc, m->r2x_model, m->r2y_model, m->sdec);
    if(q2y->order > 0){
      cutoff = UPLSGetPCModelCutOff(q2y);

      q2_ = 0.f;
      cc = 0;
      for(i = 0; i < q2y->m[cutoff]->row; i++){
        for(j = 0; j < q2y->m[cutoff]->col; j++){
          q2_ += getTensorValue(q2y, cutoff, i, j);
          cc++;
        }
      }
      q2_ /=cc;

      r2_ = 0.f;
      cc = 0;
      for(i = 0; i < m->r2y_model->m[cutoff]->row; i++){
        for(j = 0; j < m->r2y_model->m[cutoff]->col; j++){
          r2_ += getTensorValue(m->r2y_model, cutoff, i, j);
          cc++;
        }
      }
      r2_ /=cc;
    }
    else{
      r2_ = q2_ = 0.f;
    }

    (*q2) = q2_;
    (*r2) = r2_;

    DelUPLSModel(&m);
    DelDVector(&r2x);
    DelTensor(&q2y);
    DelTensor(&sdep);
  }
  DelTensor(&subax);
}

void PrintUPLSModel(UPLSMODEL* m)
{
  puts("t Scores");
  PrintMatrix(m->xscores);

  puts("u Scores");
  PrintMatrix(m->yscores);

  puts("---------------");
  puts("X Loadings");
  PrintTensor(m->xloadings);

  puts("Y Loadings");
  PrintTensor(m->yloadings);


  puts("X Weights");
  PrintTensor(m->xweights);

  puts("---------------");

  puts("Regression Coefficients");
  PrintDVector(m->b);

  puts("X Variance Explained");
  PrintDVector(m->xvarexp);

  puts("Recalculated Y");
  PrintTensor(m->recalculated_y);

  puts("Recalculated Residuals");
  PrintTensor(m->recalc_residuals);
}

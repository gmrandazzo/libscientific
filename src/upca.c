/* Implements UPCA algorithms.
 * Copyright (C) 2016-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "list.h"
#include "memwrapper.h"
#include "upca.h"
#include "tensor.h"
#include "matrix.h"
#include "vector.h"
#include "pca.h"
#include "preprocessing.h"
#include "numeric.h"


void NewUPCAModel(UPCAMODEL** m)
{
  (*m) = xmalloc(sizeof(UPCAMODEL));
  initMatrix(&(*m)->scores);
  initTensor(&(*m)->loadings);
  initDVector(&(*m)->varexp);
  initDVectorList(&(*m)->colaverage);
  initDVectorList(&(*m)->colscaling);
}

void DelUPCAModel(UPCAMODEL** m)
{
  DelMatrix(&(*m)->scores);
  DelTensor(&(*m)->loadings);
  DelDVector(&(*m)->varexp);
  DelDVectorList(&(*m)->colaverage);
  DelDVectorList(&(*m)->colscaling);
  xfree((*m));
}

int CheckTensor(tensor *X_){
  size_t i, j;
  for(i = 0; i < X_->order; i++){
    for(j = i+1; j < X_->order; j++){
      if(X_->m[i]->row == X_->m[j]->row
        && X_->m[i]->col == X_->m[j]->col)
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
 * X_ = input Tensor
 * X = tensor MeanCentered and Autoscaled
 * P = Matrix
 *
 * step 1 : Preprocessing X_ data: X_ is mean centered and after autoscaled.
 *
 * step 2: Initialize the component index a = 0, then for each compoent do:
 *
 * step 3: a = a +1
 *
 * step 4: t-start = column in X with largest Variance
 *
 *
 * step 5: P = t'X    ||P|| = 1 -> Loadings
 *
 * step 5.1:  additional NIPALS loop that decomposes P into u_1r_1'+ u_2r_2'. This loop should the be run through twice.
 *     a) h = 0; u-start is the "largest" column in P;
 *     b) r' = u'P/u'u
 *     c) norm to ||r|| = 1
 *     d) u = Pr
 *     e) h = h + 1; if h < 2 then back to "b"
 *     f) set P = ur' and procede with "step 6"
 *
 * step 6: t = X..P/||P||^2
 *
 * step 7: check convergente:
 *     if d = (t_new - t_old)'(t_new-t_old)/N*t_new't_new) > 10^-10
 *        back to "step 5"
 *
 * step 8: After convergence, residuals: E = X - t <scalar_product>
 * Then set X = E and extimante the next component by returning to step 3
 *
 */

void UPCA(tensor *X_, size_t npc, size_t autoscaling, UPCAMODEL *m, ssignal *s)
{
  if(CheckTensor(X_) == 0){
    size_t pc, i, j, k, order, column;
    double a, ss;
    tensor *X;
    dvector *colvar;
    dvector *t_old; /* row vector with size X->m->row */
    dvector *t_new;
    dvector *eval; /* t't */

    matrix *P;
    dvector *r;
    dvector *u;

    /* step 1 center and autoscale matrix */

    NewTensor(&X, X_->order);
    for(k = 0; k < X_->order; k++){
      NewTensorMatrix(X, k, X_->m[k]->row, X_->m[k]->col);
    }
    TensorPreprocess(X_, autoscaling, m->colaverage, m->colscaling, X);

   /* END Disable */

    if(npc > X->m[0]->col)
      npc = X->m[0]->col;

    ss = 0.f;
    for(k = 0; k < X->order; k++)
      for(i = 0; i < X->m[k]->row; i++)
        for(j = 0; j < X->m[k]->col; j++)
          ss += square(getTensorValue(X, k, i, j));

    #ifdef DEBUG
    puts("Processed Tensor");
    PrintTensor(X);
    #endif

    NewDVector(&t_old, X->m[0]->row); /* the object number is the same for all the matrix in tensor */
    NewDVector(&t_new, X->m[0]->row); /* the object number is the same for all the matrix in tensor */
    NewMatrix(&P, X->m[0]->col, X->order); /*Check that all the matrix of tensor have the same row and column */
    NewDVector(&eval, npc);

    /* step 2 and 3 */
    for(pc = 0; pc < npc; pc++){
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* step 4 select the column vector t with the largest column variance */
        order = column = 0;
        for(i = 0; i < X->order; i++){
          initDVector(&colvar);
          MatrixColVar(X->m[i], colvar);

          #ifdef DEBUG
          puts("Column Variance for chosing the best t-start column vector");
          PrintDVector(colvar);
          #endif

          /* check the column in tensor with the largest column variance */
          k = 0;
          for(j = 1; j < X->m[i]->col; j++){
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


        /* copy the best column vector into t */
        for(i = 0; i < X->m[order]->row; i++){
          setDVectorValue(t_old, i, getTensorValue(X, order, i, column));
        }

        #ifdef DEBUG
        puts("The best t-start vector");
        PrintDVector(t_old);
        #endif

        while(1){
          /*step 5*/
          MatrixSet(P, 0.f);
          DvectorTensorDotProduct(X, t_old, P);

          MatrixNorm(P, P);

          #ifdef DEBUG
          puts("P Loadings Matrix");
          PrintMatrix(P);
          #endif

          /*step 5.1*/
          #ifdef DEBUG
          puts("Nipals SHUNT");
          #endif
          NewDVector(&r, P->col);
          NewDVector(&u, P->row);

          /*u-start is the larget column in P*/
          initDVector(&colvar);
          MatrixColAverage(P, colvar);
          j = 0;
          for(i = 1; i < colvar->size; i++){
            if(getDVectorValue(colvar, i) > getDVectorValue(colvar, j))
              j = i;
            else
              continue;
          }
          DelDVector(&colvar);

          /*copy the best column vector from P to u*/
          for(i = 0; i < P->row; i++)
            setDVectorValue(u, i, getMatrixValue(P, i, j));

          #ifdef DEBUG
          puts("U-start");
          PrintDVector(u);
          #endif

          for(i = 0; i < 2; i++){ /* for the simplest case */

            /* 2.  p' := (t'X)/(t't)        Project the matrix E onto t in order to find the corresponding loading p
              * 3.  p' := p'/|p'|            Normalize the loading vector p to length 1
              * 4.  t_old := t
              *     t := (Xp)/(p'p)          Store the score vector t into uold and project the matrix E onto p in order to find corresponding score vector t
              */


            /* r' = u'P/u'u */
            DVectorSet(r, 0.f);
            DVectorMatrixDotProduct(P, u, r);
            a = DVectorDVectorDotProd(u, u);

            for(j = 0; j < r->size; j++)
              setDVectorValue(r, j, getDVectorValue(r, j) / a);

            /*||r|| = 1; */
            DVectNorm(r, r);

            #ifdef DEBUG
            puts("r vector");
            PrintDVector(r);
            #endif

            /*u = Pr */
            DVectorSet(u, 0.f); /* reset the u vector*/
            MatrixDVectorDotProduct(P, r, u);

            /*
            a = DVectorDVectorDotProd(r, r);

            for(j = 0; j < u->size; j++)
              setDVectorValue(u, j, getDVectorValue(u, j) / a);

            */
            #ifdef DEBUG
            puts("u vector");
            PrintDVector(u);
            #endif
          }

          /*step 5= f  set P = ur' */
          MatrixSet(P, 0.f);
          RowColOuterProduct(u, r, P);

          #ifdef DEBUG
          puts("Recomputed P Loadings Matrix from ur'");
          PrintMatrix(P);
          #endif


          DelDVector(&r);
          DelDVector(&u);

          /* step 7
          #ifdef DEBUG
          puts("old t score vector");
          PrintDVector(t_old);
          #endif
          */

          DVectorSet(t_new, 0.f);
          TensorMatrixDotProduct(X, P, t_new);

          a = Matrixnorm(P);
          a = square(a);

          for(i = 0; i < t_new->size; i++)
            setDVectorValue(t_new, i, getDVectorValue(t_new, i)/a);


          #ifdef DEBUG
          puts("new t score vector");
          PrintDVector(t_new);
          #endif

          if(calcConvergence(t_new, t_old) < UPCACONVERGENCE){
            /* storing results */
            MatrixAppendCol(m->scores, t_new);
            TensorAppendMatrix(m->loadings, P);
            setDVectorValue(eval, pc, DVectorDVectorDotProd(t_new, t_new));
            break;
          }
          else if(_isnan_(DVectorDVectorDotProd(t_new, t_new))){
            fprintf(stderr, "UPCA Error! The Solver Engine was Unable to Converge! Please Check your data.\n");
            fflush(stderr);
            /*   abort(); */
          }
          else{
            DVectorCopy(t_new, t_old);
          }
        }

        /* step 8 remove computed component */
        for(k = 0; k < X->order; k++){
          for(i = 0; i < X->m[k]->row; i++){
            for(j = 0; j < X->m[k]->col; j++){
              setTensorValue(X, k, i, j, getTensorValue(X, k, i, j) - (getDVectorValue(t_new, i)*getMatrixValue(P, j, k)));
            }
          }
        }

        /*Reset all*/
        MatrixSet(P, 0.f);
        DVectorSet(t_new, 0.f);
        DVectorSet(t_old, 0.f);
      }
    }

    calcVarExpressed(ss, eval, m->varexp);

    DelDVector(&eval);
    DelMatrix(&P);
    DelDVector(&t_old);
    DelDVector(&t_new);
    DelTensor(&X);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute Multi Way PCA. The (row,col) size of tensor submatrix are different!\n");
    /* abort(); */
  }
}

void UPCAScorePredictor(tensor *X_,
                        UPCAMODEL *model,
                        size_t npc,
                        matrix *pscores)
{
  size_t pc, i, j , k;
  tensor *X;
  dvector *t;

  double a;

  NewTensor(&X, X_->order);
  for(k = 0; k < X_->order; k++){
    NewTensorMatrix(X, k, X_->m[k]->row, X_->m[k]->col);
  }

  NewDVector(&t, X->m[0]->row);

  for(k = 0; k < X_->order; k++){
    for(j = 0; j < X_->m[k]->col; j++){
      /* Mean Centering */
      if(model->colaverage->size > 0){
        for(i = 0; i < X_->m[k]->row; i++){
          X->m[k]->data[i][j] = X->m[k]->data[i][j] - model->colaverage->d[k]->data[j];
        }
      }
      else{
        for(i = 0; i < X_->m[k]->row; i++){
          setTensorValue(X, k, i, j, getTensorValue(X_, k, i, j));
        }
      }

      if(model->colscaling->size > 0){
        /* Scaling to Column SDEV */
        for(i = 0; i < X->m[k]->row; i++){
          X->m[k]->data[i][j] = X->m[k]->data[i][j]/model->colscaling->d[k]->data[j];
        }
      }
    }
  }

  if(model->loadings->order < npc)
    npc = model->loadings->order;

  for(pc = 0; pc < npc; pc++){
    DVectorSet(t, 0.f);
    TensorMatrixDotProduct(X, model->loadings->m[pc], t);



    a = Matrixnorm(model->loadings->m[pc]);
    a = square(a);

    for(i = 0; i < t->size; i++)
      setDVectorValue(t, i, getDVectorValue(t, i)/a);

    MatrixAppendCol(pscores, t);


    for(k = 0; k < X->order; k++){
      for(i = 0; i < X->m[k]->row; i++){
        for(j = 0; j < X->m[k]->col; j++){
          setTensorValue(X, k, i, j, getTensorValue(X, k, i, j) - (getDVectorValue(t, i)*getMatrixValue(model->loadings->m[pc], j, k)));
        }
      }
    }
  }

  DelDVector(&t);
  DelTensor(&X);
}


/*
 * X = TP' + ....
 *
 */
void UPCAIndVarPredictor(matrix *T,
                         tensor *P,
                         dvectorlist *colaverage,
                         dvectorlist *colscaling,
                         size_t npc,
                         tensor *X)
{
  size_t pc, i, j, k;

  /*Allocate the output tensor*/
  for(i = 0; i < P->m[0]->col; i++)
    AddTensorMatrix(X, T->row, P->m[0]->row);
    /*AddTensorMatrix(X, T->row, P->m[0]->col);*/

  for(pc = 0; pc < npc; pc++){
    for(k = 0; k < X->order; k++){
      for(i = 0; i < X->m[k]->row; i++){
        for(j = 0; j < X->m[k]->col; j++){
          setTensorValue(X, k, i, j, getTensorValue(X, k, i, j) + (getMatrixValue(T, i, pc) * getTensorValue(P, pc, j, k)));
        }
      }
    }
  }

  /*WARNING da testare se tornano i valori predetti*/

 /*
  * xcolaverage->row and xcolscaling->row are the number of variable for each submatrix in X_, so these are equal to X_->m[i]->col
  *
  * xcolaverage->col and xcolscaling->col are the number of order present in X_, so these are equal to X_->order
  */

  if(colaverage->size > 0){
    for(k = 0; k < X->order; k++){
      for(i = 0; i < X->m[k]->row; i++){
        for(j = 0; j < X->m[k]->col; j++){
          if(colscaling->size > 0){
            X->m[k]->data[i][j] = (X->m[k]->data[i][j] + colscaling->d[k]->data[j])  + colaverage->d[k]->data[j];
          }
          else{
            X->m[k]->data[i][j] = (X->m[k]->data[i][j] + colscaling->d[k]->data[j]);
          }
        }
      }
    }
  }

}


void PrintUPCAModel(UPCAMODEL* m)
{
  puts("Scores");
  PrintMatrix(m->scores);

  puts("Loadings");
  PrintTensor(m->loadings);

  puts("Variance Explained");
  PrintDVector(m->varexp);
}

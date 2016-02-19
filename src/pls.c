/* pls.c
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
#include "pls.h"
#include "pca.h" /*Using: MatrixAutoScaling(); and calcVarExpressed(); */
#include "numeric.h" /* Using:  if(FLOAT_EQ(NumOne, NumTwo));*/
#include "metricspace.h"
#include <math.h>
#include <pthread.h>

void NewPLSModel(PLSMODEL** m)
{
  (*m) = xmalloc(sizeof(PLSMODEL));
  initMatrix(&(*m)->xscores);
  initMatrix(&(*m)->xloadings);
  initMatrix(&(*m)->xweights);
  initMatrix(&(*m)->yscores);
  initMatrix(&(*m)->yloadings);
  initDVector(&(*m)->b);
  initDVector(&(*m)->xvarexp);
  /*initDVector(&(*m)->yvarexp);*/
  initDVector(&(*m)->xcolaverage);
  initDVector(&(*m)->xcolscaling);
  initDVector(&(*m)->ycolaverage);
  initDVector(&(*m)->ycolscaling);
  initMatrix(&(*m)->r2y_model);
  initMatrix(&(*m)->r2y_validation);
  initMatrix(&(*m)->recalc_residuals);
  initMatrix(&(*m)->q2y);
  initMatrix(&(*m)->sdep);
  initMatrix(&(*m)->sdec);
  initMatrix(&(*m)->bias);
  initMatrix(&(*m)->recalculated_y);
  initMatrix(&(*m)->predicted_y);
  initMatrix(&(*m)->pred_residuals);
  initMatrix(&(*m)->r2q2scrambling);
  initMatrix(&(*m)->q2_sample_validation);
  initMatrix(&(*m)->sdep_sample_validation);
  initMatrix(&(*m)->q2_sample_validation_surface);
  initMatrix(&(*m)->sdep_sample_validation_surface);
}

void DelPLSModel(PLSMODEL** m)
{
  DelMatrix(&(*m)->xscores);
  DelMatrix(&(*m)->xloadings);
  DelMatrix(&(*m)->xweights);
  DelMatrix(&(*m)->yscores);
  DelMatrix(&(*m)->yloadings);
  DelDVector(&(*m)->b);
  DelDVector(&(*m)->xvarexp);
  /*DelDVector(&(*m)->yvarexp); */
  DelDVector(&(*m)->xcolaverage);
  DelDVector(&(*m)->xcolscaling);
  DelDVector(&(*m)->ycolaverage);
  DelDVector(&(*m)->ycolscaling);
  DelMatrix(&(*m)->r2y_model);
  DelMatrix(&(*m)->r2y_validation);
  DelMatrix(&(*m)->recalc_residuals);
  DelMatrix(&(*m)->q2y);
  DelMatrix(&(*m)->sdep);
  DelMatrix(&(*m)->sdec);
  DelMatrix(&(*m)->bias);
  DelMatrix(&(*m)->recalculated_y);
  DelMatrix(&(*m)->predicted_y);
  DelMatrix(&(*m)->pred_residuals);
  DelMatrix(&(*m)->r2q2scrambling);
  DelMatrix(&(*m)->q2_sample_validation);
  DelMatrix(&(*m)->sdep_sample_validation);
  DelMatrix(&(*m)->q2_sample_validation_surface);
  DelMatrix(&(*m)->sdep_sample_validation_surface);
  xfree((*m));
}

  /*
  *  Algorithm Geladi
  *
  * 1) take u  = some y[j] from Y
  *
  * 2) w' = u'X/u'u
  *
  * 3) w' = w'/||w||
  *
  * 4) t = Xw/w'w
  *
  * 5) q' = t'Y/t't
  *
  * 6) q' = q'/||q||
  *
  * 7) u = Yq/q'q
  *
  * 8) Compare the t in step 4 with th eone from the preceding iteration.
  * If they are equal (within a certain rounding error) go to step 9, else go to step 2.
  *
  * 9) p' = t'X/t't
  *
  * 10) p_new = p_old'/||p_old'||
  *
  * 11) t_new = t_old*||p_old'||
  *
  * 12) w_new' = w_old'*||p_old'||
  *
  */

void PLS(matrix *mx, matrix *my, size_t nlv, size_t xautoscaling, size_t yautoscaling, PLSMODEL* model, ssignal *s)
{
  if(nlv > 0){
    #ifdef DEBUG
    size_t step;
    #endif

    size_t i, j, pc, loop; /* pc = principal component */

    matrix *X; /* data matrix of mean centred object  for mx */
    matrix *Y; /* data matrix of mean centred object  for my */

    dvector *Y_avg;

    dvector *t; /* X score vector*/
    dvector *t_old; /* X score vector*/
    dvector *p; /* X loading vector */
    dvector *w; /* X weights vector */

    dvector *u; /* Y score vector*/
    dvector *q; /* Y loadings vector*/

    dvector *xeval;
    dvector *yeval;

    double min, max, dot_t, deltat = 0.f, dot_u, dot_q, dot_w, mod_p_old, ssx;

    if(mx->row == my->row){
      if(nlv > mx->col) /* if the number of principal component selected is major of the permitted */
        nlv = mx->col;

      NewMatrix(&X, mx->row, mx->col);
      /* MeanCenteredMatrix(mx, X);
      * We need to do the centered matrix manually here because these parameters could be useful for prediction.
      */
      MatrixCheck(mx);
      MatrixColAverage(mx, &(model->xcolaverage));

      for(j = 0; j < mx->col; j++){
        for(i = 0; i < mx->row; i++){
          X->data[i][j] = mx->data[i][j] - model->xcolaverage->data[j];
        }
      }

     /* AUTOSCALING
      * For autoscaling see:
      * Centering scaling and trasfomrations: improving the biological information content of metabolomics dataset
      * A van den Berg
      * BMC Genomics 2006, 7:142  doi:101 186/147-214-7-142
      */
      if(xautoscaling > 0){
        if(xautoscaling == 1){ /* SDEV Autoscaling */
          MatrixColSDEV(mx, &(model->xcolscaling));
        }
        else if(xautoscaling == 2){ /* RMS Autoscaling */
          MatrixColRMS(mx, &(model->xcolscaling));
        }
        else if(xautoscaling == 3){ /* PARETO Autoscaling */
          MatrixColSDEV(mx, &(model->xcolscaling));
          for(i = 0; i < model->xcolscaling->size; i++){
            model->xcolscaling->data[i] = sqrt(model->xcolscaling->data[i]);
          }
        }
        else if(xautoscaling == 4){ /* Range Scaling */
          for(i = 0; i < mx->col; i++){
            MatrixColumnMinMax(mx, i, &min, &max);
            DVectorAppend(&model->xcolscaling, (max - min));
          }
        }
        else if(xautoscaling == 5){ /* Level Scaling  */
          DVectorCopy(model->xcolaverage, &model->xcolscaling);
        }
        else{
          for(int i = 0; i < model->xcolaverage->size; i++){
            DVectorAppend(&model->xcolscaling, 1.0);
          }
        }


        for(j = 0; j < X->col; j++){
          if(model->xcolscaling->data[j] == 0){
            for(i = 0; i< X->row; i++){
              X->data[i][j] = 0.f;
            }
          }
          else{
            for(i = 0; i < X->row; i++){
              X->data[i][j] /= model->xcolscaling->data[j];
            }
          }
        }
      }

      NewMatrix(&Y, my->row, my->col);

      /*mean centering the y matrix*/
      MatrixColAverage(my, &(model->ycolaverage));
      for(j = 0; j < my->col; j++){
        for(i = 0; i < my->row; i++){
          Y->data[i][j] = my->data[i][j] - model->ycolaverage->data[j];
        }
      }


      if(yautoscaling > 0){
        if(yautoscaling == 1){ /* RMS Scaling: Divite for the Standar Deviation Column */
          MatrixColSDEV(my, &(model->ycolscaling));
        }
        else if(yautoscaling == 2){ /* RMS Scaling: Divite for the Root Mean Square of the column */
          MatrixColRMS(my, &(model->ycolscaling));
        }
        else if(yautoscaling == 3){ /* PARETO Autoscaling */
          MatrixColSDEV(my, &(model->ycolscaling));
          for(i = 0; i < model->ycolscaling->size; i++){
            model->ycolscaling->data[i] = sqrt(model->ycolscaling->data[i]);
          }
        }
        else if(yautoscaling == 4){ /* Range Scaling */
          for(i = 0; i < my->col; i++){
            MatrixColumnMinMax(my, i, &min, &max);
            DVectorAppend(&model->ycolscaling, (max - min));
          }
        }
        else if(yautoscaling == 5){ /* Level Scaling  */
          DVectorCopy(model->ycolaverage, &model->ycolscaling);
        }
        else{
          for(int i = 0; i < model->ycolaverage->size; i++){
            DVectorAppend(&model->ycolscaling, 1.0);
          }
        }

        for(j = 0; j < Y->col; j++){
          if(model->ycolscaling->data[j] == 0){
            for(i = 0; i< Y->row; i++){
              Y->data[i][j] = 0.f;
            }
          }
          else{
            for(i = 0; i < Y->row; i++){
              Y->data[i][j] /= model->ycolscaling->data[j];
            }
          }
        }
      }

      /* calculating the sum of squares for x and y in order to extimate
      * the % of variance explained
      */
      ssx = 0.f;
      for(i = 0; i < X->row; i++){
        for(j = 0; j < X->col; j++){
          ssx += square(X->data[i][j]);
        }
      }

      /*
      EXPLAIN THE VARIANCE FOR Y
      ssy = 0.f;
      for(i = 0; i < Y->row; i++)
        for(j = 0; j < Y->col; j++)
          ssy += square(getMatrixValue(Y, i, j));

      */

      NewDVector(&t, X->row);
      NewDVector(&t_old, X->row);
      NewDVector(&p, X->col);
      NewDVector(&w, X->col);

      NewDVector(&u, Y->row);
      NewDVector(&q, Y->col);

      NewDVector(&xeval, nlv);
      NewDVector(&yeval, nlv);


      /*RESIZE THE OUTPUT MATRIX*/
      ResizeMatrix(&(model->xscores), X->row, nlv);
      ResizeMatrix(&(model->xloadings), X->col, nlv);
      ResizeMatrix(&(model->xweights), X->col, nlv);
      ResizeMatrix(&(model->yscores), Y->row, nlv);
      ResizeMatrix(&(model->yloadings), Y->col, nlv);

      #ifdef DEBUG
      puts("X mean centered");
      PrintMatrix(X);
      puts("Y mean centered");
      PrintMatrix(Y);
      #endif

      for(pc = 0; pc < nlv; pc++){
        if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
          break;
        }
        else{
          /* Step 1: select the column vector u with the largest column average from Y  take u = some y_j */
          if(Y->col > 1){
            initDVector(&Y_avg);
            MatrixColVar(Y, &Y_avg);
            j = 0;
            for(i = 1; i < Y_avg->size; i++){
              if(Y_avg->data[i] > Y_avg->data[j])
                j = i;
              else
                continue;
            }
            DelDVector(&Y_avg);
          }
          else{ /* only one y dependent variable */
            j = 0;
          }

          /* copy the vector to the score u for computing weights w vector  */
          for(i = 0; i < u->size; i++)
            u->data[i] = Y->data[i][j];
          /* End Step 1 */

          #ifdef DEBUG
          step = 0;
          #endif
          loop = 0;
          while(1){
            #ifdef DEBUG
            printf("######### Step %u\n", (unsigned int)step);
            step++;
            #endif
            /* Step 2: Compute a weight vector w = u'X/u'u   (NB.: w is size of X_t->col = X->row ) */
            DVectorSet(w, 0.f); /* Reset the w vector */
            DVectorMatrixDotProduct(X, u, w);

            dot_u = DVectorDVectorDotProd(u,u);

            for(i = 0; i < w->size; i++){
              w->data[i] /= dot_u;
            }
            /* End Step2 */

            /* Step 3: w = w/||w1||   Normalize */
            DVectNorm(w, w);

            /* End Step 3 */

            #ifdef DEBUG
            printf("\n w Weight Vector for X matrix \n");
            PrintDVector(w);
            #endif

            /* Step 4:   t = Xw/w'w Compute a t score vector. t i size of X->row */
            DVectorSet(t, 0.f); /* Reset the t vector */
            MatrixDVectorDotProduct(X, w, t);

            dot_w = DVectorDVectorDotProd(w, w);

            for(i = 0; i < t->size; i++){
              t->data[i] /= dot_w;
            }

            #ifdef DEBUG
            printf("\n t Score Vector for X matrix \n");
            PrintDVector(t);
            #endif
            /* End Step 4 */

            if(Y->col > 1){
              /* Step 5: Compute the loadings Y vector q' = t'Y/t't  and normalize q = q/||q||*/
              DVectorSet(q, 0.f); /* Reset the c vector */
              DVectorMatrixDotProduct(Y, t, q);

              dot_t = DVectorDVectorDotProd(t, t);

              for(i = 0; i < q->size; i++){
                q->data[i] /= dot_t;
              }

              DVectNorm(q, q);
              /* End Step 5 */

              /* Step 6. Update the Y Scores u  u = Yq/q'q */
              DVectorSet(u, 0.f);
              MatrixDVectorDotProduct(Y, q, u);
              dot_q = DVectorDVectorDotProd(q, q);

              for(i = 0; i < u->size; i++){
                u->data[i] /= dot_q;
              }

              /* End step 6 */

              #ifdef DEBUG
              printf("\n New Score vector u for Y\n");
              PrintDVector(u);
              #endif
            }
            else{
              /*If the Y block has only one variable, step 5-8 can be omitted by putting q = 1 and no more iteration is necessary.*/
              DVectorSet(q, 1);
              //break;
            }

            /* Step 8 if t_old == t new with a precision PLSCONVERGENCE then stop iteration else restart from step 2 */

            if(loop == 0){
              for(i = 0; i < t->size; i++)
                t_old->data[i] = t->data[i];
            }
            else{
              deltat = 0.f;
              for(i = 0; i < t->size; i++){
                deltat = (t->data[i] - t_old->data[i]) * (t->data[i] - t_old->data[i]);
                t_old->data[i] = t->data[i];
              }
              deltat = sqrt(deltat);
              if(deltat < PLSCONVERGENCE){
                break;
              }
            }
            loop++;
            /* End step 8 */
          }

          /* Step 9 compute the loading vector for X: p' = t'X/t't  and y */
          DVectorMatrixDotProduct(X, t, p);

          dot_t = DVectorDVectorDotProd(t, t);
          for(i = 0; i < p->size; i++){
            p->data[i] /= dot_t;
          }

          #ifdef DEBUG
          printf("X Loadings\n");
          PrintDVector(p);
          #endif
          /* End Step 9*/

          mod_p_old = DvectorModule(p);

          /* Step 10  p_new = p_old'/||p_old'|| */
          DVectNorm(p, p);

          /*Step 11 t = t||p'||*/
          for(i = 0; i < t->size; i++){
            t->data[i] *= mod_p_old;
          }

          /*Step 12 */
          for(i = 0; i < w->size; i++){
            w->data[i] *= mod_p_old;
          }

          /* Step 13 Find the Regression Coefficient b:  b = u't/t't */
          dot_t = DVectorDVectorDotProd(t,t);

          DVectorAppend(&(model->b), (DVectorDVectorDotProd(u,t))/dot_t);

          /*End Step 13*/

          /* Step 14  Adjust X for what has been found: Xnew=X-tp'  and Y fort what has been found Ynew = Y - btp' */
          /* X */
          for(i = 0; i < X->row; i++){
            for(j = 0; j < X->col; j++){
              X->data[i][j] -= t->data[i]*p->data[j];
            }
          }

          /* Y */
          for(i = 0; i < Y->row; i++){
            for(j = 0; j < Y->col; j++){
              Y->data[i][j] -= (model->b->data[pc] * t->data[i] * q->data[j]);
            }
          }

          #ifdef DEBUG
          printf("\n X new \n");
          PrintMatrix(X);
          printf("\n Y new \n");
          PrintMatrix(Y);
          #endif


          /* Calculating eigenvalue for X and Y in order to estimate the explained variance for each component
          * module of u vector and module of t vector are the eigenvalue
          */
          dot_u = DVectorDVectorDotProd(u,u); /* Eigen Value for the component of Y */
          /* mod_t eigenvalue was previously calculated */
          /* setDVectorValue(yeval, pc, mod_u); NO MAKE SENSE EXPLAIN THE VARIANCE FOR Y*/
          xeval->data[pc] = dot_t; /* Adding the t eigenvalue */


          /* If more pairs (t,u) are needed go to 1. with X=Xnew and Y=Ynew.
          * Storing scores and loadings
          */

          for(i = 0; i < t->size; i++){
            model->xscores->data[i][pc] = t->data[i];
            model->yscores->data[i][pc] = u->data[i];
          }

          for(i = 0; i < p->size; i++){
            model->xloadings->data[i][pc] = p->data[i];
            model->xweights->data[i][pc] = w->data[i];
          }

          for(i = 0; i < q->size; i++){
            model->yloadings->data[i][pc] = q->data[i];
          }

          /*
          MatrixAppendCol(&(model->xscores), t);
          MatrixAppendCol(&(model->xloadings), p);
          MatrixAppendCol(&(model->xweights), w);
          MatrixAppendCol(&(model->yscores), u);
          MatrixAppendCol(&(model->yloadings), q);
          */
          mod_p_old = dot_q = dot_t = dot_u = dot_w = 0.f;
        }
      }

      calcVarExpressed(ssx, xeval, &(model->xvarexp));
      /*calcVarExpressed(ssy, yeval, &(model->yvarexp)); */

      /* PLS Recaculated */
      for(i = 1; i <= nlv; i++){
        matrix *recalculated_y;
        initMatrix(&recalculated_y);
        PLSYPredictor(model->xscores, model, i, &recalculated_y);

        for(j = 0; j < recalculated_y->col; j++){
          dvector *v = getMatrixColumn(recalculated_y, j);
          MatrixAppendCol(&model->recalculated_y, v);
          DelDVector(&v);
        }
        DelMatrix(&recalculated_y);
      }

      /* compute residuals */
      ResizeMatrix(&model->recalc_residuals, my->row, my->col*nlv);
      for(i = 0; i < model->recalculated_y->row; i++){
        for(j = 0; j < model->recalculated_y->col; j++){
          model->recalc_residuals->data[i][j] = model->recalculated_y->data[i][j] - my->data[i][(size_t)floor(j/nlv)];
        }
      }

      DelMatrix(&X);
      DelMatrix(&Y);

      DelDVector(&t);
      DelDVector(&t_old);
      DelDVector(&p);
      DelDVector(&w);

      DelDVector(&u);
      DelDVector(&q);

      DelDVector(&xeval);
      DelDVector(&yeval);
    }
    else{
      fprintf(stderr, "Unable to run PLS Calculation.\n The number of experiments differ (%u != %u)\n", (unsigned int)mx->row, (unsigned int)my->row);
    }
  }
  else{
    fprintf(stderr, "Unable to run PLS Calculation.\n The number of principal component are <= 0\n");
  }
}

void PLSBetasCoeff(PLSMODEL *model, size_t nlv, dvector **betas)
{
  /* compute beta coefficient
   Wstar = W *(P'*W)^(-1);
   B = Wstar*b’;
  */
  matrix *W, *P_, *B_;
  NewMatrix(&W, model->xweights->row, nlv);
  NewMatrix(&P_, nlv, model->xweights->row);
  NewMatrix(&B_, nlv, 1);;

  for(size_t j = 0; j < nlv; j++){
    for(size_t i = 0; i < model->xweights->row; i++){
      W->data[i][j] = model->xweights->data[i][j];
      P_->data[j][i] = model->xloadings->data[i][j];
    }
    B_->data[j][0] = model->b->data[j];
  }

  matrix *PW;
  NewMatrix(&PW, nlv, nlv);
  MatrixDotProduct(P_, W, PW);
  DelMatrix(&P_);

  matrix *PWinv;
  initMatrix(&PWinv);
  MatrixInversion(PW, &PWinv);
  DelMatrix(&PW);

  matrix *WStar;
  NewMatrix(&WStar, W->row, nlv);
  MatrixDotProduct(W, PWinv, WStar);
  DelMatrix(&PWinv);

  matrix *betas_;
  NewMatrix(&betas_, WStar->row, 1);
  MatrixDotProduct(WStar, B_, betas_);

  DVectorResize(betas, model->xweights->row);
  for(size_t i = 0; i < betas_->row; i++){
    (*betas)->data[i] = betas_->data[i][0];
  }
  // PrintDVector((*betas));
  DelMatrix(&WStar);
  DelMatrix(&betas_);

  DelMatrix(&B_);
  DelMatrix(&W);
}
/*
 * For prediction we need p', q', w', b coefficient, the column average and the column sdev (for mean centering and scaling the matrix), from the PLS Calibration.
 *
 * For estimate the score the procedure is this:
 *
 * Mean centering the X matrix with the previously column average made from the PLS Model
 * Autoscaling if is needed the X Matrix by the previously column standard deviation made from the PLS Model
 *
 * t = Xw
 * X = X - tp'
 *
 * For estimate the dependent variable
 * y = Sum btq'
 */
void PLSScorePredictor(matrix *mx, PLSMODEL *model, size_t nlv, matrix **xscores)
{
  size_t pc, i ,j;
  matrix *X;
  dvector *t;
  dvector *w;

  NewMatrix(&X, mx->row, mx->col);
  if(model->xcolaverage->size > 0){
    for(i = 0; i < mx->row; i++){
      for(j = 0; j < mx->col; j++){
        X->data[i][j] = mx->data[i][j] - model->xcolaverage->data[j];
      }
    }
  }
  else{
    for(i = 0; i < mx->row; i++){
      for(j = 0; j < mx->col; j++){
        X->data[i][j] = mx->data[i][j];
      }
    }
  }

  if(model->xcolscaling != NULL && model->xcolscaling->data != NULL
    && model->xcolscaling->size > 0){
    for(j = 0; j < X->col; j++){
      if(model->xcolscaling->data[j] == 0){
        for(i = 0; i< X->row; i++){
          X->data[i][j] = 0.f;
        }
      }
      else{
        /* Autoscaling for the column j of data by its standar deviation */
        for(i = 0; i < X->row; i++){
          X->data[i][j] /= model->xcolscaling->data[j];
        }
      }
    }
  }

  #ifdef DEBUG
  puts("Preprocessed Matrix to predict");
  PrintMatrix(X);
  #endif

  if(nlv > model->xweights->col)
    nlv = model->xweights->col;

  NewDVector(&t, X->row);
  NewDVector(&w, model->xweights->row);

  ResizeMatrix(xscores, X->row, nlv);

  for(pc = 0; pc < nlv; pc++){
    for(i = 0; i < model->xweights->row; i++){
      w->data[i] = model->xweights->data[i][pc];
    }

    #ifdef DEBUG
    puts("Selected w");
    PrintDVector(w);
    #endif

    MatrixDVectorDotProduct(X, w, t);

    #ifdef DEBUG
    puts("Recalculated t");
    PrintDVector(t);
    puts("---------------------");
    #endif

    /*  X = X - tp' Remove the estimated PLS component from X. */
    for(i = 0; i < X->row; i++){
      for(j = 0; j < X->col; j++){
        X->data[i][j] -= (t->data[i]* model->xloadings->data[j][pc]);
      }
    }

    for(i = 0; i < t->size; i++){
      (*xscores)->data[i][pc] = t->data[i];
      t->data[i] = 0.f;
    }
  }

  DelDVector(&w);
  DelMatrix(&X);
  DelDVector(&t);
}


/* Estimate the y dependent variable
 *
 * y = matrix
 * b = coefficient
 * t = score vector
 * q = loadings vector
 *
 * Y = sum b*t*q'
 */
void PLSYPredictor(matrix *tscore, PLSMODEL *model, size_t nlv, matrix **y)
{
  size_t pc, i, j;
  double _b;


  if(nlv > tscore->col)
    nlv = tscore->col;

  /*Allocate the y results*/
  ResizeMatrix(y, tscore->row, model->yloadings->row);

  for(pc = 0; pc < nlv; pc++){
    _b = model->b->data[pc];

    for(i = 0; i < tscore->row; i++){
      for(j = 0; j < model->yloadings->row; j++){
        (*y)->data[i][j] += _b * tscore->data[i][pc] * model->yloadings->data[j][pc];
      }
    }
  }


  if(model->ycolaverage->size > 0){
    for(j = 0; j < model->ycolaverage->size; j++){
      if(model->ycolscaling->size > 0 && model->ycolscaling != NULL && model->ycolscaling->data != NULL){
        for(i = 0; i < (*y)->row; i++){
          (*y)->data[i][j] *= model->ycolscaling->data[j];
        }
      }
      for(i = 0; i < (*y)->row; i++){
        (*y)->data[i][j] += model->ycolaverage->data[j];
      }
    }
  }
}

void PLSRSquared(matrix *mx, matrix *my, PLSMODEL *model, size_t nlv, matrix** r2y, matrix **sdec)
{
  size_t pc, i, j;
  matrix *recalcy;
  matrix *recalcscores;
  dvector *r2y_;
  dvector *sdec_;
  double n, d;

  if(nlv > model->b->size)
    nlv = model->b->size;

  for(pc = 1; pc <= nlv; pc++){
    initMatrix(&recalcy);
    initMatrix(&recalcscores);
    initDVector(&r2y_);
    initDVector(&sdec_);

    PLSScorePredictor(mx, model, pc, &recalcscores);

    PLSYPredictor(recalcscores, model, pc, &recalcy);

    /* r2y is not cumulative as x and is calculated for each y
    *
    * R2Y = 1 -  Sum((Y_real - Y_Extimated)^2) / Sum((Y_real - Y_med>)^2)
    */
    for(j = 0; j < recalcy->col; j++){
      n = d = 0.f;
      for(i = 0; i < recalcy->row; i++){
        n += square(my->data[i][j] - recalcy->data[i][j]);
        d += square(my->data[i][j] - model->ycolaverage->data[j]);
      }
      DVectorAppend(&r2y_, 1 - (n/d));
      DVectorAppend(&sdec_, sqrt(n/my->row));
    }

    MatrixAppendRow(sdec, sdec_);
    MatrixAppendRow(r2y, r2y_);

    DelDVector(&sdec_);
    DelDVector(&r2y_);
    DelMatrix(&recalcy);
    DelMatrix(&recalcscores);
  }
}

/* Used During Cross Validation LOO E Random Group
 * yss_err and xss_err are the PRESS predicted residual sums of squares
 */
void PLSRSquared_SSErr_SSTot(matrix *mx, matrix *my, PLSMODEL *model, size_t nlv, dvector** xss_err, dvector **xss_tot, matrix** yss_err, matrix** yss_tot, matrix** pred_y)
{
  size_t pc, i, j;
  matrix *recalcy;
  matrix *recalcx;
  matrix *recalcscores;
  dvector *tmp;
  double n, d;

  if(nlv > model->b->size){
    nlv = model->b->size;
  }

  if(((*yss_err)->row == 0 ||
    (*yss_err)->row != nlv ||
    (*yss_err)->col != my->col ||
    (*yss_tot)->row == 0 ||
    (*yss_tot)->row != nlv ||
    (*yss_tot)->col != mx->col)){
    ResizeMatrix(yss_err, nlv, my->col);
    ResizeMatrix(yss_tot, nlv, my->col);
  }
  else{
    for(i = 0; i < (*yss_err)->row; i++){
      for(j = 0; j < (*yss_err)->col; j++){
        (*yss_err)->data[i][j] = (*yss_tot)->data[i][j] = 0.f;
      }
    }
  }


  if(xss_err != NULL && xss_tot != NULL){
    if(((*xss_err)->size != my->row || (*xss_tot)->size != my->row)){
      DVectorResize(xss_err, nlv);
      DVectorResize(xss_tot, nlv);
    }
    else{
      for(i = 0; i < (*xss_err)->size; i++){
        (*xss_err)->data[i] = (*xss_tot)->data[i] = 0.f;
      }
    }
  }

  for(pc = 1; pc <= nlv; pc++){
    initMatrix(&recalcy);
    initMatrix(&recalcscores);

    PLSScorePredictor(mx, model, nlv, &recalcscores);

    PLSYPredictor(recalcscores, model, pc, &recalcy);

    /* r2y is not cumulative as x and is calculated for each y
    *
    * R2Y = 1 -  Sum((Y_real - Y_Extimated)^2) / Sum((Y_real - Y_med>)^2)
    */
    for(j = 0; j < recalcy->col; j++){
      n = d = 0.f;
      for(i = 0; i < recalcy->row; i++){
        n += square(my->data[i][j] - recalcy->data[i][j]);
        d += square(my->data[i][j] - model->ycolaverage->data[j]);
      }

      (*yss_err)->data[pc-1][j] = n;
      (*yss_tot)->data[pc-1][j] = d;

      if(pred_y != NULL){ /*for each component we have j y values*/
        tmp = getMatrixColumn(recalcy, j);
        MatrixAppendCol(pred_y, tmp);
        DelDVector(&tmp);
      }

    }

    /* r2x is cumulative because the X is all the matrix */
    if(xss_err != NULL && xss_tot != NULL){
      initMatrix(&recalcx);
      PCAIndVarPredictor(recalcscores, model->xloadings, model->xcolaverage, model->xcolscaling, pc, &recalcx);

      n = d = 0.f;
      for(j = 0; j < recalcx->col; j++){
        for(i = 0; i < recalcx->row; i++){
          n += square(recalcx->data[i][j] - mx->data[i][j]);
          d += square(mx->data[i][j] - model->xcolaverage->data[j]);
        }
      }

      (*xss_err)->data[pc-1] = n;
      (*xss_tot)->data[pc-1] = d;

      DelMatrix(&recalcx);
    }

    DelMatrix(&recalcy);
    DelMatrix(&recalcscores);
  }
}

void PLSYScrambling(matrix *mx, matrix *my,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t block,
                        size_t valtype, size_t rgcv_group, size_t rgcv_iterations,
                        matrix **r2q2scrambling, size_t nthreads, ssignal *s)
{
  size_t scrambiterations, iterations_, i, j, k, n, y_, blocksize;
  int id;
  double temp;
  matrix *randomY, *sorted_y_id, *sorty, *gid;
  matrix *tmpq2;
  PLSMODEL *tmpmod;

  dvector *yaverage;

  initDVector(&yaverage);
  MatrixColAverage(my, &yaverage);

  srand(mx->row*mx->col*my->col*block);
  NewMatrix(&randomY, my->row, my->col);
  NewMatrix(&sorted_y_id, my->row, my->col);

  NewMatrix(&sorty, my->row, 2);
  for(j = 0; j < my->col; j++){
    for(i = 0; i < my->row; i++){
      sorty->data[i][0] = my->data[i][j];
      sorty->data[i][1] = i;
    }
    MatrixSort(sorty, 0);

    for(i = 0; i < my->row; i++){
      sorted_y_id->data[i][j] = sorty->data[i][1];
    }
  }
  DelMatrix(&sorty);


  /*calcualte the block size for the rotate matrix*/
  blocksize = (size_t)ceil(mx->row/(double)block);
  blocksize += (size_t)ceil((float)((blocksize*block) - mx->row)/  block);

  NewMatrix(&gid, block, blocksize);
  MatrixSet(gid, -2);
  /* Crate the boxes to fill -2 means no value to fill, -1 means value to fill*/
  for(i = 0, j = 0, k = 0; i < mx->row; i++){
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

  /*Create a the r2q2scrambling matrix */
  ResizeMatrix(r2q2scrambling, 1+scrambiterations, my->col*3); /* each row is a model. first my->col columns are the r2 scrambled/real, second are the r2 scrabled/scrambled and third the q2 scrabled/scrabmled */

  /*First row is the model not scrambled...*/
  initMatrix(&tmpq2);

  if(valtype == 0){
    PLSLOOCV(mx, my, xautoscaling, yautoscaling, nlv, &tmpq2, NULL, NULL, NULL, NULL, nthreads, s);
  }
  else{
    PLSRandomGroupsCV(mx, my, xautoscaling, yautoscaling, nlv, rgcv_group, rgcv_iterations, &tmpq2, NULL, NULL, NULL, NULL, nthreads, s);
  }

  NewPLSModel(&tmpmod);
  PLS(mx, my, nlv, xautoscaling, yautoscaling, tmpmod, s);
  PLSRSquared(mx, my, tmpmod, nlv, &(tmpmod->r2y_model), &(tmpmod->sdec));

  /* Calculate y real vs yscrambled and add other r2 q2 */
  size_t r2cutoff = GetLVCCutoff(tmpmod->r2y_model);
  size_t q2cutoff = GetLVCCutoff(tmpq2);
  for(j = 0; j < my->col; j++){
    double rss = 0.f, tss = 0.f;
    for(i = 0; i < my->row; i++){
      rss += square(my->data[i][j] - my->data[i][j]);
      tss += square(my->data[i][j] - yaverage->data[j]);
    }
    (*r2q2scrambling)->data[0][j] = 1 - (rss/tss);
    (*r2q2scrambling)->data[0][j+my->col] = tmpmod->r2y_model->data[r2cutoff][j];
    (*r2q2scrambling)->data[0][j+my->col+my->col] = tmpq2->data[q2cutoff][j];
  }

  DelPLSModel(&tmpmod);
  DelMatrix(&tmpq2);

  for(y_ = 0; y_ < sorted_y_id->col; y_++){

    /* START WITH THE ORDERED Y_*/
    k = 0;
    for(i = 0; i < gid->row; i++){
      for(j = 0; j < gid->col; j++){
        if(gid->data[i][j] >= -1){
          gid->data[i][j] = sorted_y_id->data[k][y_];
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

        /*
        printf("Shuffled Y ID %d\n", (int)iterations_);
        PrintMatrix(gid);
        */

        /*Fill the shifted y*/
        n = 0;
        for(i = 0; i < gid->row; i++){
          for(j = 0; j < gid->col; j++){
            id = gid->data[i][j];
            if(id > -1){
              for(k = 0; k < my->col; k++){
                randomY->data[n][k] = my->data[id][k];
              }
              n++;
            }
            else{
              continue;
            }
          }
        }

        /*
        puts("MX");
        PrintMatrix(mx);
        puts("MY");
        PrintMatrix(my);
        puts("RANDOMY");
        PrintMatrix(randomY);
        */
        /* Calculate calculate Q2 for y predicted...*/

        initMatrix(&tmpq2);

        if(valtype == 0){
          PLSLOOCV(mx, randomY, xautoscaling, yautoscaling, nlv, &tmpq2, NULL, NULL, NULL, NULL, nthreads, s);
        }
        else{
          PLSRandomGroupsCV(mx, randomY, xautoscaling, yautoscaling, nlv, rgcv_group, rgcv_iterations, &tmpq2, NULL, NULL, NULL, NULL, nthreads, s);
        }

        NewPLSModel(&tmpmod);
        PLS(mx, randomY, nlv, xautoscaling, yautoscaling, tmpmod, s);
        PLSRSquared(mx, randomY, tmpmod, nlv, &(tmpmod->r2y_model), &(tmpmod->sdec));

        /* Calculate y real vs yscrambled and add other r2 q2 */
        size_t r2cutoff = GetLVCCutoff(tmpmod->r2y_model);
        size_t q2cutoff = GetLVCCutoff(tmpq2);
        for(j = 0; j < my->col; j++){
          double rss = 0.f, tss = 0.f;
          for(i = 0; i < my->row; i++){
            rss += square(my->data[i][j] - randomY->data[i][j]);
            tss += square(my->data[i][j] - yaverage->data[j]);
          }
          (*r2q2scrambling)->data[iterations_+1][j] = 1 - (rss/tss);
          (*r2q2scrambling)->data[iterations_+1][j+my->col] = tmpmod->r2y_model->data[r2cutoff][j];
          (*r2q2scrambling)->data[iterations_+1][j+my->col+my->col] = tmpq2->data[q2cutoff][j];
        }

        DelPLSModel(&tmpmod);
        DelMatrix(&tmpq2);
        iterations_++;
      }
    }
  }

  DelDVector(&yaverage);
  DelMatrix(&sorted_y_id);
  DelMatrix(&randomY);
  DelMatrix(&gid);
}

/* Cross Validation
 *
 * 1) Divide the dataset in "g" random group
 * 2) for each randomization run:
 *    for each group thake this out and run a PLS model with all the remains groups
 * 3) Predict the response value for the out group and compute PRESS (y_pred - y_r)^2 and the Sum of square of (y_r - y_mid)^2
 * 4) for each randomization divide all the sum by the number of group and sum this untill all the computation are done.
 */


typedef struct{

  matrix *mx, *my; /*INPUT*/

  matrix *predicted_y;  /*OUPUT*/
  uivector *predictioncounter; /*OUPUT*/

  size_t xautoscaling, yautoscaling, nlv, group; /*INPUT*/
  unsigned int srand_init;
} rgcv_th_arg;


void *RandomGroupCVModel(void *arg_)
{
  size_t i, j, k, n, g, lv;
  rgcv_th_arg *arg;
  matrix *gid; /* randomization and storing id for each random group into a matrix */

  /* Matrix for compute the PLS models for groups */
  matrix *subX;
  matrix *subY;
  PLSMODEL *subm;

  /*matrix to predict*/
  matrix *predictX;
  matrix *realY;
  matrix *predicty;


  arg = (rgcv_th_arg*) arg_;

  NewMatrix(&gid, arg->group, (size_t)ceil(arg->mx->row/(double)arg->group));

  /* Divide in group  all the Dataset */
  MatrixSet(gid, -1);

  /* step 1 generate the random groups */
  k = 0;
  for(i = 0; i <  gid->row; i++){
    for(j = 0; j <  gid->col; j++){
      do{
        /*n = randInt(0, arg->mx->row);*/
        n = (size_t)myrand_r(&arg->srand_init) % (arg->mx->row);
      } while(ValInMatrix(gid, n) == 1 && k < (arg->mx->row));
      if(k < arg->mx->row){
        gid->data[i][j] = n;
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
      /*
    printf("Excuded the group number %u\n", (unsigned int)g);
    puts("Sub Model\nX:");
    PrintArray(subX);
    puts("Y:");
    PrintArray(subY);

    puts("\n\nPredict Group\nX:");
    PrintArray(predictX);
    puts("RealY:");
    PrintArray(realY);
    */
  /*step 2*/
  for(g = 0; g < gid->row; g++){ /*For aeach group */
    /* Estimate how many objects are inside the sub model without the group "g" */
    n = 0;
    for(i = 0; i < gid->row; i++){
      if(i != g){
        for(j = 0; j < gid->col; j++){
          if((int)gid->data[i][j] != -1)
            n++;
          else
            continue;
        }
      }
      else
        continue;
    }

    /*Allocate the submodel*/
    NewMatrix(&subX, n, arg->mx->col);
    NewMatrix(&subY, n, arg->my->col);

    /* Estimate how many objects are inside the group "g" to predict*/
    n = 0;
    for(j = 0; j < gid->col; j++){
      if((int)gid->data[g][j] != -1)
        n++;
      else
        continue;
    }


    /*Allocate the */
    NewMatrix(&predictX, n, arg->mx->col);
    NewMatrix(&realY, n, arg->my->col);

    /* copy the submodel values */

    for(i = 0, k = 0; i < gid->row; i++){
      if(i != g){
        for(j = 0; j < gid->col; j++){
          size_t a =  (size_t)gid->data[i][j]; /* get the row index */
          if(a != -1){
            for(n = 0; n < arg->mx->col; n++){
              subX->data[k][n] = arg->mx->data[a][n];
            }
            for(n = 0; n < arg->my->col; n++){
              subY->data[k][n] = arg->my->data[a][n];
            }
            k++;
          }
          else{
            continue;
          }
        }
      }
      else{
        continue;
      }
    }

    /* copy the objects to predict into predictmx*/
    for(j = 0, k = 0; j < gid->col; j++){
      size_t a = (size_t)gid->data[g][j];
      if(a != -1){
        for(n = 0; n < arg->mx->col; n++){
          predictX->data[k][n] = arg->mx->data[a][n];
        }
        for(n = 0; n < arg->my->col; n++){
          realY->data[k][n] = arg->my->data[a][n];
        }
        k++;
      }
      else{
        continue;
      }
    }

    NewPLSModel(&subm);

    PLS(subX, subY, arg->nlv, arg->xautoscaling, arg->yautoscaling, subm, NULL);

    /* Predict Y for each latent variable */
    initMatrix(&predicty);

    for(lv = 1; lv <= arg->nlv; lv++){
      matrix *recalcy;
      matrix *recalcscores;

      initMatrix(&recalcy);
      initMatrix(&recalcscores);

      PLSScorePredictor(predictX, subm, lv, &recalcscores);
      PLSYPredictor(recalcscores, subm, lv, &recalcy);

      for(j = 0; j < recalcy->col; j++){
        dvector *tmp = getMatrixColumn(recalcy, j);
        MatrixAppendCol(&predicty, tmp);
        DelDVector(&tmp);
      }

      DelMatrix(&recalcy);
      DelMatrix(&recalcscores);
    }

    for(j = 0, k = 0; j < gid->col; j++){
      size_t a = (size_t)gid->data[g][j]; /*riga dell'oggetto....*/
      if(a != -1){
        arg->predictioncounter->data[a] += 1; /* this object was visited */
        /* updating y */
        for(n = 0; n < predicty->col; n++){
          arg->predicted_y->data[a][n] +=  predicty->data[k][n];
        }

        k++;
      }
      else{
        continue;
      }
    }

    DelMatrix(&predicty);
    DelPLSModel(&subm);
    DelMatrix(&subX);
    DelMatrix(&subY);
    DelMatrix(&predictX);
    DelMatrix(&realY);
  }
  DelMatrix(&gid);
  return 0;
}

void PLSRandomGroupsCV(matrix *mx, matrix *my, size_t xautoscaling, size_t yautoscaling, size_t nlv, size_t group, size_t iterations, matrix **q2y, matrix **sdep, matrix **bias, matrix **predicted_y, matrix **pred_residuals, size_t nthreads, ssignal *s)
{
  if(nlv > 0 && mx->row == my->row && group > 0 && iterations > 0){
    size_t th, iterations_, i, j;

    pthread_t *threads;

    dvector *ymean;
    uivector *predictcounter;
    matrix *sum_ypredictions;

    if(nlv > mx->col){
      nlv = mx->col;
    }

    NewMatrix(&sum_ypredictions, my->row, my->col*nlv); /* each component have my->col ypsilon */
    NewUIVector(&predictcounter, my->row);

    /* each thread have its argument type */
    threads = xmalloc(sizeof(pthread_t)*nthreads);
    rgcv_th_arg *arg;
    arg = xmalloc(sizeof(rgcv_th_arg)*nthreads);

    /*
    iterations_ = 0;
    while(iterations_ <  iterations){*/
    for(iterations_ = 0; iterations_ < iterations; iterations_ += nthreads){
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* Create independent threads each of which will execute function */
        for(th = 0; th < nthreads; th++){
          arg[th].mx = mx;
          arg[th].my = my;
          arg[th].group = group;
          arg[th].nlv = nlv;
          arg[th].xautoscaling = xautoscaling;
          arg[th].yautoscaling = yautoscaling;
          arg[th].srand_init = (unsigned int) group + mx->row + my->col + iterations + th + iterations_;
          NewMatrix(&arg[th].predicted_y, my->row, my->col*nlv); /* each component have my->col ypsilon */
          NewUIVector(&arg[th].predictioncounter, my->row);
          pthread_create(&threads[th], NULL, RandomGroupCVModel, (void*) &arg[th]);
        }

        /* Wait till threads are complete before main continues. Unless we  */
        /* wait we run the risk of executing an exit which will terminate   */
        /* the process and all threads before the threads have completed.   */
        for(th = 0; th < nthreads; th++){
          pthread_join(threads[th], NULL);
        }

        /* finalize thread outputs and free the memory.....*/
        for(th = 0; th < nthreads; th++){
          for(i = 0; i < arg[th].predicted_y->row; i++){
            for(j = 0; j < arg[th].predicted_y->col; j++){
              sum_ypredictions->data[i][j] += arg[th].predicted_y->data[i][j];
            }

            predictcounter->data[i] += arg[th].predictioncounter->data[i];
          }

          DelUIVector(&arg[th].predictioncounter);
          DelMatrix(&arg[th].predicted_y);
        }
      }
    }

    /*Finalize the output by dividing for the number of times that the object was predicted*/

    if(predicted_y != NULL){
      ResizeMatrix(predicted_y, sum_ypredictions->row, sum_ypredictions->col);
    }

    if(pred_residuals != NULL){
      ResizeMatrix(pred_residuals, my->row, my->col*nlv); /* each component have my->col ypsilon */
    }

    for(i = 0; i < sum_ypredictions->row; i++){
      for(j = 0; j < sum_ypredictions->col; j++){
        sum_ypredictions->data[i][j] /= (double)predictcounter->data[i];

        if(predicted_y != NULL)
          (*predicted_y)->data[i][j] = sum_ypredictions->data[i][j];

        if(pred_residuals != NULL)
          (*pred_residuals)->data[i][j] = sum_ypredictions->data[i][j] - my->data[i][(size_t)floor(j/nlv)];
      }
    }

    /*Calculate the Q2 and SDEP and Bias */
    if(q2y != NULL)
      ResizeMatrix(q2y, nlv, my->col);

    if(sdep != NULL)
      ResizeMatrix(sdep, nlv, my->col);

    if(bias != NULL)
      ResizeMatrix(bias, nlv, my->col);

    if(q2y != NULL || sdep != NULL){
      initDVector(&ymean);
      MatrixColAverage(my, &ymean);

      for(size_t lv = 0; lv < nlv; lv++){
        for(j = 0; j < my->col; j++){
          double ssreg = 0.f;
          double sstot = 0.f;

          /*y = m x + k: bias is the m angular coefficient */
          /*double ypredaverage = 0.f; used to calculate the k */
          for(i = 0; i < my->row; i++){
            ssreg += square(my->data[i][j] - sum_ypredictions->data[i][my->col*lv+j]);
            sstot += square(my->data[i][j] - ymean->data[j]);
            /*ypredaverage = sum_ypredictions->data[i][j+lv]; */
          }

          if(bias != NULL){
            double sum_xi = 0.f, sum_yi = 0.f;
            /*ypredaverage /= (double)my->row;*/
            for(i = 0; i < my->row; i++){
              sum_yi+=(sum_ypredictions->data[i][j+lv]*(my->data[i][j]-ymean->data[j]));
              sum_xi+=(my->data[i][j]*(my->data[i][j]-ymean->data[j]));
            }

            (*bias)->data[lv][j] = fabs(1 - sum_yi/sum_xi);
            /*k = Y-(X*b);*/
          }

          if(q2y != NULL)
            (*q2y)->data[lv][j] = 1.f - (ssreg/sstot);

          if(sdep != NULL)
            (*sdep)->data[lv][j] = sqrt(ssreg/(double)my->row);


        }
      }

      DelDVector(&ymean);
    }

    DelMatrix(&sum_ypredictions);
    DelUIVector(&predictcounter);
    xfree(threads);
    xfree(arg);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute PLS Random Group Cross Validation!!\n");
  }
}


/* Leave One Ouut Validation
 *
 * 1) remove one object
 * 2) calculate the model
 * 3) predict the removed object and so the r2 and q2
 */

typedef struct{
  matrix *submx, *submy, *predmx, *predmy, *pred_y;
  size_t nlv, xautoscaling, yautoscaling;
} loocv_th_arg;

void *PLSLOOModel(void *arg_)
{
  size_t lv, j;
  loocv_th_arg *arg;
  arg = (loocv_th_arg*) arg_;

  PLSMODEL *subm;
  NewPLSModel(&subm);

  PLS(arg->submx, arg->submy, arg->nlv, arg->xautoscaling, arg->yautoscaling, subm, NULL);

  /* Predict Y for each latent variable */
  for(lv = 1; lv <= arg->nlv; lv++){
    matrix *recalcy;
    matrix *recalcscores;

    initMatrix(&recalcy);
    initMatrix(&recalcscores);

    PLSScorePredictor(arg->predmx, subm, lv, &recalcscores);
    PLSYPredictor(recalcscores, subm, lv, &recalcy);

    for(j = 0; j < recalcy->col; j++){
      dvector *tmp = getMatrixColumn(recalcy, j);
      MatrixAppendCol(&arg->pred_y, tmp);
      DelDVector(&tmp);
    }

    DelMatrix(&recalcy);
    DelMatrix(&recalcscores);
  }

  DelPLSModel(&subm);
  return 0;
}

void PLSLOOCV(matrix* mx, matrix* my, size_t xautoscaling, size_t yautoscaling, size_t nlv, matrix** q2y, matrix** sdep, matrix **bias, matrix** predicted_y, matrix **pred_residuals, size_t nthreads, ssignal *s)
{
 if(nlv > 0 && mx->row == my->row){
    size_t i, j, k, l, th, model;
    pthread_t *threads;
    loocv_th_arg *arg;

    if(nlv > mx->col){
      nlv = mx->col;
    }

    dvector *ymean;
    matrix *loopredictedy;
    NewMatrix(&loopredictedy, my->row, my->col*nlv);

    threads = xmalloc(sizeof(pthread_t)*nthreads);
    arg = xmalloc(sizeof(loocv_th_arg)*nthreads);

    /* initialize threads arguments.. */
    for(th = 0; th < nthreads; th++){
      arg[th].nlv = nlv;
      arg[th].xautoscaling = xautoscaling;
      arg[th].yautoscaling = yautoscaling;
      NewMatrix(&arg[th].submx, mx->row-1, mx->col);
      NewMatrix(&arg[th].submy, my->row-1, my->col);
      NewMatrix(&arg[th].predmx, 1, mx->col);
      NewMatrix(&arg[th].predmy, 1, my->col);
    }


    for(model = 0; model < mx->row; model += nthreads){ /* we compute mx->row models  */
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{

        /* copy data into subX, subY, predictX and predictY for each thread argument
         * and run the thread
         */
        for(th = 0; th < nthreads; th++){
          if(th+model < mx->row){
            l = 0;
            for(j = 0; j < mx->row; j++){
              if(j != model+th){
                for(k = 0; k  < mx->col; k++){
                  /*setMatrixValue(arg[th].submx, l, k, getMatrixValue(mx, j, k));*/
                  arg[th].submx->data[l][k] = mx->data[j][k];
                }
                for(k = 0; k < my->col; k++){
                  /*setMatrixValue(arg[th].submy, l, k, getMatrixValue(my, j, k));*/
                  arg[th].submy->data[l][k] = my->data[j][k];
                }
                l++;
              }
              else{
                for(k = 0; k < mx->col; k++){
                  arg[th].predmx->data[0][k] = mx->data[j][k];
                }
                for(k = 0; k < my->col; k++){
                  arg[th].predmy->data[0][k] = my->data[j][k];
                }
              }
            }

            initMatrix(&arg[th].pred_y);
            pthread_create(&threads[th], NULL, PLSLOOModel, (void*) &arg[th]);
          }
          else{
            continue;
          }
        }

        /* Wait till threads are complete before main continues. Unless we  */
        /* wait we run the risk of executing an exit which will terminate   */
        /* the process and all threads before the threads have completed.   */
        for(th = 0; th < nthreads; th++){
          if(th+model < mx->row){
            pthread_join(threads[th], NULL);
          }
          else{
            continue;
          }
        }

        /*Collapse the threads output*/
        for(th = 0; th < nthreads; th++){
          if(th+model < mx->row){
            for(j = 0; j < arg[th].pred_y->col; j++){
              loopredictedy->data[model+th][j] = arg[th].pred_y->data[0][j];
            }
            DelMatrix(&arg[th].pred_y);
          }
        }
      }
    }

    /*Delete thread arguments*/

    for(th = 0; th < nthreads; th++){
      DelMatrix(&arg[th].submx);
      DelMatrix(&arg[th].submy);
      DelMatrix(&arg[th].predmx);
      DelMatrix(&arg[th].predmy);
    }

    /*Finalize the output by dividing for the number of models*/
    if(predicted_y != NULL){
      ResizeMatrix(predicted_y, loopredictedy->row, loopredictedy->col);
      MatrixCopy(loopredictedy, predicted_y);
    }

    if(pred_residuals != NULL){
      ResizeMatrix(pred_residuals, my->row, my->col*nlv); /* each component have my->col ypsilon */
      for(i = 0; i < loopredictedy->row; i++){
        for(j = 0; j < loopredictedy->col; j++){
          (*pred_residuals)->data[i][j] = loopredictedy->data[i][j] - my->data[i][(size_t)floor(j/nlv)];
        }
      }
    }

    /*Calculate the Q2 and SDEP */
    if(q2y != NULL)
      ResizeMatrix(q2y, nlv, my->col);

    if(sdep != NULL)
      ResizeMatrix(sdep, nlv, my->col);

    if(bias != NULL)
      ResizeMatrix(bias, nlv, my->col);

    if(sdep != NULL || q2y != NULL){
      initDVector(&ymean);
      MatrixColAverage(my, &ymean);

      for(size_t lv = 0; lv < nlv; lv++){
        for(j = 0; j < my->col; j++){
          double ssreg = 0.f;
          double sstot = 0.f;
          for(i = 0; i < my->row; i++){
            ssreg += square(loopredictedy->data[i][my->col*lv+j] - my->data[i][j]);
            sstot += square(my->data[i][j] - ymean->data[j]);
          }

          if(bias != NULL){
            double sum_yi = 0.f, sum_xi = 0.f;
            /*ypredaverage /= (double)my->row;*/
            for(i = 0; i < my->row; i++){
              sum_yi+=(loopredictedy->data[i][my->col*lv+j]*(my->data[i][j]-ymean->data[j]));
              sum_xi+=(my->data[i][j]*(my->data[i][j]-ymean->data[j]));
            }
            /*sum_yi/sum_xi = m */
            (*bias)->data[lv][j] = fabs(1 - sum_yi/sum_xi);
            /*k = Y-(X*b);*/
          }

          if(q2y != NULL)
            (*q2y)->data[lv][j] = 1.f - (ssreg/sstot);

          if(sdep != NULL)
            (*sdep)->data[lv][j] = sqrt(ssreg/(double)my->row);
        }

      }
      DelDVector(&ymean);
    }

    DelMatrix(&loopredictedy);
    xfree(arg);
    xfree(threads);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute PLS Leave One Out Validation!!\n");
  }
}


void PLSStaticSampleValidator(matrix *mx, matrix *my, uivector *obj_class,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t sample_size, size_t niters,
                        size_t rgcv_group, size_t rgcv_iterations, size_t nthreads,
                        matrix **q2_distr, matrix **sdep_distr, uivector **bestid, ssignal *s)
{
  if(obj_class->size == mx->row){
    size_t iter, i, j, nclass = 0, srand_init;
    matrix *xsample, *ysample, *tmpq2, *tmpsdep;
    uivector *tmpid;
    double bestq2 = -9999;
    for(i = 0; i < obj_class->size; i++){
      if(obj_class->data[i] > nclass){
        nclass = obj_class->data[i];
      }
      else{
        continue;
      }
    }

    initMatrix(&tmpq2);
    initMatrix(&tmpsdep);

    srand_init = mx->row+mx->col+my->col+sample_size+niters+nclass;

    NewMatrix(&xsample, sample_size, mx->col);
    NewMatrix(&ysample, sample_size, my->col);
    if(bestid != NULL)
      UIVectorResize(bestid, sample_size);
    NewUIVector(&tmpid, sample_size);
    for(iter = 0; iter < niters; iter++){
      if(nclass == 0){ /* simple bootstrap on all objects */
        for(i = 0; i < sample_size; i++){
          do{
            size_t n = (size_t)myrand_r((unsigned *)&srand_init) % (mx->row);
            tmpid->data[i] = n;
            if(n > obj_class->size){
              continue;
            }
            else{
              for(j = 0; j < mx->col; j++){
                xsample->data[i][j] = mx->data[n][j];
              }

              for(j = 0; j < my->col; j++){
                ysample->data[i][j] = my->data[n][j];
              }
              break;
            }
          }while(1);
        }
      }
      else{ /* select a random object from each class */
        size_t a_class = 0;
        for(i = 0; i < sample_size; i++){
          do{
            size_t n = (size_t)myrand_r((unsigned *)&srand_init) % (mx->row);
            tmpid->data[i] = n;
            if(n > obj_class->size){
              continue;
            }
            else{
              if(obj_class->data[n] == a_class){
                if(a_class == nclass){
                  a_class =  (size_t)myrand_r((unsigned *)&srand_init) % (nclass); /*all class are already completed then restart random...*/
                }
                else{
                  a_class++;
                }

                for(j = 0; j < mx->col; j++){
                  xsample->data[i][j] = mx->data[n][j];
                }

                for(j = 0; j < my->col; j++){
                  ysample->data[i][j] = my->data[n][j];
                }
                break;
              }
              else{
                continue;
              }
            }
          }while(1);
        }
      }

      /*Now Compute validation...*/
      matrix *q2y;
      matrix *sdep;
      initMatrix(&q2y);
      initMatrix(&sdep);

      if(rgcv_group == 0 || rgcv_iterations == 0){
        PLSLOOCV(xsample, ysample, xautoscaling, yautoscaling, nlv,
                  &q2y,
                  &sdep,
                  NULL,
                  NULL,
                  NULL, nthreads, s);
      }
      else{
        PLSRandomGroupsCV(xsample, ysample, xautoscaling, yautoscaling, nlv, rgcv_group, rgcv_iterations,
                      &q2y,
                      &sdep,
                      NULL,
                      NULL,
                      NULL, nthreads, s);
      }

      if(bestid != NULL){
        size_t pc = GetLVCCutoff(q2y);
        double cumq2 = 0.f;
        for(i = 0; i < q2y->col; i++)
          cumq2 += q2y->data[pc][i];

        if(cumq2 > bestq2){
          bestq2 = cumq2;
          for(i = 0; i < tmpid->size; i++)
            (*bestid)->data[i] = tmpid->data[i];
        }
      }

      //Search for the best Q2 Model and Select ID to return
      for(i = 0; i < q2y->col; i++){
        dvector *q2row = getMatrixColumn(q2y, i);
        MatrixAppendRow(&tmpq2, q2row);
        DelDVector(&q2row);

        dvector *sdep_row = getMatrixColumn(sdep, i);
        MatrixAppendRow(&tmpsdep, sdep_row);
        DelDVector(&sdep_row);
      }

      DelMatrix(&q2y);
      DelMatrix(&sdep);
    }
    DelMatrix(&xsample);
    DelMatrix(&ysample);

    /*Finalize output */
    ResizeMatrix(q2_distr, tmpq2->col, tmpq2->row);
    ResizeMatrix(sdep_distr, tmpsdep->col, tmpsdep->row);

    MatrixTranspose(tmpq2, (*q2_distr));
    MatrixTranspose(tmpsdep, (*sdep_distr));
    DelMatrix(&tmpsdep);
    DelMatrix(&tmpq2);
    DelUIVector(&tmpid);
  }
  else{
    fprintf(stderr, "Error!! The number of objects differ from the number of class objects\n");
  }
}

void PLSDynamicSampleValidator(matrix *mx, matrix *my,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t niters,
                        uivector *obj_class, size_t deltaobj, size_t maxobj,
                        size_t rgcv_group, size_t rgcv_iterations, size_t nthreads,
                        matrix **q2_surface, matrix **sdep_surface, uivector **bestid, ssignal *s)
{
  if(obj_class->size == mx->row && mx->row > 10){
    size_t iter, i, j, k, l, nclass = 0, srand_init, incobj;
    matrix *xsample, *ysample;
    uivector *tmpflag, *tmpbestflag;
    double bestq2 = -9999;
    for(i = 0; i < obj_class->size; i++){
      if(obj_class->data[i] > nclass){
        nclass = obj_class->data[i];
      }
      else{
        continue;
      }
    }

    if(maxobj > mx->row)
      maxobj = mx->row;

    NewMatrix(q2_surface, ceil(maxobj/deltaobj)*niters*nlv, 2+my->col);
    NewMatrix(sdep_surface, ceil(maxobj/deltaobj)*niters*nlv, 2+my->col);

    /*first colum: nobj
     second column: nlv
     from third to the end each serie of my->col is a q2 for each y*/

    srand_init = mx->row+mx->col+my->col+niters+nclass+maxobj;

    if(bestid != NULL){
      NewUIVector(&tmpflag, mx->row);
      NewUIVector(&tmpbestflag, mx->row);
    }

    l = 0;
    for(iter = 0; iter < niters; iter++){
      incobj = 10;
      for(k = 0; k < ceil(maxobj/deltaobj); k++){
        NewMatrix(&xsample, incobj, mx->col);
        NewMatrix(&ysample, incobj, my->col);

        if(nclass == 0){ /* simple bootstrap on all objects */
          for(i = 0; i < incobj; i++){
            do{
              size_t n = (size_t)myrand_r((unsigned *)&srand_init) % (mx->row);
              if(n > mx->row){
                continue;
              }
              else{
                if(bestid != NULL)
                  tmpflag->data[n] = 1;

                for(j = 0; j < mx->col; j++){
                  xsample->data[i][j] = mx->data[n][j];
                }

                for(j = 0; j < my->col; j++){
                  ysample->data[i][j] = my->data[n][j];
                }
                break;
              }
            }while(1);
          }
        }
        else{ /* select a random object from each class */
          size_t a_class = 0;
          for(i = 0; i < incobj; i++){
            do{
              size_t n = (size_t)myrand_r((unsigned *)&srand_init) % (mx->row);
              if(n > mx->row){
                continue;
              }
              else{
                if(obj_class->data[n] == a_class){
                  if(a_class == nclass){
                    a_class =  (size_t)myrand_r((unsigned *)&srand_init) % (nclass); /*all class are already completed then restart random...*/
                  }
                  else{
                    a_class++;
                  }

                  if(bestid != NULL)
                    tmpflag->data[n] = 1;

                  for(j = 0; j < mx->col; j++){
                    xsample->data[i][j] = mx->data[n][j];
                  }

                  for(j = 0; j < my->col; j++){
                    ysample->data[i][j] = my->data[n][j];
                  }
                  break;
                }
                else{
                  continue;
                }
              }
            }while(1);
          }
        }

        /*Now Compute validation...*/
        matrix *q2y;
        matrix *sdep;
        initMatrix(&q2y);
        initMatrix(&sdep);

        if(rgcv_group == 0 || rgcv_iterations == 0){
          PLSLOOCV(xsample, ysample, xautoscaling, yautoscaling, nlv,
                    &q2y,
                    &sdep,
                    NULL,
                    NULL,
                    NULL, nthreads, s);
        }
        else{
          PLSRandomGroupsCV(xsample, ysample, xautoscaling, yautoscaling, nlv, rgcv_group, rgcv_iterations,
                        &q2y,
                        &sdep,
                        NULL,
                        NULL,
                        NULL, nthreads, s);
        }

        if(bestid != NULL){
          size_t pc = GetLVCCutoff(q2y);
          double cumq2 = 0.f;
          for(i = 0; i < q2y->col; i++)
            cumq2 += q2y->data[pc][i];

          if(cumq2 > bestq2){
            bestq2 = cumq2;
            for(i = 0; i < tmpflag->size; i++){
              tmpbestflag->data[i] = tmpflag->data[i];
              tmpflag->data[i] = 0;
            }
          }
        }

        for(i = 0; i < q2y->row; i++){
          (*q2_surface)->data[l][0] = incobj;
          (*q2_surface)->data[l][1] = i+1;
          (*sdep_surface)->data[l][0] = incobj;
          (*sdep_surface)->data[l][1] = i+1;
          for(j = 0; j < q2y->col; j++){
            if(q2y->data[i][j] < -1){
              (*q2_surface)->data[l][j+2] = 0;
            }
            else
              (*q2_surface)->data[l][j+2] = q2y->data[i][j];
            (*sdep_surface)->data[l][j+2] = sdep->data[i][j];
          }
          l++;
        }

        if(incobj + deltaobj < mx->row)
          incobj += deltaobj;
        else
          incobj = mx->row;

        DelMatrix(&q2y);
        DelMatrix(&sdep);
        DelMatrix(&xsample);
        DelMatrix(&ysample);
      }
    }

    if(bestid != NULL){

      for(i = 0, j = 0; i < tmpbestflag->size; i++){
        if(tmpbestflag->data[i] == 1)
          j++;
        else
          continue;
      }

      UIVectorResize(bestid, j);
      for(i = 0, j = 0; i < tmpbestflag->size; i++){
        if(tmpbestflag->data[i] == 1){
          (*bestid)->data[j] = i;
          j++;
        }
        else
          continue;
      }

      DelUIVector(&tmpflag);
      DelUIVector(&tmpbestflag);
    }
  }
  else{
    if(mx->row > 10)
      fprintf(stderr, "Error!! Insufficient number of objects\n");
    else
      fprintf(stderr, "Error!! The number of objects differ from the number of class objects\n");
  }
}


int GetLVCCutoff(matrix *rq2y){
  size_t i, j, cutoff = 0;
  double prev, next, max = -9999;

  for(i = 0; i < rq2y->row-1; i += 2){
    prev = next = 0.f;
    for(j = 0; j < rq2y->col; j++){
      prev += rq2y->data[i][j];
      next += rq2y->data[i+1][j];
    }
    if(next > prev && next > max && (fabs(next-max)/next) > 0.03){ // do not select the component if the differences with the previous best is less than 3%!!
      cutoff = i+1;
      max = next;
    }
    else{
      if(prev > next && prev > max && (fabs(prev-max)/prev) > 0.03){
        cutoff = i;
        max = prev;
        break;
      }
      else{
        break;
      }
    }
  }
  return cutoff;
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

void PLSCostPopulation(matrix *mx, matrix *my,
                      matrix *px, matrix *py,
                      uivector *popvector,
                      size_t xautoscaling, size_t yautoscaling, size_t nlv,
                      int validation_type, size_t ngroup, size_t niter,
                      double *r2, double *q2, double *bias_, size_t nthreads, ssignal *s)
{
  size_t i, j, k, subsize, cutoff;
  double r2_, q2_, _bias_;
  matrix *submx, *subpx, *r2y, *sdec, *q2y, *sdep, *bias;
  PLSMODEL *m;


  subsize = 0;
  if(popvector != NULL){
    for(i = 0; i < popvector->size; i++){
      if(getUIVectorValue(popvector, i) == 1){
        subsize++;
      }
      else{
        continue;
      }
    }
    NewMatrix(&submx, mx->row, subsize);
    for(i = 0; i < mx->row; i++){
      k = 0;
      for(j = 0; j < mx->col; j++){
        if(getUIVectorValue(popvector, j) == 1){
          setMatrixValue(submx, i, k, getMatrixValue(mx, i, j));
          k++;
        }
        else{
          continue;
        }
      }
    }
  }
  else{
    initMatrix(&submx);
    MatrixCopy(mx, &submx);
    subsize = mx->col;
  }

  if(px != NULL && py != NULL && validation_type == 0){ /* cost = R^2*/
    NewPLSModel(&m);
    PLS(submx, my, nlv, xautoscaling, yautoscaling, m, s);

    NewMatrix(&subpx, px->row, subsize);

    k = 0;
    for(i = 0; i < px->row; i++){
      for(j = 0; j < px->col; j++){
        if(getUIVectorValue(popvector, j) == 1){
          setMatrixValue(subpx, i, k, getMatrixValue(px, i, j));
          k++;
        }
        else{
          continue;
        }
      }
    }

    initMatrix(&r2y);
    initMatrix(&sdec);
    PLSRSquared(subpx, py, m, nlv, &r2y, &sdec);

    cutoff = GetLVCCutoff(r2y);
//     MatrixGetMaxValue(r2y, &cutoff, NULL);
    r2_ = 0.f;
    for(j = 0; j < r2y->col; j++){
      r2_ += getMatrixValue(r2y, cutoff, j);
    }

    if(r2_ != 0){
      r2_ /= (double)r2y->col;
    }
    else{
      r2_ = 0.f;
    }

    if(r2!= NULL)
      (*r2) = (*q2) = r2_;
    else
      (*q2) = r2_;
    if(bias_ != NULL)
      (*bias_) = 0.f;

    DelMatrix(&sdec);
    DelMatrix(&r2y);
    DelMatrix(&subpx);
    DelPLSModel(&m);
  }
  else if(validation_type == 1){ /*LEAVE ONE OUT*/
    initMatrix(&q2y);
    initMatrix(&sdep);
    initMatrix(&bias);
    PLSLOOCV(submx, my, xautoscaling, yautoscaling, nlv,
                      &q2y,
                      &sdep,
                      &bias,
                      NULL, NULL, nthreads, s);

    NewPLSModel(&m);
    PLS(submx, my, nlv, xautoscaling, yautoscaling, m, s);

    PLSRSquared(submx, my,  m, nlv, &m->r2y_model, &m->sdec);

    cutoff = GetLVCCutoff(q2y);
//     MatrixGetMaxValue(q2y, &cutoff, NULL);

    q2_ = 0.f;
    for(j = 0; j < q2y->col; j++){
      q2_ += getMatrixValue(q2y, cutoff, j);
    }

    _bias_ = 0.f;
    for(j = 0; j < bias->col; j++){
      _bias_ += getMatrixValue(bias, cutoff, j);
    }

    r2_ = 0.f;
    for(j = 0; j < m->r2y_model->col; j++){
      r2_ += getMatrixValue(m->r2y_model, cutoff, j);
    }

    if(q2_ != 0){
      q2_ /= (double)q2y->col;
      r2_ /= (double)m->r2y_model->col;
      _bias_ /= (double)bias->col;
    }
    else{
      q2_ = 0.f;
      r2_ = 0.f;
      _bias_ = 99;
    }

    (*q2) = q2_;
    if(r2 != NULL)
      (*r2) = r2_;
    if(bias_ != NULL)
      (*bias_) = _bias_;

    DelPLSModel(&m);
    DelMatrix(&q2y);
    DelMatrix(&sdep);
    DelMatrix(&bias);
  }
  else{ /*CROSS VALIDATION*/
    initMatrix(&q2y);
    initMatrix(&sdep);
    initMatrix(&bias);
    PLSRandomGroupsCV(submx, my, xautoscaling, yautoscaling, nlv, ngroup, niter,
                  &q2y,
                  &sdep,
                  &bias,
                  NULL, NULL, nthreads, s);

    NewPLSModel(&m);
    PLS(submx, my, nlv, xautoscaling, yautoscaling, m, s);
    PLSRSquared(submx, my,  m, nlv, &m->r2y_model, &m->sdec);

    cutoff = GetLVCCutoff(q2y);
//     MatrixGetMaxValue(q2y, &cutoff, NULL);

    q2_ = 0.f;
    for(j = 0; j < q2y->col; j++){
      q2_ += getMatrixValue(q2y, cutoff, j);
    }

    _bias_ = 0.f;
    for(j = 0; j < bias->col; j++){
      _bias_ += getMatrixValue(bias, cutoff, j);
    }

    r2_ = 0.f;
    for(j = 0; j < m->r2y_model->col; j++){
      r2_ += getMatrixValue(m->r2y_model, cutoff, j);
    }

    if(q2_ != 0){
      q2_ /= (double)q2y->col;
      r2_ /= (double)m->r2y_model->col;
      _bias_ /= (double)bias->col;
    }
    else{
      q2_ = 0.f;
      r2_ = 0.f;
      _bias_ = 99;
    }

    (*q2) = q2_;
    if(r2 != NULL)
      (*r2) = r2_;

    if(bias_ != NULL)
      (*bias_) = _bias_;

    DelPLSModel(&m);
    DelMatrix(&q2y);
    DelMatrix(&sdep);
    DelMatrix(&bias);
  }

  DelMatrix(&submx);
}


void PSLPSOVariableSelection(matrix *mx, matrix *my, matrix *px, matrix *py,
                       size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type, size_t ngroup, size_t niter,
                       size_t population_size, double randomness,
                       uivector **varselected, matrix **map, uivector **vardistribution, size_t nthreads, ssignal *s)
{
  size_t i, j, iter, new_position, population_size_, tmp_var_id, nvarON;
  double model_r2, cost_g_best, cost_best, tmp_model_r2, tmp_cost_best, tmp_g_best, new_v, new_pos, rand1, rand2;
  matrix *models, *position, *velocity;
  uivector *popvector, *best_model;
  dvector *b_position, *g_position, *rowmap;

  /*1) Initialize the population*/

  population_size_ = population_size + 1; /* + 1 is the model with all variables */

  NewMatrix(&models, population_size_, mx->col);
  NewMatrix(&position, population_size_, mx->col); /*position chagned stored if are -1 then the position is not stored.... else is  a number >= 0 and < mx->col*/
  NewMatrix(&velocity, population_size_, mx->col);
  NewDVector(&b_position, mx->col);
  NewDVector(&g_position, mx->col);
  NewUIVector(&popvector, mx->col);
  NewUIVector(&best_model, mx->col);
  NewDVector(&rowmap, 4);

  UIVectorSet(popvector, 1); /*get all variables*/

  PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &model_r2, &cost_g_best, NULL, nthreads, s);

  setDVectorValue(rowmap, 0, model_r2);
  setDVectorValue(rowmap, 1, cost_g_best);
  setDVectorValue(rowmap, 2, (cost_g_best  * (population_size - mx->col -1)) /  (mx->col*(1-cost_g_best )));
  setDVectorValue(rowmap, 3, mx->col);

  MatrixAppendRow(map, rowmap);

  for(i = 0; i < mx->col; i++){
    UIVectorAppend(vardistribution, 1);
  }

  /*the best and global position are the same in the first model that have all the variable setted ON*/
  for(i = 0; i < mx->col; i++){
    /*All variable are global and best...*/
    setDVectorValue(g_position, i, i);
    setDVectorValue(b_position, i, i);
    setMatrixValue(velocity, 0, i, 0.f);
  }

  /*
  printf("model 0 with all variable: %f\n", cost_g_best);
  */

  /*2) Initialize the models particle.*/
  MatrixSet(models, 1.0); // all variables are setted ON!
  MatrixSet(position, -1.0); // all variables on model are unvariate...
  MatrixSet(velocity, 0); // all variables have no speed...

  srand(mx->row + mx->col + my->col + population_size_);

  for(i = 1; i < population_size_; i++){

    for(j = 0; j < (size_t)floor(mx->col * randomness); j++){ /* % of object to randomize in this case is the 20 % of variable*/
      setMatrixValue(velocity, i, j, (double)rand()/RAND_MAX); /*velocity is a number between 0 an 1*/

      new_position = rand() % mx->col;

      setMatrixValue(position, i, j, new_position);

      if((size_t)getMatrixValue(models, i, new_position) == 1){ /* set to OFF*/
        setMatrixValue(models, i, new_position, 0.0);
      }
      else{ /* Set to ON*/
        setMatrixValue(models, i, new_position, 1.0);
      }
    }

    /*COMPUTE THE COST OF ACTUAL MODEL RANDOM MODEL*/
    for(j = 0; j < models->col; j++){
      setUIVectorValue(popvector, j, (size_t)getMatrixValue(models, i, j));
    }

    nvarON = 0;
    for(j = 0; j < models->col; j++){
      if(getUIVectorValue(popvector, j) == 1){
        nvarON++;
      }
    }

    if(nvarON > 0){
      PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &tmp_model_r2, &cost_best, NULL, nthreads, s);
    }
    else{
      tmp_model_r2 = cost_best = 0;
    }

    for(j = 0; j < models->col; j++){
      setDVectorValue(b_position, j, getMatrixValue(position, i, j));
    }

    if(cost_best > cost_g_best){
      cost_g_best = cost_best;
      model_r2 = tmp_model_r2;

      for(j = 0; j < models->col; j++){
        setDVectorValue(g_position, j, getMatrixValue(position, i, j));
        setUIVectorValue(best_model, j, getUIVectorValue(popvector, j));
      }

      setDVectorValue(rowmap, 0, model_r2);
      setDVectorValue(rowmap, 1, cost_g_best);
      setDVectorValue(rowmap, 2, (cost_g_best  * (population_size - mx->col -1)) /  (mx->col*(1-cost_g_best )));
      setDVectorValue(rowmap, 3, nvarON);

      MatrixAppendRow(map, rowmap);
      DVectorSet(rowmap, 0.f);
    }

    /*
    printf("Modello %u Generato\n", (size_t)i);
    PrintUIVector(popvector);
    */
  }


  /*
  puts("modelli");
  PrintMatrix(models);

  puts("Best Position");
  PrintDVector(b_position);
  puts("Great Best Position");
  PrintDVector(g_position);

  puts("Second step!");
  */


/* While(R^2 or Q^2 is not good enough...):
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
 */

  iter = 0;
  do{
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      /*
      puts("Velocity");
      PrintMatrix(velocity);
      puts("Positions");
      PrintMatrix(position);
      while (!getchar()=='\n');
      */

      tmp_g_best = cost_g_best;

      for(i = 0; i < models->row; i++){
        /*
        * update velocity and position and set variable on/off in model
        */
        for(j = 0; j < velocity->col; j++){
            rand1 = (double)rand()/(double)RAND_MAX;
            rand2 = (double)rand()/(double)RAND_MAX;
            new_v = getMatrixValue(velocity, i, j);
            new_v += (1 * rand1 * ((double)getDVectorValue(b_position, j) - getMatrixValue(position, i, j))); /*1 is a constant weight and could be variable */
            new_v += (1 * rand2 * ((double)getDVectorValue(g_position, j) - getMatrixValue(position, i, j))); /*1 is a constant weight and could be variable */
            setMatrixValue(velocity, i, j, new_v);

            /*
            * p(t+1) = p(t) + v(t); this value must be positive...
            */
            new_pos = fabs(floor(getMatrixValue(position, i, j) + new_v));
            setMatrixValue(position, i, j, new_pos);
            tmp_var_id = (size_t)getMatrixValue(position, i, j);

            if(new_pos > models->col-1 || new_pos < 0){
              /*
              printf("Error on setting the variable. Out of range: %f\n", getMatrixValue(position, i, j));
              */
              continue;
            }
            else{
              if((size_t)getMatrixValue(models, i, tmp_var_id) == 1){ /* set to OFF*/
                setMatrixValue(models, i, tmp_var_id, 0.0);
              }
              else{ /* Set to ON*/
                setMatrixValue(models, i, tmp_var_id, 1.0);
              }
            }
        }

        for(j = 0; j < models->col; j++){
          setUIVectorValue(popvector, j, (size_t)getMatrixValue(models, i, j));
          setUIVectorValue((*vardistribution), j, (getUIVectorValue((*vardistribution), j)+(size_t)getMatrixValue(models, i, j)));
        }


        nvarON = 0;
        for(j = 0; j < position->col; j++){
          if(getUIVectorValue(popvector, j) == 1){
            nvarON++;
          }
        }

        if(nvarON > 0){
          PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &tmp_model_r2, &tmp_cost_best, NULL, nthreads, s);
        }
        else{
          tmp_model_r2 = tmp_cost_best = 0;
        }

        if(tmp_cost_best > cost_best){
          cost_best = tmp_cost_best;

          for(j = 0; j < position->col; j++){
            setDVectorValue(b_position, j, fabs(floor(getMatrixValue(position, i, j))));
          }

          if(cost_best > cost_g_best){
            cost_g_best = cost_best;
            model_r2 = tmp_model_r2;

            for(j = 0; j < position->col; j++){
              setDVectorValue(g_position, j, fabs(floor(getMatrixValue(position, i, j))));
              setUIVectorValue(best_model, j, getUIVectorValue(popvector, j));
            }

            setDVectorValue(rowmap, 0, model_r2);
            setDVectorValue(rowmap, 1, cost_g_best);
            setDVectorValue(rowmap, 2, (cost_g_best  * (population_size - mx->col -1)) /  (mx->col*(1-cost_g_best )));
            setDVectorValue(rowmap, 3, nvarON);

            MatrixAppendRow(map, rowmap);

            DVectorSet(rowmap, 0.f);

            iter = 0;
          }
        }

        /*
        printf("iteration %u best_val %f global_best_value %f\n", (size_t) iter, cost_best, cost_g_best);
        */
        iter++;
      }
    }
  }while(FLOAT_EQ(cost_g_best, tmp_g_best, PLSCONVERGENCE) && iter < 100);

  /*
  PrintDVector(g_position);

  printf("Selected Model: %u\n", (size_t)mod_id);
  PrintMatrix(models);
  */

  for(i = 0; i < best_model->size; i++){
    UIVectorAppend(varselected, getUIVectorValue(best_model, i));
  }

  DelUIVector(&best_model);
  DelDVector(&rowmap);
  DelDVector(&g_position);
  DelDVector(&b_position);
  DelUIVector(&popvector);
  DelMatrix(&velocity);
  DelMatrix(&position);
  DelMatrix(&models);
}


/* Variable selection using Genetic algorithm
 *
 *
 */

int FitnessMaxsimized(dvector* modelsfitness, dvector* modelsfitness_old, double populationconvergence)
{
  size_t i, popsize = 0;

  /*
  puts("OLD");
  PrintDVector(modelsfitness_old);
  puts("NEW");
  PrintDVector(modelsfitness);

  sleep(1);
  */

  for(i = 0; i < modelsfitness->size; i++){
    if(FLOAT_EQ(getDVectorValue(modelsfitness, i), getDVectorValue(modelsfitness_old, i), 1e-3)){
      continue;
    }
    else{
      popsize++;
    }
  }

  if((double)popsize/(double)modelsfitness->size < populationconvergence){
    return 1;
  }
  else{
    return 0;
  }
}

void PLSGAVariableSelection(matrix *mx, matrix *my, matrix *px, matrix *py,
                       size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type, size_t ngroup, size_t niter,
                       size_t population_size, double fraction_of_population, double mutation_rate, size_t crossovertype, double nswapping, double populationconvergence,
                      uivector **varselected, matrix **map, uivector **vardistribution, size_t nthreads, ssignal *s)
{
  size_t i, j, k, position, mutatebit, bitswap, totbitswapp, copy_selection, crossover_selection, mutation_selection, init, a, b, nvarON, blockitercount;
  double best_fitness = 0, best_bias = 99, model_r2, tmp_fitness, tmp_bias, tmp_model_r2;
  matrix *models1, *models2;
  dvector *modelsfitness, *best_modelsfitness, *model_f_distribution, *best_model_f_distribution, *rowmap;
  uivector *popvector, *best_model, *crosspoints, *selectedid;

  /* Initialize with rando models..*/
  NewMatrix(&models1, population_size, mx->col); /*models of the generation k */
  NewMatrix(&models2, population_size, mx->col); /*models of the generation k + 1*/
  NewDVector(&modelsfitness, population_size); /* Stored the Q^2 or R^2 of models*/
  NewDVector(&best_modelsfitness, population_size); /* Stored the Q^2 or R^2 of models*/
  NewDVector(&model_f_distribution, population_size); /* Stored the Q^2 or R^2 of models*/
  NewDVector(&best_model_f_distribution, population_size); /* Stored the Q^2 or R^2 of models*/
  NewUIVector(&best_model, mx->col);
  NewUIVector(&popvector, mx->col);

  for(i = 0; i < mx->col; i++){
    UIVectorAppend(vardistribution, 0);
  }

  NewDVector(&rowmap, 4); /* 1 = R2, 2 = Q2, 3 = f-test, 4 = number of variables ON */
  best_fitness = 0.f;

  copy_selection = (size_t)floor((1-fraction_of_population) * population_size);
  crossover_selection = (size_t)floor(fraction_of_population * population_size);
  mutation_selection = (size_t)floor(mutation_rate*population_size);

  if(crossovertype == 0){ /*single point crossover */
    totbitswapp = 1;
  }
  else if(crossovertype == 1){ /* two point crossover */
    totbitswapp = 2;
  }
  else{ /*uniform crossover*/
    if(nswapping > 1){
      totbitswapp = (size_t)floor(((double)models1->col * nswapping)/ 100.f);
    }
    else{
      totbitswapp = (size_t)floor((double)models1->col * nswapping);
    }

    if(totbitswapp > 0){
      if(totbitswapp % 2 != 0){
        totbitswapp +=1;
      }
    }
    else{
      totbitswapp = 1;
    }
  }

  NewUIVector(&crosspoints, totbitswapp);

  /*printf("copy %u cross %u muta %u bit to swap: %u\n", (size_t)copy_selection, (size_t)crossover_selection, (size_t)mutation_selection, (size_t)totbitswapp);*/
  init = (size_t)(mx->row + mx->col + my->col + population_size + copy_selection + mutation_selection);
  srand(init);

  DVectorSet(modelsfitness, 0.f);
  DVectorSet(best_modelsfitness, 0.f);
  DVectorSet(model_f_distribution, 0.f);
  DVectorSet(best_model_f_distribution, 0.f);
  MatrixSet(models1, 0.f);
  MatrixSet(models2, 0.f);

  for(i = 0; i < population_size; i++){

    for(j = 0; j < mx->col; j++){
      position = rand() % mx->col;
      if((size_t)getMatrixValue(models1, i, position) == 1){ /* set to OFF*/
        setMatrixValue(models1, i, position, 0.0);
      }
      else{ /* Set to ON*/
        setMatrixValue(models1, i, position, 1.0);
      }
    }

    /*COMPUTE THE COST OF ACTUAL MODEL RANDOM MODEL*/
    for(j = 0; j < models1->col; j++){
      setUIVectorValue(popvector, j, (size_t)getMatrixValue(models1, i, j));
      setUIVectorValue((*vardistribution), j, getUIVectorValue((*vardistribution),j)+(size_t)getMatrixValue(models1, i, j));
    }

    nvarON = 0;
    for(j = 0; j < popvector->size; j++){
      if(getUIVectorValue(popvector, j) == 1){
        nvarON++;
      }
    }

    if(nvarON > 0){
      PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &tmp_model_r2, &tmp_fitness, &tmp_bias, nthreads, s);
    }
    else{
      tmp_model_r2 = tmp_fitness = 0;
    }

    setDVectorValue(modelsfitness, i, tmp_fitness);
    setDVectorValue(model_f_distribution, i, (tmp_fitness * (population_size - mx->col -1)) /  (mx->col*(1-tmp_fitness)));

    if(tmp_fitness > best_fitness && tmp_bias < best_bias){
      best_fitness = tmp_fitness;
      model_r2 = tmp_model_r2;
      best_bias = tmp_bias;

      for(j = 0; j < models1->col; j++){
        setUIVectorValue(best_model, j, getUIVectorValue(popvector, j));
      }

      setDVectorValue(rowmap, 0, model_r2);
      setDVectorValue(rowmap, 1, best_fitness);
      setDVectorValue(rowmap, 2, getDVectorValue(model_f_distribution, i));
      setDVectorValue(rowmap, 3, (double)nvarON);

      MatrixAppendRow(map, rowmap);
      DVectorSet(rowmap, 0.f);
    }


    /*
    printf("Modello %u Generato\n", (unsigned int)i);
    PrintUIVector(popvector);
    */
  }

  /*
  puts("first best model");
  PrintUIVector(best_model);

  PrintMatrix((*map));

  puts("Models1");
  PrintMatrix(models1);
  puts("Models2");
  PrintMatrix(models2);
  */

  blockitercount = 0;
  do{
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      /*printf("blockitercount %u\n", (unsigned int) blockitercount);*/

      tmp_fitness = tmp_model_r2 = 0.f;
      tmp_bias = 99;
      /*Create generation k+1*/
      /*1. Copy
        *
        * Select (1-fraction_of_population) * population_size members of models1 and insert into the new population models2
        * through the roulette whell selection
        */
      if(crossover_selection % 2 != 0){
        crossover_selection -=1;
        copy_selection +=1;
      }

      initUIVector(&selectedid);

      StochasticUniversalSample(modelsfitness, copy_selection, init, &selectedid);
      /*RouletteWheelselection(modelsfitness, copy_selection, init, &selectedid);*/

      k = 0;
      for(i = 0; i < selectedid->size; i++, k++){
        position = getUIVectorValue(selectedid, i);
        for(j = 0; j < models1->col; j++){
          setMatrixValue(models2, i, j, getMatrixValue(models1, position, j));
        }
      }
      DelUIVector(&selectedid);


      /*
      * 2. Crossover
      *
      * Select fraction_of_population * population_size members of models1; pair them up; produce offspring; inser the offspring into models2
      *
      * The selection is made with the roulette whell selection
      *
      *
      */

      initUIVector(&selectedid);

      StochasticUniversalSample(modelsfitness, crossover_selection, init, &selectedid);
      /*RouletteWheelselection(modelsfitness, crossover_selection, init, &selectedid);*/

      for(i = 0; i < selectedid->size; i+=2){

        position = getUIVectorValue(selectedid, i);

        for(bitswap = 0; bitswap < totbitswapp; bitswap++){
          setUIVectorValue(crosspoints, bitswap, rand() % models1->col);
        }

        SortUIVector(crosspoints);

        if(totbitswapp >=2){
          /*offspring 1*/
          for(j = 0; j < models1->col; j++){
            for(bitswap = 0; bitswap < totbitswapp; bitswap += 2){
              if(j > getUIVectorValue(crosspoints, bitswap) && j < getUIVectorValue(crosspoints, bitswap+1)){ // set pos1
                setMatrixValue(models2, k, j, getMatrixValue(models1, position, j));
              }
              else{ // set pos 2
                setMatrixValue(models2, k, j, getMatrixValue(models1, position+1, j));
              }
            }
          }

          /*offspring 2*/
          k++;
          for(j = 0; j < models1->col; j++){
            for(bitswap = 0; bitswap < totbitswapp; bitswap += 2){
              if(j > getUIVectorValue(crosspoints, bitswap) && j < getUIVectorValue(crosspoints, bitswap+1)){ // set pos1
                setMatrixValue(models2, k, j, getMatrixValue(models1, position+1, j));
              }
              else{ // set pos 2
                setMatrixValue(models2, k, j, getMatrixValue(models1, position, j));
              }
            }
          }
        }
        else{
          /*offspring 1*/
          for(j = 0; j < models1->col; j++){
            if(j < getUIVectorValue(crosspoints, 0)){
              setMatrixValue(models2, k, j, getMatrixValue(models1, position, j));
            }
            else{
              setMatrixValue(models2, k, j, getMatrixValue(models1, position+1, j));
            }
          }

          /*offspring 2*/
          k++;
          for(j = 0; j < models1->col; j++){
            if(j < getUIVectorValue(crosspoints, 0)){
              setMatrixValue(models2, k, j, getMatrixValue(models1, position+1, j));
            }
            else{
              setMatrixValue(models2, k, j, getMatrixValue(models1, position, j));
            }
          }
        }
        k++;

      }
      DelUIVector(&selectedid);

    /*3. Mutate
      * Select a mutation_rate * population_size elements of models2 and invert a randomly selected bit in each;
      * The selection is always made by the roulette Whell Selection algorithm
      */

      initUIVector(&selectedid);

      StochasticUniversalSample(modelsfitness, mutation_selection, init, &selectedid);
  /*      RouletteWheelselection(modelsfitness, mutation_selection, init, &selectedid);*/

      for(i = 0; i < selectedid->size; i++){

        position = getUIVectorValue(selectedid, i);
        mutatebit = rand() % models1->col;

        if((size_t)getMatrixValue(models2, position, mutatebit) == 1){ /*set OFF*/
          setMatrixValue(models2, position, mutatebit, 0.0);
        }
        else{ /*set ON*/
          setMatrixValue(models2, position, mutatebit, 1.0);
        }
      }
      DelUIVector(&selectedid);


    /* Evaluate fittness for the new generation.....
      *
      * Copy the old fitness to
      *
      */

      for(i = 0; i < modelsfitness->size; i++){
        if(getDVectorValue(modelsfitness, i) > getDVectorValue(best_modelsfitness, i)){
          setDVectorValue(best_modelsfitness, i, getDVectorValue(modelsfitness, i));
        }

        if(getDVectorValue(model_f_distribution, i) > getDVectorValue(best_model_f_distribution, i)){
          setDVectorValue(best_model_f_distribution, i, getDVectorValue(model_f_distribution, i));
        }
      }

      for(i = 0; i < models2->row; i++){
        for(j = 0; j < models2->col; j++){
          setUIVectorValue(popvector, j, (size_t)getMatrixValue(models2, i, j));
          setUIVectorValue((*vardistribution), j, getUIVectorValue((*vardistribution),j)+(size_t)getMatrixValue(models1, i, j));
        }

        nvarON = 0;
        for(j = 0; j < models2->col; j++){
          if(getUIVectorValue(popvector, j) == 1){
            nvarON++;
          }
        }

        if(nvarON > 0){
          PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &tmp_model_r2, &tmp_fitness, &tmp_bias, nthreads, s);
        }
        else{
          tmp_model_r2 = tmp_fitness = 0;
          tmp_bias = 99;
        }

        setDVectorValue(modelsfitness, i, tmp_fitness); /* Stored the Q^2 or R^2 of models*/
        setDVectorValue(model_f_distribution, i, (tmp_fitness * (population_size - mx->col -1)) /  (mx->col*(1-tmp_fitness)));

        if(tmp_fitness > best_fitness && tmp_bias < best_bias){

          best_bias = tmp_bias;
          best_fitness = tmp_fitness;
          model_r2 = tmp_model_r2;

          for(j = 0; j < models2->col; j++){
            setUIVectorValue(best_model, j, getUIVectorValue(popvector, j));
          }

          /*
          puts("POPVECTOR");
          PrintUIVector(popvector);

          puts("Models1");
          PrintMatrix(models1);
          puts("Models2");
          PrintMatrix(models2);

          sleep(1);
          */
          setDVectorValue(rowmap, 0, model_r2);
          setDVectorValue(rowmap, 1, best_fitness);
          setDVectorValue(rowmap, 2, getDVectorValue(model_f_distribution, i));
          setDVectorValue(rowmap, 3, nvarON);

          MatrixAppendRow(map, rowmap);
          DVectorSet(rowmap, 0.f);
        }
      }

      /* copy the new generation... */
      for(i = 0; i < models2->row; i++){
        for(j = 0; j < models2->col; j++){
          setMatrixValue(models1, i, j, getMatrixValue(models2, i, j));
        }
      }

      a = FitnessMaxsimized(modelsfitness, best_modelsfitness, populationconvergence);
      b = FitnessMaxsimized(model_f_distribution, best_model_f_distribution, populationconvergence);

      /*
      printf("a %d  b %d\n", (unsigned int)a, (unsigned int)b);
      */

      if(a == b && a == 1){
        if(blockitercount > 1000){
          /*No more good solutions found...
          * decrease the population size and the population convergence
          */
          break;
        }
        else{
          blockitercount++;
        }
      }
      else{
        blockitercount = 0;
      }
    }
  }while((a != 0 && b != 0));


  /*
  puts("BEST_MODEL");
  PrintUIVector(best_model);
  */

  for(i = 0; i < best_model->size; i++){
    UIVectorAppend(varselected, getUIVectorValue(best_model,  i));
  }

  /*
  puts("VARSELEVTED");
  PrintUIVector((*varselected));
  puts("----------");
  */

  DelDVector(&rowmap);
  DelUIVector(&crosspoints);
  DelDVector(&best_modelsfitness);
  DelDVector(&model_f_distribution);
  DelDVector(&best_model_f_distribution);
  DelDVector(&modelsfitness);
  DelUIVector(&popvector);
  DelUIVector(&best_model);
  DelMatrix(&models2);
  DelMatrix(&models1);
}

void PLSSpearmannVariableSelection(matrix* mx, matrix* my, matrix* px, matrix* py,
                                   size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type, size_t ngroup, size_t niter,
                                   double threshold,
                                   uivector** varselected, matrix** map, uivector **vardistribution, size_t nthreads, ssignal* s){
  size_t i, j, k, nstep, nvarON;
  double step, r2, q2, best_q2 = 0, rangemax;

  matrix *mxy, *spearmanncorrelmx;
  uivector *popvector, *bestmodel;
  dvector *rowmap;

  NewDVector(&rowmap, 4);
  NewUIVector(&popvector, mx->col);
  NewUIVector(&bestmodel, mx->col);
  NewMatrix(&mxy, mx->row, mx->col+my->col);
  UIVectorResize(vardistribution, mx->col);

  for(i = 0; i < mx->row; i++){
    k = 0;
    for(j = 0; j < mx->col; j++){
      setMatrixValue(mxy, i, k, getMatrixValue(mx, i, j));
      k++;
    }

    for(j = 0; j < my->col; j++){
      setMatrixValue(mxy, i, k, getMatrixValue(my, i, j));
      k++;
    }
  }

  initMatrix(&spearmanncorrelmx);

  SpearmanCorrelMatrix(mxy, spearmanncorrelmx);

  rangemax = fabs(getMatrixValue(spearmanncorrelmx, 0, mx->col)); /* the first y value in the spearmanncorrelmx */
  for(j = 0; j < my->col; j++){
    for(i = 1; i < mx->col; i++){
      if(rangemax < getMatrixValue(spearmanncorrelmx, i, mx->col+j)){
        rangemax = getMatrixValue(spearmanncorrelmx, i, mx->col+j);
      }
      else{
        continue;
      }
    }
  }

  if(threshold < 0 || threshold > 1){
    threshold = 0.1; /* 1 / 10 */
  }

  threshold *= rangemax;
  nstep = (size_t) ceil(rangemax/threshold);
  step = rangemax - threshold; /* From maximum of relaction to minimum */
  for(i = 0; i < nstep; i++){
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      nvarON = 0;
      for(k = 0; k < mx->col; k++){
        double nyes = 0.f;
        for(j = 0; j < my->col; j++){
          if(fabs(getMatrixValue(spearmanncorrelmx, mx->col+j, k)) > step){
            nyes += 1;
          }
          else{
            continue;
          }
        }

        if((double)((my->col - nyes)/(double)my->col) > 0.5){ /* if the variable is no more good for the 50% of y reject*/
          setUIVectorValue(popvector, k, 0);
        }
        else{
          setUIVectorValue(popvector, k, 1);
          setUIVectorValue((*vardistribution), k, getUIVectorValue((*vardistribution), k)+1);
          nvarON++;
        }
      }

      if(nvarON > 0){
        PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &r2, &q2, NULL, nthreads, s);


        if(q2 > best_q2){
          for(j = 0; j < popvector->size; j++){
            setUIVectorValue(bestmodel, j, getUIVectorValue(popvector, j));
          }
          best_q2 = q2;
        }

        setDVectorValue(rowmap, 0, r2);
        setDVectorValue(rowmap, 1, q2);
        setDVectorValue(rowmap, 2, step);
        setDVectorValue(rowmap, 3, nvarON); /* number of variables */
        MatrixAppendRow(map, rowmap);
        DVectorSet(rowmap, 0);
        UIVectorSet(popvector, 0);
      }
      step -= threshold;
    }
  }

  for(i = 0; i < bestmodel->size; i++){
    UIVectorAppend(varselected, getUIVectorValue(bestmodel,  i));
  }

  DelUIVector(&bestmodel);
  DelDVector(&rowmap);
  DelUIVector(&popvector);
  DelMatrix(&spearmanncorrelmx);
  DelMatrix(&mxy);
}

void PrintPLSModel(PLSMODEL* model)
{
  puts("T Scores");
  PrintMatrix(model->xscores);
  puts("U Scores");
  PrintMatrix(model->yscores);

  puts("P Loadings");
  PrintMatrix(model->xloadings);
  puts("Q Loadings");
  PrintMatrix(model->yloadings);

  puts("W Weight");
  PrintMatrix(model->xweights);

  puts("Coefficient Regression for the relaction u = b * t");
  PrintDVector(model->b);

  puts("X Variance Explained");
  PrintDVector(model->xvarexp);

  puts("Recalculated y");
  PrintMatrix(model->recalculated_y);

  puts("Recalculated Residuals");
  PrintMatrix(model->recalc_residuals);
}

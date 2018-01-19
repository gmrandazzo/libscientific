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
#include "array.h"
#include "matrix.h"
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
void LVCalc(matrix **X, matrix **Y, dvector **t, dvector **u, dvector **p, dvector **q, dvector **w, double *bcoef)
{
  size_t i, j, loop;
  double mod_p_old, dot_q, dot_t, dot_u, dot_w, deltat;
  mod_p_old = dot_q = dot_t = dot_u = dot_w = 0.f;

  dvector *t_old;
  NewDVector(&t_old, (*t)->size);
  /* Step 1: select the column vector u with the largest column average from Y  take u = some y_j */
  if((*Y)->col > 1){
    dvector *Y_avg;
    initDVector(&Y_avg);
    MatrixColVar((*Y), &Y_avg);
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
  for(i = 0; i < (*u)->size; i++)
    (*u)->data[i] = (*Y)->data[i][j];
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
    /* Step 2: Compute a weight vector w = u'X/u'u   (NB.: w is size of X_t->col = (*X)->row ) */
    DVectorSet((*w), 0.f); /* Reset the w vector */
    DVectorMatrixDotProduct((*X), (*u), (*w));

    dot_u = DVectorDVectorDotProd((*u), (*u));

    for(i = 0; i < (*w)->size; i++){
      (*w)->data[i] /= dot_u;
    }
    /* End Step2 */

    /* Step 3: w = w/||w1||   Normalize */
    DVectNorm((*w), (*w));

    /* End Step 3 */

    #ifdef DEBUG
    printf("\n w Weight Vector for X matrix \n");
    PrintDVector((*w));
    #endif

    /* Step 4:   t = Xw/w'w Compute a t score vector. t i size of (*X)->row */
    DVectorSet((*t), 0.f); /* Reset the t vector */
    MatrixDVectorDotProduct((*X), (*w), (*t));

    dot_w = DVectorDVectorDotProd((*w), (*w));

    for(i = 0; i < (*t)->size; i++){
      (*t)->data[i] /= dot_w;
    }

    #ifdef DEBUG
    printf("\n t Score Vector for X matrix \n");
    PrintDVector(t);
    #endif
    /* End Step 4 */

    if((*Y)->col > 1){
      /* Step 5: Compute the loadings Y vector q' = t'Y/t't  and normalize q = q/||q||*/
      DVectorSet((*q), 0.f); /* Reset the c vector */
      DVectorMatrixDotProduct((*Y), (*t), (*q));

      dot_t = DVectorDVectorDotProd((*t), (*t));

      for(i = 0; i < (*q)->size; i++){
        (*q)->data[i] /= dot_t;
      }

      DVectNorm((*q), (*q));
      /* End Step 5 */

      /* Step 6. Update the Y Scores u  u = Yq/q'q */
      DVectorSet((*u), 0.f);
      MatrixDVectorDotProduct((*Y), (*q), (*u));
      dot_q = DVectorDVectorDotProd((*q), (*q));

      for(i = 0; i < (*u)->size; i++){
        (*u)->data[i] /= dot_q;
      }

      /* End step 6 */

      #ifdef DEBUG
      printf("\n New Score vector u for Y\n");
      PrintDVector(u);
      #endif
    }
    else{
      /*If the Y block has only one variable, step 5-8 can be omitted by putting q = 1 and no more iteration is necessary.*/
      DVectorSet((*q), 1);
      //break;
    }

    /* Step 8 if t_old == t new with a precision PLSCONVERGENCE then stop iteration else restart from step 2 */

    if(loop == 0){
      for(i = 0; i < (*t)->size; i++)
        t_old->data[i] = (*t)->data[i];
    }
    else{
      deltat = 0.f;
      for(i = 0; i < (*t)->size; i++){
        deltat = ((*t)->data[i] - t_old->data[i]) * ((*t)->data[i] - t_old->data[i]);
        t_old->data[i] = (*t)->data[i];
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
  DVectorMatrixDotProduct((*X), (*t), (*p));

  dot_t = DVectorDVectorDotProd((*t), (*t));
  for(i = 0; i < (*p)->size; i++){
    (*p)->data[i] /= dot_t;
  }

  #ifdef DEBUG
  printf("X Loadings\n");
  PrintDVector(p);
  #endif
  /* End Step 9*/

  mod_p_old = DvectorModule((*p));

  /* Step 10  p_new = p_old'/||p_old'|| */
  DVectNorm((*p), (*p));

  /*Step 11 t = t||p'||*/
  for(i = 0; i < (*t)->size; i++){
    (*t)->data[i] *= mod_p_old;
  }

  /*Step 12 */
  for(i = 0; i < (*w)->size; i++){
    (*w)->data[i] *= mod_p_old;
  }

  /* Step 13 Find the Regression Coefficient b:  b = u't/t't */
  dot_t = DVectorDVectorDotProd((*t), (*t));

  (*bcoef) = DVectorDVectorDotProd((*u), (*t))/dot_t;

  /*End Step 13*/

  /* Step 14  Adjust X for what has been found: Xnew=X-tp'  and Y fort what has been found Ynew = Y - btp' */
  /* X */
  for(i = 0; i < (*X)->row; i++){
    for(j = 0; j < (*X)->col; j++){
      (*X)->data[i][j] -= (*t)->data[i]*(*p)->data[j];
    }
  }

  /* Y */
  for(i = 0; i < (*Y)->row; i++){
    for(j = 0; j < (*Y)->col; j++){
      (*Y)->data[i][j] -= ((*bcoef) * (*t)->data[i] * (*q)->data[j]);
    }
  }

  #ifdef DEBUG
  printf("\n Deflated X \n");
  PrintMatrix(X);
  printf("\n Deflated Y\n");
  PrintMatrix(Y);
  #endif
  DelDVector(&t_old);
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

    size_t i, j, pc; /* pc = principal component */

    matrix *X; /* data matrix of mean centred object  for mx */
    matrix *Y; /* data matrix of mean centred object  for my */

    dvector *t; /* X score vector*/
    dvector *p; /* X loading vector */
    dvector *w; /* X weights vector */

    dvector *u; /* Y score vector*/
    dvector *q; /* Y loadings vector*/

    dvector *xeval;
    dvector *yeval;

    double min, max, ssx;

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
          double bcoef = 0.f;
          /* Calculate the Latent Variable (LV) according the NIPALS algorithm */
          LVCalc(&X, &Y, &t, &u, &p, &q, &w, &bcoef);

          #ifdef DEBUG
          printf("\n Deflated X\n");
          PrintMatrix(X);
          printf("\n Deflated Y\n");
          PrintMatrix(Y);
          #endif

          /* Calculating eigenvalue for X and Y in order to estimate the explained variance for each component
          * module of u vector and module of t vector are the eigenvalue
          */
          yeval->data[pc] = DVectorDVectorDotProd(u,u); /* NO MAKE SENSE EXPLAIN THE VARIANCE FOR Y*/
          xeval->data[pc] = DVectorDVectorDotProd(t,t); /* Adding the t eigenvalue */

          /* If more pairs (t,u) are needed go to 1. with X=Xnew and Y=Ynew.
          * Storing scores, loadings, weights, bcoefficients
          */
          DVectorAppend(&(model->b), bcoef);

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


/* Estimate the y dependent variable at nlv
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
  size_t lv, i, j;
  double _b;


  if(nlv > tscore->col)
    nlv = tscore->col;

  /*Allocate the y results*/
  ResizeMatrix(y, tscore->row, model->yloadings->row);

  for(lv = 0; lv < nlv; lv++){
    _b = model->b->data[lv];

    for(i = 0; i < tscore->row; i++){
      for(j = 0; j < model->yloadings->row; j++){
        (*y)->data[i][j] += _b * tscore->data[i][lv] * model->yloadings->data[j][lv];
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

void PLSYPredictorAllLV(matrix *mx, PLSMODEL *model, matrix **tscores, matrix **y)
{
  size_t lv, i, j, n_y, nlv;
  matrix *predicted_y;
  matrix *predicted_xscores;

  if(tscores == NULL){
    initMatrix(&predicted_xscores);
  }
  else{
    predicted_xscores = (*tscores);
  }

  nlv = model->b->size;

  n_y = model->yloadings->row;
  ResizeMatrix(y, mx->row, n_y*nlv);

  PLSScorePredictor(mx, model, nlv, &predicted_xscores);

  for(lv = 0; lv < nlv; lv++){
    initMatrix(&predicted_y);
    PLSYPredictor(predicted_xscores, model, lv+1, &predicted_y);

    for(i = 0; i < mx->row; i++){
      for(j = 0; j < n_y; j++){
        (*y)->data[i][n_y*lv+j] = predicted_y->data[i][j];
      }
    }
    DelMatrix(&predicted_y);
  }

  if(tscores == NULL){
    DelMatrix(&predicted_xscores);
  }
}

void PLSRegressionStatistics(matrix *my_true, matrix *my_pred, matrix** ccoeff, matrix **stdev, matrix **bias)
{
  size_t lv, i, j;
  dvector *ymean;

  size_t nlv = (size_t) my_pred->col/my_true->col;
  /*Calculate the Q2 and SDEP */
  if(ccoeff != NULL)
    ResizeMatrix(ccoeff, nlv, my_true->col);

  if(stdev != NULL)
    ResizeMatrix(stdev, nlv, my_true->col);

  if(bias != NULL)
    ResizeMatrix(bias, nlv, my_true->col);

  initDVector(&ymean);
  MatrixColAverage(my_true, &ymean);

  for(lv = 0; lv < nlv; lv++){
    for(j = 0; j < my_true->col; j++){
      double ssreg = 0.f;
      double sstot = 0.f;
      for(i = 0; i < my_true->row; i++){
        ssreg += square(my_pred->data[i][my_true->col*lv+j] - my_true->data[i][j]);
        sstot += square(my_true->data[i][j] - ymean->data[j]);
      }

      if(bias != NULL){
        double sum_yi = 0.f, sum_xi = 0.f;
        /*ypredaverage /= (double)my->row;*/
        for(i = 0; i < my_true->row; i++){
          sum_yi+=(my_pred->data[i][my_true->col*lv+j]*(my_true->data[i][j]-ymean->data[j]));
          sum_xi+=(my_true->data[i][j]*(my_true->data[i][j]-ymean->data[j]));
        }
        /*sum_yi/sum_xi = m */
        (*bias)->data[lv][j] = fabs(1 - sum_yi/sum_xi);
        /*k = Y-(X*b);*/
      }

      if(ccoeff != NULL)
        (*ccoeff)->data[lv][j] = 1.f - (ssreg/sstot);

      if(stdev != NULL)
        (*stdev)->data[lv][j] = sqrt(ssreg/(double)my_true->row);
    }

  }
  DelDVector(&ymean);
}

void PLSDiscriminantAnalysisStatistics(matrix *my_true, matrix *my_score, array **roc, matrix **roc_auc, array **precision_recall, matrix **precision_recall_ap)
{
  size_t nlv, lv, i, j, k, n_y;
  matrix *roc_;
  matrix *pr_;
  dvector *y_true;
  dvector *y_score;
  dvector *auc_row;
  dvector *ap_row;
  double auc, ap;

  nlv = (size_t) my_score->col/my_true->col;
  n_y = my_true->col;

  NewDVector(&y_true, my_true->row);
  NewDVector(&y_score, my_true->row);

  for(lv = 0; lv < nlv; lv++){
    if(roc != NULL)
      AddArrayMatrix(roc, my_true->row, my_true->col*2);

    if(precision_recall != NULL)
      AddArrayMatrix(precision_recall, my_true->row, my_true->col*2);


    initDVector(&auc_row);
    initDVector(&ap_row);

    k = 0;
    for(j = 0; j < my_true->col; j++){
      for(i = 0; i < my_true->row; i++){
        y_true->data[i] = my_true->data[i][j];
        y_score->data[i] = my_score->data[i][n_y*lv+j];
      }

      initMatrix(&roc_);
      ROC(y_true, y_score,  &roc_, &auc);
      initMatrix(&pr_);
      PrecisionRecall(y_true, y_score,  &pr_, &ap);
      DVectorAppend(&auc_row, auc);
      DVectorAppend(&ap_row, ap);

      for(i = 0; i < my_true->row; i++){
        if(roc != NULL){
          (*roc)->m[lv]->data[i][k] = roc_->data[i][0];
          (*roc)->m[lv]->data[i][k+1] = roc_->data[i][1];
        }

        if(precision_recall != NULL){
          (*precision_recall)->m[lv]->data[i][k] = pr_->data[i][0];
          (*precision_recall)->m[lv]->data[i][k+1] = pr_->data[i][1];
        }
      }
      k+=2;

      DelMatrix(&roc_);
      DelMatrix(&pr_);
    }

    if(roc_auc != NULL)
      MatrixAppendRow(roc_auc, auc_row);

    if(precision_recall_ap != NULL)
      MatrixAppendRow(precision_recall_ap, ap_row);

    DelDVector(&auc_row);
    DelDVector(&ap_row);
  }

  DelDVector(&y_true);
  DelDVector(&y_score);
}

/*
 * VIP are calculated according to the formula:
 * VIP[j][pc] = sqrt(n_predvars * Sum((b[pc]^2*t[pc]*t^[pc]) * (w[j][pc]/||w[pc]||)^2) / Sum((b[pc]^2*t[pc]*t^[pc])))
 */
void PLSVIP(PLSMODEL *model, matrix **vip)
{
  size_t nlv = model->xscores->col;
  size_t npred = model->xloadings->row;
  size_t i, j, k;
  ResizeMatrix(vip, npred, nlv);
  for(i = 0; i < nlv; i++){
    for(j = 0; j < i; j++){
      double n = 0.f;
      double d = 0.f;
      for(k = 0; k < npred; k++){
        /*double tmp = model->b[j]^2*t[pc]*t^[pc];
        n += */
        d = k;
        n += d;
      }
    }
  }

}

int GetLVCCutoff_(matrix *rq2y){
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

int GetLVCCutoff(matrix *coeff){
  size_t i, j, cutoff = 0;
  double prec;
  dvector *coeff_avg;
  NewDVector(&coeff_avg, coeff->row);
  for(i = 0; i < coeff->row; i++){
    for(j = 0; j < coeff->col; j++){
      coeff_avg->data[i] += coeff->data[i][j];
    }
  }

  prec = coeff_avg->data[0];
  for(i = 1; i < coeff_avg->size-1; i++){
    if(prec < coeff_avg->data[i] && coeff_avg->data[i] < coeff_avg->data[i+1]){
      prec = coeff_avg->data[i];
      cutoff = i+1;
      continue;
    }
    else{
      break;
    }
  }
  DelDVector(&coeff_avg);
  return cutoff;
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

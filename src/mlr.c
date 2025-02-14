/* mlr.c
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

#include "mlr.h"
#include "memwrapper.h"
#include "matrix.h"

#include "algebra.h"
#include "numeric.h"
#include "statistic.h"

#include <math.h>


void NewMLRModel(MLRMODEL** m)
{
  (*m) = xmalloc(sizeof(MLRMODEL));
  initMatrix(&(*m)->b);
  initMatrix(&(*m)->recalculated_y);
  initMatrix(&(*m)->predicted_y);
  initMatrix(&(*m)->recalc_residuals);
  initMatrix(&(*m)->pred_residuals);
  initDVector(&(*m)->ymean);
  initDVector(&(*m)->r2y_model);
  initDVector(&(*m)->q2y);
  initDVector(&(*m)->sdep);/* Standard Deviation over Prediction */
  initDVector(&(*m)->sdec); /* Standard Deviation over Recalculating */
  initDVector(&(*m)->bias);
  initMatrix(&(*m)->r2q2scrambling);
}

void DelMLRModel(MLRMODEL** m)
{
  DelMatrix(&(*m)->b);
  DelMatrix(&(*m)->recalculated_y);
  DelMatrix(&(*m)->predicted_y);
  DelMatrix(&(*m)->recalc_residuals);
  DelMatrix(&(*m)->pred_residuals);
  DelDVector(&(*m)->ymean);
  DelDVector(&(*m)->r2y_model);
  DelDVector(&(*m)->q2y);
  DelDVector(&(*m)->sdep);/* Standard Deviation over Prediction */
  DelDVector(&(*m)->sdec); /* Standard Deviation over Recalculating */
  DelDVector(&(*m)->bias);
  DelMatrix(&(*m)->r2q2scrambling);
  xfree((*m));
}

void MLR(matrix* mx, matrix* my, MLRMODEL* model, ssignal *s)
{
  size_t i, j;
  matrix *mx_;
  dvector *y_;
  dvector *b_;
  NewMatrix(&mx_, mx->row, mx->col+1);
  NewDVector(&y_, mx->row);

  /* Preparate the matrix to correlate with y trough OLS */
  MatrixCheck(mx);

  /*Adding constant coefcient for polinomial equation */
  for(j = 1; j < mx_->col; j++){
    for(i = 0; i < mx_->row; i++){
      mx_->data[i][j] = mx->data[i][j-1]; 
    }
  }

  for(i = 0; i < mx_->row; i++){
    mx_->data[i][0] = 1.f;
  }

  /*build models with more than one y*/
  for(j = 0; j < my->col; j++){
    for(i = 0; i < my->row; i++){
      y_->data[i] = my->data[i][j];
    }

    initDVector(&b_);

    OrdinaryLeastSquares(mx_, y_, b_);

    MatrixAppendCol(model->b, b_);
    DelDVector(&b_);
  }

  /* Recalculate y from model */
  MatrixColAverage(my, model->ymean);
  MLRPredictY(mx, my, model, model->recalculated_y, model->recalc_residuals, model->r2y_model, model->sdec);

  DelDVector(&y_);
  DelMatrix(&mx_);
}

/* Predict y using the equation:
 *
 * y = c + ax1 + bx2 + cx3 + dx4 + ...
 *
 */
void MLRPredictY(matrix* mx,
                 matrix *my,
                 MLRMODEL *model,
                 matrix *predicted_y,
                 matrix *predicted_residuals,
                 dvector *r2y,
                 dvector *sdep)
{
  if(mx->col == (model->b->row-1)){
    size_t i, j, k;
    double ypred, rss, tss;

    if(predicted_y->row != mx->row && predicted_y->col != model->b->col){
      ResizeMatrix(predicted_y, mx->row, model->b->col);
    }

    for(k = 0; k < model->b->col; k++){ /*FOR EACH Y*/
      for(i = 0; i < mx->row; i++){ /* FOR EACH OBJECT */
        ypred = model->b->data[0][k]; /* CONSTANT c + ...*/
        for(j = 0; j < mx->col; j++){
          ypred += mx->data[i][j]*model->b->data[j+1][k];
        }
        predicted_y->data[i][k] = ypred;
      }
    }

    if(my != 0){
      if(predicted_residuals != 0){
        ResizeMatrix(predicted_residuals, my->row, my->col);
      }
      /*Calculate Residual Error
      *      ei = y^i - yi.
      * where:
      *   yi = real y
      *   y^i = predicted y
      *
      * Calculate R^2 and SDEC
      */

      for(j = 0; j < my->col; j++){
        rss = tss = 0.f;
        for(i = 0; i < my->row; i++){
          if(predicted_residuals != 0){
            predicted_residuals->data[i][j] = predicted_y->data[i][j] - my->data[i][j];
          }
          rss += square(my->data[i][j] - predicted_y->data[i][j]);
          tss += square(my->data[i][j] - model->ymean->data[j]);
        }

        if(r2y != 0){
          DVectorAppend(r2y, 1-(rss/tss));
        }

        if(sdep != 0){
          DVectorAppend(sdep, sqrt(rss/my->row));
        }
      }
    }
  }
  else{
    fprintf(stderr, "Error!! Unable to compute MLR Prediction.\n The number of variable differ with the variables in model.\n");
  }
}

void MLRRegressionStatistics(matrix *my_true,
                             matrix *my_pred,
                             dvector *ccoeff,
                             dvector *rmse,
                             dvector *bias)
{
  size_t i, j;
  dvector *yt;
  dvector *yp;

  if(ccoeff != NULL)
    DVectorResize(ccoeff,my_true->col);

  if(rmse != NULL)
    DVectorResize(rmse, my_true->col);

  if(bias != NULL)
    DVectorResize(bias, my_true->col);

  for(j = 0; j < my_true->col; j++){
    initDVector(&yt);
    initDVector(&yp);
    for(i = 0; i < my_true->row; i++){
      DVectorAppend(yt, my_true->data[i][j]);
      DVectorAppend(yp, my_pred->data[i][j]);
    }

    if(bias != NULL){
      bias->data[j] = BIAS(yt, yp);
    }

    if(ccoeff != NULL){
      ccoeff->data[j] = R2(yt, yp);
    }

    if(rmse != NULL){
      rmse->data[j] = RMSE(yt, yp);
    }
    DelDVector(&yt);
    DelDVector(&yp);
  }
}

void PrintMLR(MLRMODEL *m)
{
  puts("b Coeffcicient");
  PrintMatrix(m->b);

  puts("R^2");
  PrintDVector(m->r2y_model);

  puts("SDEC");
  PrintDVector(m->sdec);

  puts("Recalculated y");
  PrintMatrix(m->recalculated_y);

  puts("Recalculated Residuals");
  PrintMatrix(m->recalc_residuals);

  if(m->q2y->size > 0){
    puts("Q^2");
    PrintDVector(m->q2y);

    puts("SDEP");
    PrintDVector(m->sdep);

    puts("BIAS");
    PrintDVector(m->bias);

    puts("Predicted y");
    PrintMatrix(m->predicted_y);

    puts("Predicted Residuals");
    PrintMatrix(m->pred_residuals);
  }

  if(m->r2q2scrambling->row > 0){
    puts("R2 Q^2 Y Scrambling");
    PrintMatrix(m->r2q2scrambling);
  }
}

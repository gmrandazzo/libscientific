/* statistic.c
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

#include "statistic.h"
#include "numeric.h"
#include "memwrapper.h"
#include <math.h>
#include <float.h>

/*
 * Calculate R^2
 */
double R2(dvector *ytrue, dvector *ypred)
{
  size_t i, ny;
  double ssreg, sstot, yavg;
  ssreg = sstot = yavg = 0.f;
  ny = 0;
  for(i = 0; i < ytrue->size; i++){
    if(FLOAT_EQ(ytrue->data[i], MISSING, 1e-1)){
      continue;
    }
    else{
      yavg += ytrue->data[i];
      ny+=1;
    }
  }
  yavg /= (double)ny;

  for(i = 0; i < ytrue->size; i++){
    if(FLOAT_EQ(ytrue->data[i], MISSING, 1e-1)){
      continue;
    }
    else{
      ssreg += square(ypred->data[i] - ytrue->data[i]);
      sstot += square(ytrue->data[i] - yavg);
    }
  }
  return 1.f - (ssreg/sstot);
}

/*
 * Calculate mean absolute error
 */
double MAE(dvector *ytrue, dvector *ypred)
{
  size_t i, ny;
  double mae;
  mae = 0.f;
  ny = 0;
  for(i = 0; i < ytrue->size; i++){
    if(FLOAT_EQ(ytrue->data[i], MISSING, 1e-1)){
      continue;
    }
    else{
      mae += fabs(ypred->data[i] - ytrue->data[i]);
      ny += 1;
    }
  }
  mae /= (double)ny;
  return mae;
}

/*
 * Calculate mean absolute error
 */
double MSE(dvector *ytrue, dvector *ypred)
{
  size_t i, ny;
  double mse;
  mse = 0.f;
  ny = 0;
  for(i = 0; i < ytrue->size; i++){
    if(FLOAT_EQ(ytrue->data[i], MISSING, 1e-1)){
      continue;
    }
    else{
      mse += square(ypred->data[i] - ytrue->data[i]);
      ny += 1;
    }
  }
  mse /= (double)ny;
  return mse;
}

/**
 * Blue (1978) ACM Transactions on Mathematical Software 4, 15-23.
 * Calculate MSE avoiding overflow and underflow
 */
double mse_blue(dvector *ytrue, dvector *ypred)
{
  size_t i, ny;
  double sum_small, sum_medium, sum_large;
  double b1, b2;
  double diff, abs_diff;

  sum_small = sum_medium = sum_large = 0.0;
  ny = 0;

  /* Constants for Blue's algorithm
   * b1: values smaller than this may underflow when squared
   * b2: values larger than this may overflow when squared
   */
  b1 = sqrt(DBL_MIN / DBL_EPSILON);
  b2 = sqrt(DBL_MAX * DBL_EPSILON);

  for(i = 0; i < ytrue->size; i++){
    if(FLOAT_EQ(ytrue->data[i], MISSING, 1e-1)){
      continue;
    }
    else{
      diff = ypred->data[i] - ytrue->data[i];
      abs_diff = fabs(diff);
      ny += 1;

      if (abs_diff < b1) {
        if (abs_diff > 0) {
            sum_small += square(diff / b1);
        }
      } else if (abs_diff > b2) {
        sum_large += square(diff / b2);
      } else {
        sum_medium += square(diff);
      }
    }
  }

  if (ny == 0) return 0.0;

  /* Combine the sums robustly */
  if (sum_large > 0) {
    return square(b2) * (sum_large / (double)ny + ((sum_medium + square(b1) * sum_small) / square(b2)) / (double)ny);
  } else if (sum_medium > 0) {
    return (sum_medium / (double)ny + (square(b1) * sum_small / (double)ny));
  } else {
    return square(b1) * (sum_small / (double)ny);
  }
}

double RMSE(dvector *ytrue, dvector *ypred)
{
  return sqrt(MSE(ytrue, ypred));
}

/*
 * bias = 1 - slope(ypred,ytrue)
 */
double BIAS(dvector *ytrue, dvector *ypred)
{
  size_t i, ny;
  double sum_yi, sum_xi;
  double yavg;
  sum_yi = sum_xi = yavg = 0.f;
  ny = 0;
  for(i = 0; i < ytrue->size; i++){
    if(FLOAT_EQ(ytrue->data[i], MISSING, 1e-1)){
      continue;
    }
    else{
      yavg += ytrue->data[i];
      ny+=1;
    }
  }
  yavg /= (double)ny;

  /*ypredaverage /= (double)my->row;*/
  for(i = 0; i < ytrue->size; i++){
    if(FLOAT_EQ(ytrue->data[i], MISSING, 1e-1)){
      continue;
    }
    else{
      sum_yi+=(ypred->data[i]*(ytrue->data[i]-yavg));
      sum_xi+=(ytrue->data[i]*(ytrue->data[i]-yavg));
    }
  }
  /*sum_yi/sum_xi = m */
  return fabs(1 - sum_yi/sum_xi);
}

/*k = Y-(X*b);*/
/*
 * TP = True PositivePredictedValue
 * FN = False Negative
 *
 * Sensitivity = TP / TP + FN
 */
void Sensitivity(dvector *dtp,
                 double thmin,
                 double thmax,
                 double thstep,
                 matrix *s)
{
  size_t i, row, tp, fn;
  double dx;

  ResizeMatrix(s, (size_t) ceil((thmax-thmin)/thstep), 2); /* two column: x and y */

  row = 0;
  dx = thstep;

  for(row = 0; row < s->row; row++){
    tp = fn = 0;

    for(i = 0; i < dtp->size; i++){
      if(FLOAT_EQ(dtp->data[i], MISSING, 1e-1)){
        continue;
      }
      else{
        if(dtp->data[i] < dx)
          tp++;
        else if(FLOAT_EQ(dtp->data[i], dx, EPSILON))
          tp++;
        else
          fn++;
      }
    }

    if(tp == 0){
      s->data[row][0] = dx;
      s->data[row][1] = 0.f;
    }
    else{
      s->data[row][0] = dx;
      s->data[row][1] = (double)tp/(double)(tp+fn);
    }

    dx+=thstep;
  }

}

/*
 * TP = True Positive
 * FP = False Positive
 *
 * PPV = TP / TP + FP
 */
void PositivePredictedValue(dvector* dtp,
                            dvector* dtn,
                            double thmin,
                            double thmax,
                            double thstep,
                            matrix *p)
{
  size_t i, row, tp, fp;
  double dx;

  ResizeMatrix(p, (size_t) ceil((thmax-thmin)/thstep), 2); /* two column: x and y */

  row = 0;
  dx = thstep;

  for(row = 0; row < p->row; row++){
    tp = fp = 0;

    for(i = 0; i < dtp->size; i++){
      if(FLOAT_EQ(dtp->data[i], MISSING, 1e-1)){
        continue;
      }
      else{
        if(dtp->data[i] < dx)
          tp++;
        else if(FLOAT_EQ(dtp->data[i], dx, EPSILON))
          tp++;
      }
    }

    for(i = 0; i < dtn->size; i++){
      if(FLOAT_EQ(dtn->data[i], MISSING, 1e-1)){
        continue;
      }
      else{
        if(getDVectorValue(dtn, i) < dx)
          fp++;
        else if(FLOAT_EQ(getDVectorValue(dtn, i), dx, EPSILON))
          fp++;
      }
    }

    if(tp == 0){
      p->data[row][0] = dx;
      p->data[row][1] = 0;
    }
    else{
      p->data[row][0] = dx;
      p->data[row][1] = (double)tp/(double)(tp+fp);
    }

    dx+=thstep;
  }
}


/*
 * This function Code the matrix by using the 0, 1 scale
 * according to the function
 *
 * x = (val - mid) / step
 */
void MatrixCode(matrix* inmx, matrix* outmx)
{
  size_t i, j;
  double min, nmin, max, nmax, mid, step, tmp;
  ResizeMatrix(outmx, inmx->row, inmx->col);

  for(j = 0; j < inmx->col; j++){
    min = max = nmin = nmax = getMatrixValue(inmx, 0, j);
    for(i = 0; i < inmx->row; i++){
      tmp = getMatrixValue(inmx, i, j);
      if(FLOAT_EQ(tmp, MISSING, 1e-1)){
        continue;
      }
      else{
        if(tmp < min){
          min = tmp;
        }
        if(tmp > max){
          max = getMatrixValue(inmx, i, j);
        }
      }
    }

    mid = min + ((max - min) / 2);

    for(i = 0; i < inmx->row; i++){
      tmp = getMatrixValue(inmx, i, j);
      if(FLOAT_EQ(tmp, MISSING, 1e-1)){
        continue;
      }
      else{
        if(tmp > min && tmp < mid && tmp < nmin){
          nmin = tmp;
        }

        if(tmp < max && tmp > mid && tmp > nmax){
          nmax = getMatrixValue(inmx, i, j);
        }
      }
    }

    if(nmax != nmin){
      step = (nmax- nmin) / 2;
      mid = nmin + step;
    }
    else{
      step = (max - min) / 2;
      mid = min + step;
    }

    for(i = 0; i < inmx->row; i++){
      tmp = getMatrixValue(inmx, i, j);
      if(FLOAT_EQ(tmp, MISSING, 1e-1)){
        continue;
      }
      else{
        setMatrixValue(outmx, i, j, (tmp - mid) / step);
      }
    }
  }
}

/*
 * This function expand the matrix to Bifactorial cross product
 * a, b, c
 * a*a, b*b, c*c,
 * a*b, a*c, b*c
 */
void BifactorialMatrixExpansion(matrix* inmx, matrix* outmx)
{
  size_t i, j, k;
  size_t n = inmx->col, n_k, comb = 1;

  if(inmx->col-2 != 0){
    /*n! / (n-k)!*k!*/
    n = Factorial(inmx->col);
    n_k = Factorial(inmx->col-2);
    comb = (n / (n_k *2));
  }

  ResizeMatrix(outmx, inmx->row, 2*inmx->col + comb + 1);


  /* Generate the matrix */
  for(i = 0; i < inmx->row; i++){
    setMatrixValue(outmx, i, 0, 1);
  }

  for(i = 0; i < inmx->row; i++){
    for(j = 0; j < inmx->col; j++){
      setMatrixValue(outmx, i, j+1, getMatrixValue(inmx, i, j));
      setMatrixValue(outmx, i, inmx->col+1+j+comb, square(getMatrixValue(inmx, i, j)));
    }
  }

  for(i = 0; i < inmx->row; i++){
    for(j = 0; j < inmx->col; j++){
      for(k = j+1; k < inmx->col; k++){
        setMatrixValue(outmx, i, inmx->col+j+k, getMatrixValue(inmx, i, j)*getMatrixValue(inmx, i, k));
      }
    }
  }

}

/* Algorithm from:
 * An introduction to ROC analysis
 * Tom Fawcett
 * Pattern Recognition Letters 27 (2006) 861–874
 * doi: 10.1016/j.patrec.2005.10.010
 *
 * Calculate the roc curve to plot and it's AUC
 */
void ROC(dvector *y_true, dvector *y_score,  matrix *roc, double *auc)
{
  size_t i;
  matrix *yy;
  initMatrix(&yy);
  MatrixAppendCol(yy, y_true);
  MatrixAppendCol(yy, y_score);
  MatrixReverseSort(yy, 1); /* sort by y_score*/

  dvector *roc_row;
  NewDVector(&roc_row, 2);
  DVectorSet(roc_row, 0.0);
  MatrixAppendRow(roc, roc_row);
  /*Calculate the number of tp and tn*/
  size_t n_tp = 0;
  size_t n_tn = 0;
  for(i = 0; i < yy->row; i++){
    if(FLOAT_EQ(yy->data[i][0], MISSING, 1e-1)){
      continue;
    }
    else{
      if(FLOAT_EQ(yy->data[i][0], 1.0, 1e-1) == 1){
        n_tp += 1;
      }
      else
        n_tn += 1;
    }
  }

  size_t tp = 0;
  size_t fp = 0;
  for(i = 0; i < yy->row; i++){
    if(FLOAT_EQ(yy->data[i][0], MISSING, 1e-1)){
      continue;
    }
    else{
      if(FLOAT_EQ(yy->data[i][0], 1.0, 1e-1) == 1){
        tp += 1;
      }
      else{
        fp += 1;
      }
      roc_row->data[0] = (double)fp/(double)n_tn;
      roc_row->data[1] = (double)tp/(double)n_tp;
      MatrixAppendRow(roc, roc_row);
    }
  }
  DelDVector(&roc_row);
  DelMatrix(&yy);
  (*auc) = curve_area(roc, 0);
  if(_isnan_((*auc)))
    (*auc) = 0.f;
}

/* Algorithm from:
 * An introduction to ROC analysis
 * Tom Fawcett
 * Pattern Recognition Letters 27 (2006) 861–874
 * doi: 10.1016/j.patrec.2005.10.010
 *
 * Calculate the precision-recall curve to plot and it's average precision
 */
void PrecisionRecall(dvector *y_true,
                     dvector *y_score,
                     matrix *pr,
                     double *ap)
{
  size_t i;
  matrix *yy;
  initMatrix(&yy);
  MatrixAppendCol(yy, y_true);
  MatrixAppendCol(yy, y_score);
  MatrixReverseSort(yy, 1); /* sort by y_score*/

  dvector *pr_row;
  NewDVector(&pr_row, 2);
  pr_row->data[0] = 0.f; // recall
  pr_row->data[1] = 1.f; // precision
  MatrixAppendRow(pr, pr_row);
  /*Calculate the number of tp and tn*/
  size_t n_tp = 0;

  for(i = 0; i < yy->row; i++){
    if(FLOAT_EQ(yy->data[i][0], MISSING, 1e-1)){
      continue;
    }
    else{
      if(FLOAT_EQ(yy->data[i][0], 1.0, 1e-1) == 1){
        n_tp += 1;
      }
    }
  }

  size_t tp = 0;
  size_t fp = 0;
  for(i = 0; i < yy->row; i++){
    if(FLOAT_EQ(yy->data[i][0], MISSING, 1e-1)){
      continue;
    }
    else{
      if(FLOAT_EQ(yy->data[i][0], 1.0, 1e-1) == 1){
        tp += 1;
      }
      else{
        fp += 1;
      }
      pr_row->data[0] = (double)tp/(double)n_tp;
      pr_row->data[1] = (double)tp/(double)(tp+fp);
      MatrixAppendRow(pr, pr_row);
    }
  }
  DelDVector(&pr_row);
  DelMatrix(&yy);
  (*ap) = curve_area(pr, 0);
  if(_isnan_((*ap)))
    (*ap) = 0.f;
}

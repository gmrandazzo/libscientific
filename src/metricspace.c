/* metricspace.c
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

#include "metricspace.h"
#include "numeric.h"
#include <math.h>


void EuclideanDistance(matrix* m1, matrix* m2, matrix **distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double dist;

    ResizeMatrix(distances, m2->row, m1->row); /* each column is a distance that correspond to m1->row */

    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        dist = 0.f;
        for(j = 0; j < m1->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          dist += square(m1->data[i][j] - m2->data[k][j]);
        }
        (*distances)->data[k][i] = sqrt(dist);
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Euclidean Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void SquaredEuclideanDistance(matrix *m1, matrix *m2, matrix **distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double dist;

    ResizeMatrix(distances, m2->row, m1->row); /* each column is a distance that correspond to m1->row */

    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        dist = 0.f;
        for(j = 0; j < m1->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          dist += square(m1->data[i][j] - m2->data[k][j]);
        }
        (*distances)->data[k][i] = dist;
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Euclidean Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void ManhattanDistance(matrix* m1, matrix* m2, matrix** distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double dist;

    ResizeMatrix(distances, m2->row, m1->row); /* each column is a distance that correspond to m1->row */

    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        dist = 0.f;
        for(j = 0; j < m1->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          dist += fabs(m1->data[i][j] - m2->data[k][j]);
        }
        (*distances)->data[k][i] = dist;
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Manhattan Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void CosineDistance(matrix* m1, matrix* m2, matrix** distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double n, d_a, d_b;

    ResizeMatrix(distances, m2->row, m1->row); /* each column is a distance that correspond to m1->row */

    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        n = 0.f; d_a = 0.f; d_b = 0.f;
        for(j = 0; j < m1->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          n += m1->data[i][j] * m2->data[k][j];
          d_a += square(m1->data[i][j]);
          d_b += square(m2->data[k][j]);
        }
        (*distances)->data[k][i] = n/(sqrt(d_a)*sqrt(d_b));
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Cosine Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

double MahalanobisDistance(matrix* g1, matrix* g2)
{
  if(g1->col == g2->col){
    size_t i, j;
    double dist, a, b;
    dvector *mean1, *mean2, *xdiff, *x;
    matrix *cov1, *cov2, *invcov;

    initDVector(&mean1);
    MatrixColAverage(g1, &mean1);

    initDVector(&mean2);
    MatrixColAverage(g2, &mean2);

    initMatrix(&cov1);
    MatrixCovariance(g1, &cov1);

    initMatrix(&cov2);
    MatrixCovariance(g2, &cov2);

    /* pooled covariance matrix in cov1*/
    for(i = 0; i < cov1->row; i++){
      for(j = 0; j < cov1->col; j++){
        a = (double)g1->row/(double)(g1->row+g2->row);
        b = (double)g2->row/(double)(g1->row+g2->row);
        cov1->data[i][j] = (cov1->data[i][j] * a) + (cov2->data[i][j] * b);
      }
    }

    initMatrix(&invcov);

    MatrixInversion(cov1, &invcov);

    NewDVector(&xdiff, mean1->size);

    for(i = 0; i < mean1->size; i++){
      xdiff->data[i] = mean1->data[i] - mean2->data[i];
    }

    NewDVector(&x, mean1->size);
    DVectorMatrixDotProduct(invcov, xdiff, x);

    dist = sqrt(DVectorDVectorDotProd(xdiff, x));

    DelMatrix(&invcov);
    DelDVector(&mean2);
    DelDVector(&mean1);
    DelDVector(&xdiff);
    DelMatrix(&cov2);
    DelMatrix(&cov1);
    DelDVector(&x);

    return dist;
  }
  else{
    fprintf(stderr, "Unable to compute Euclidean Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
    return 0.f;
  }
}

/* algorithm taken from:
 * Natural Cubic Spline: Numerical Analysis book Richard L. Burden and J. Douglas Faires
 * pag. 149
 * ISBN-13: 978-0-538-73351-9
 *
 * N.B.: If two points have different y but share same x the algorithm will fail
 */

void cubic_spline_interpolation(matrix *xy, matrix **S)
{
  size_t i;
  int j;
  size_t np1 = xy->row;
  size_t n = np1-1;
  dvector *a, *b, *d, *h, *alpha, *c, *l, *u, *z;
  NewDVector(&a, np1);
  for(i = 0; i < np1; i++)
    a->data[i] = xy->data[i][1];

  NewDVector(&b, n);
  NewDVector(&d, n);
  NewDVector(&h, n);
  NewDVector(&alpha, n);
  for(i = 0; i < n; i++){
    h->data[i] = xy->data[i+1][0] - xy->data[i][0];
  }

  for(i = 1; i < n; i++){
    alpha->data[i] = 3.f/h->data[i]*(a->data[i+1]-a->data[i]) - 3.f/h->data[i-1]*(a->data[i]-a->data[i-1]);
  }

  /*PrintDVector(h);
  PrintDVector(alpha);*/

  NewDVector(&c, np1);
  NewDVector(&l, np1);
  NewDVector(&u, np1);
  NewDVector(&z, np1);

    l->data[0] = 1.f;
    u->data[0] = 0.f;
    z->data[0] = 0.f;
    for(i = 1; i < n; i++){
      l->data[i] = 2*(xy->data[i+1][0]-xy->data[i-1][0]) - h->data[i-1]*u->data[i-1];
      u->data[i] = h->data[i]/l->data[i];
      z->data[i] = (alpha->data[i]-h->data[i-1]*z->data[i-1])/l->data[i];
    }

    l->data[n] = 1.f;
    z->data[n] = 0.f;
    c->data[n] = 0.f;

    for(j = n-1; j > -1; j--){
        c->data[j] = z->data[j] - u->data[j]*c->data[j+1];
        b->data[j] = (a->data[j+1]-a->data[j])/h->data[j] - (h->data[j]*(c->data[j+1]+2*c->data[j]))/3.f;
        d->data[j] = (c->data[j+1]-c->data[j])/(3.f*h->data[j]);
    }

    /* Write SPLINE EQUATIONS
     * first column is the x range
     * second is a, then b,c and d.
     */
    ResizeMatrix(S, n, 5);
    for(i = 0; i < n; i++){
      (*S)->data[i][0] = xy->data[i][0];
      (*S)->data[i][1] = a->data[i];
      (*S)->data[i][2] = b->data[i];
      (*S)->data[i][3] = c->data[i];
      (*S)->data[i][4] = d->data[i];
    }

    DelDVector(&a);
    DelDVector(&b);
    DelDVector(&d);
    DelDVector(&h);
    DelDVector(&alpha);
    DelDVector(&c);
    DelDVector(&l);
    DelDVector(&u);
    DelDVector(&z);
}

void cubic_spline_predict(dvector *x_, matrix *S, dvector **y_pred)
{
  size_t i, j, n;
  double x, xi, y;
  n = S->row-1;
  DVectorResize(y_pred, x_->size);
  /* Now interpolate using the equations:
   * Sj(x) = aj + bj(x − xj) + cj(x − xj)^2 + dj(x − xj)^3
   * for xj ≤ x ≤ xj+1)
   */
  for(i = 0; i < x_->size; i++){
    x = x_->data[i];
    for(j = 0; j < n; j++){
      xi = S->data[j][0];
      y = 0.f;
      if((x > xi || FLOAT_EQ(x, xi, 1e-2)) &&
      (x < S->data[j+1][0] || FLOAT_EQ(x, S->data[j+1][0], 1e-2))){
        //printf("for %f selecting %d\n", x, j);
        y = S->data[j][1] + S->data[j][2]*(x-xi) + S->data[j][3]*(x-xi)*(x-xi) + S->data[j][4]*(x-xi)*(x-xi)*(x-xi);
        break;
      }
      else{
        continue;
      }
    }
    (*y_pred)->data[i] = y;
  }
}

void interpolate(matrix *xy, size_t npoints, matrix **interp_xy)
{
  size_t i;
  double x, dx, xmin, xmax;
  matrix *S;
  dvector *x_interp, *y_pred;
  initMatrix(&S);
  cubic_spline_interpolation(xy, &S);

  ResizeMatrix(interp_xy, npoints, 2);
  MatrixColumnMinMax(xy, 0, &xmin, &xmax);
  dx = (xmax-xmin)/(double)(npoints-1);
  x = xmin;
  NewDVector(&x_interp, npoints);
  for(i = 0; i < npoints; i++){
    x_interp->data[i] = x;
    x+=dx;
  }

  initDVector(&y_pred);
  cubic_spline_predict(x_interp, S, &y_pred);
  for(i = 0; i < x_interp->size; i++){
    (*interp_xy)->data[i][0] = x_interp->data[i];
    (*interp_xy)->data[i][1] = y_pred->data[i];
  }

  DelDVector(&x_interp);
  DelDVector(&y_pred);
  DelMatrix(&S);
}


/*
 * area of the curve via the  trapezoid rule
 */
double curve_area(matrix *xy, size_t intervals)
{
  size_t i;
  matrix *interp_xy;
  initMatrix(&interp_xy);
  /*If intervals > 0 interpolate with natural cubic splines to have more fine area */
  if(intervals > 0){
    /* If two points have different y but share same x the algorithm will fail*/
    interpolate(xy, intervals, &interp_xy);
  }
  else{
    MatrixCopy(xy, &interp_xy);
    intervals = xy->row;
  }

  double base, height;
  double area = 0.f;
  for(i = 0; i < intervals-1; i++){
      /* Trapezoidal method */
      base = interp_xy->data[i+1][0]-interp_xy->data[i][0];
      height = ((interp_xy->data[i][1]+interp_xy->data[i+1][1])/2.f);
      area += base*height;
  }
  DelMatrix(&interp_xy);
  return area;
}

/* Algorithm from:
 * An introduction to ROC analysis
 * Tom Fawcett
 * Pattern Recognition Letters 27 (2006) 861–874
 * doi: 10.1016/j.patrec.2005.10.010
 *
 * Calculate the roc curve to plot and it's AUC
 */
void ROC(dvector *y_true, dvector *y_score,  matrix **roc, double *auc)
{
  size_t i;
  matrix *yy;
  initMatrix(&yy);
  MatrixAppendCol(&yy, y_true);
  MatrixAppendCol(&yy, y_score);
  MatrixReverseSort(yy, 1); /* sort by y_score*/

  dvector *roc_row;
  NewDVector(&roc_row, 2);
  DVectorSet(roc_row, 0.0);
  MatrixAppendRow(roc, roc_row);
  /*Calculate the number of tp and tn*/
  size_t n_tp = 0;
  size_t n_tn = 0;
  for(i = 0; i < yy->row; i++){
    if(FLOAT_EQ(yy->data[i][0], 1.0, 1e-1) == 1){
      n_tp += 1;
    }
    else
      n_tn += 1;
  }

  size_t tp = 0;
  size_t fp = 0;
  for(i = 0; i < yy->row; i++){
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
  DelDVector(&roc_row);
  DelMatrix(&yy);
  (*auc) = curve_area((*roc), 0);
}

/* Algorithm from:
 * An introduction to ROC analysis
 * Tom Fawcett
 * Pattern Recognition Letters 27 (2006) 861–874
 * doi: 10.1016/j.patrec.2005.10.010
 *
 * Calculate the precision-recall curve to plot and it's average precision
 */
void PrecisionRecall(dvector *y_true, dvector *y_score,  matrix **pr, double *ap)
{
  size_t i;
  matrix *yy;
  initMatrix(&yy);
  MatrixAppendCol(&yy, y_true);
  MatrixAppendCol(&yy, y_score);
  MatrixReverseSort(yy, 1); /* sort by y_score*/

  dvector *pr_row;
  NewDVector(&pr_row, 2);
  pr_row->data[0] = 0.f; // recall
  pr_row->data[1] = 1.f; // precision
  MatrixAppendRow(pr, pr_row);
  /*Calculate the number of tp and tn*/
  size_t n_tp = 0;
  size_t n_tn = 0;
  for(i = 0; i < yy->row; i++){
    if(FLOAT_EQ(yy->data[i][0], 1.0, 1e-1) == 1){
      n_tp += 1;
    }
    else
      n_tn += 1;
  }

  size_t tp = 0;
  size_t fp = 0;
  for(i = 0; i < yy->row; i++){
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
  DelDVector(&pr_row);
  DelMatrix(&yy);
  (*ap) = curve_area((*pr), 0);
}

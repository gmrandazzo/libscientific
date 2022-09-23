#include "matrix.h"
#include "numeric.h"
#include "interpolate.h"


/* algorithm taken from:
 * Natural Cubic Spline: Numerical Analysis book Richard L. Burden and J. Douglas Faires
 * pag. 149
 * ISBN-13: 978-0-538-73351-9
 *
 * N.B.: If two points have different y but share same x the algorithm will fail
 */

void cubic_spline_interpolation(matrix *xy, matrix *S)
{
  size_t i;
  int j;
  size_t np1 = xy->row;
  size_t n = np1-1;
  double a[np1], b[n], d[n], h[n], alpha[n], c[np1], l[np1], u[np1], z[np1];

  for(i = 0; i < np1; i++)
    a[i] = xy->data[i][1];

  for(i = 0; i < n; i++){
    h[i] = xy->data[i+1][0] - xy->data[i][0];
  }

  for(i = 1; i < n; i++){
    alpha[i] = 3.f/h[i]*(a[i+1]-a[i]) - 3.f/h[i-1]*(a[i]-a[i-1]);
  }

  l[0] = 1.f;
  u[0] = 0.f;
  z[0] = 0.f;
  for(i = 1; i < n; i++){
    l[i] = 2*(xy->data[i+1][0]-xy->data[i-1][0]) - h[i-1]*u[i-1];
    u[i] = h[i]/l[i];
    z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
  }

  l[n] = 1.f;
  z[n] = 0.f;
  c[n] = 0.f;

  for(j = n-1; j > -1; j--){
    c[j] = z[j] - u[j]*c[j+1];
    b[j] = (a[j+1]-a[j])/h[j] - (h[j]*(c[j+1]+2*c[j]))/3.f;
    d[j] = (c[j+1]-c[j])/(3.f*h[j]);
  }

    /* Write SPLINE EQUATIONS
     * first column is the x range
     * second is a, then b,c and d.
     */
    ResizeMatrix(S, n, 5);
    for(i = 0; i < n; i++){
      S->data[i][0] = xy->data[i][0];
      S->data[i][1] = a[i];
      S->data[i][2] = b[i];
      S->data[i][3] = c[i];
      S->data[i][4] = d[i];
    }
}

void cubic_spline_predict(dvector *x_, matrix *S, dvector *y_pred)
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
    y = MISSING;
    for(j = 0; j < n; j++){
      xi = S->data[j][0];
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

    /* if y is still missing,
     * this means that we are extrapolating over y.
     * Then we will use the last equation
     * with the last row of coefficients in S
     */
    if(FLOAT_EQ(y, MISSING, 1e-2)){
      xi = S->data[n][0];
      j = S->row-1;
      y = S->data[j][1] + S->data[j][2]*(x-xi) + S->data[j][3]*(x-xi)*(x-xi) + S->data[j][4]*(x-xi)*(x-xi)*(x-xi);
    }

    y_pred->data[i] = y;
  }
}

void interpolate(matrix *xy, size_t npoints, matrix *interp_xy)
{
  size_t i;
  double x, dx, xmin, xmax;
  matrix *S;
  dvector *x_interp, *y_pred;
  initMatrix(&S);
  cubic_spline_interpolation(xy, S);

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
  cubic_spline_predict(x_interp, S, y_pred);
  for(i = 0; i < x_interp->size; i++){
    interp_xy->data[i][0] = x_interp->data[i];
    interp_xy->data[i][1] = y_pred->data[i];
  }

  DelDVector(&x_interp);
  DelDVector(&y_pred);
  DelMatrix(&S);
}

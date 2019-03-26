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
    y = 0.f;
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

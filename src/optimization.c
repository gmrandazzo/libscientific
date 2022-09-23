/* optimization.c
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

#include "matrix.h"
#include "vector.h"
#include "optimization.h"
#include "memwrapper.h"
#include <stdio.h>
#include <math.h>
#include <time.h>

/*
 *  NelderMeadSimplex:
 *  1) Perform k+1 steps with k the number of variable to optimise
 *  2) Rank the response from the best to the worst which is the last one x_wost
 *  3) Establish the new condition according to the function:
 *      x_new =  c - c - x_wost
 *      where c is the centroid calculate with all the best less the worst one.
 */

void gen_centroids(matrix *x, dvector *c)
{
  /* x is already ordered
  * the last is skipped;
   */
  size_t i, j;
  for(j = 0; j < x->col-1; j++){
    c->data[j] = 0.f;
    for(i = 0; i < x->row-1; i++){
      c->data[j] += x->data[i][j];
    }
    c->data[j] /= (double)(x->row-1.0);
  }
}

void replace_xnp1(matrix *x, dvector *x_new, double x_new_res)
{
  size_t i;
  for(i = 0; i < x_new->size; i++){
    x->data[x->row-1][i] = x_new->data[i];
  }
  x->data[x->row-1][x->col-1] = x_new_res;
}

void shrink(matrix *x, double delta)
{
  size_t i, j;
  for(i = 1; i < x->row; i++){
    for(j = 0; j < x->col; j++){
      x->data[i][j] = x->data[0][j] + delta*(x->data[i][j]-x->data[0][j]);
    }
  }
}

/*
 * Implementing the Nelder-Mead simplex algorithm with adaptive parameters
 * Fuchang Gao, Lixing Han
 * Comput Optim Appl
 * DOI 10.1007/s10589-010-9329-3
 */
double NelderMeadSimplex(double (*func)(),
                         dvector *x0,
                         dvector *step,
                         double xtol,
                         size_t iter,
                         dvector *best)
{
  size_t i, j, iter_;
  matrix *x;

  dvector *c, *x_r, *x_e, *x_oc, *x_ic, *x_i, *x_tmp;
  double res, x_r_res, x_e_res, x_oc_res, x_ic_res;

  double alpha = 1; /* alpha >0*/
  double beta = 1+(2.f/(double)x0->size); /* beta > 1 */
  double gamma = 0.75-1/(2.f*(double)x0->size); /* 0 < gamma < 1 */
  double delta = 1- (1.f/(double)x0->size); /* 0 < delta < 1 */

  NewDVector(&c, x0->size);
  NewDVector(&x_tmp, x0->size);
  NewDVector(&x_r, x0->size);
  NewDVector(&x_e, x0->size);
  NewDVector(&x_oc, x0->size);
  NewDVector(&x_ic, x0->size);
  NewDVector(&x_i, x0->size);

  /*1) create k+1 steps with  k+1 value in wich the last column is the fval result! */
  NewMatrix(&x, x0->size+1, x0->size+1);
  if(step != NULL){
    for(i = 0; i < x0->size+1; i++){
      for(j = 0; j < x0->size; j++){
        if(i-1 == j){
          x->data[i][j] = x0->data[j] + step->data[j];
        }
        else{
          x->data[i][j] = x0->data[j];
        }
      }
      x->data[i][x->col-1] = MISSING; /* value of function */
    }
  }
  else{
    for(i = 0; i < x0->size+1; i++){
      for(j = 0; j < x0->size; j++){
        if(i-1 == j){
          x->data[i][j] = x0->data[j] + 0.5;
        }
        else{
          x->data[i][j] = x0->data[j];
        }
      }
      x->data[i][x->col-1] = MISSING; /* value of function */
    }
  }

  /* 2 perform k+1 calculation to initialize the output function */
  for(i = 0; i < x->row; i++){
    for(j = 0; j < x->col-1; j++){
      x_tmp->data[j] = x->data[i][j];
    }

    x->data[i][ x->col-1] = func(x_tmp);
  }

  /*puts("Initial matrix");
  PrintMatrix(x);*/

  /*3 Rank the initialized output from the best to the worst */
  MatrixSort(x, x->col-1);

  iter_ = 0;
  while(iter_ < iter){

    /*4 Establish a new condition (excluding the worst) for the next experiment and start with the loop */
    /* the last is the worst one!! which is replaced with the new conditions */
    gen_centroids(x, c);

    /* compute the reflection */
    for(i = 0; i < c->size; i++){
      x_r->data[i] = c->data[i]+alpha*(c->data[i]-x->data[x->row-1][i]);
      x_tmp->data[i] = x_r->data[i];
    }

    x_r_res = func(x_tmp);
    //printf("Compute reflection: %f\n", x_r_res);
    /* if f(0) < f(r) < f(n) */
    if(x->data[0][x->col-1] < x_r_res && x_r_res < x->data[x->row-2][x->col-1]){
      /*replace xn+1 with the reflection x_r*/
      replace_xnp1(x, x_r, x_r_res);
    }
    /*else if f(r) < f(0) */
    else if(x_r_res < x->data[0][x->col-1]){
      /* compute the expansion */
      for(i = 0; i < c->size; i++){
        x_e->data[i] = c->data[i]+beta*(x_r->data[i]-c->data[i]);
        x_tmp->data[i] = x_e->data[i];
      }
      x_e_res = func(x_tmp);
      //printf("Compute expansion: %f\n", x_e_res);
      /* if f(e) < f(r) */
      if(x_e_res < x_r_res){
        /*replace xn+1 with x_e otherwise replace xn+1 with xr*/
        replace_xnp1(x, x_e, x_e_res);
      }
      else{
        replace_xnp1(x, x_r, x_r_res);
      }
    }
    /* else if f(n) <= f(r) < f(n+1) */
    else if(x->data[x->row-2][x->col-1] <= x_r_res && x_r_res < x->data[x->row-1][x->col-1]){
      for(i = 0; i < c->size; i++){
        x_oc->data[i] = c->data[i]+gamma*(x_r->data[i]-c->data[i]);
        x_tmp->data[i] = x_oc->data[i];
      }
      x_oc_res = func(x_tmp);
      //printf("Compute outside contraction: %f\n", x_oc_res);
      /*if f(oc) <= f(r) */
      if(x_oc_res <= x_r_res){
        replace_xnp1(x, x_oc, x_oc_res);
      }
      else{
        /*shrink*/
        shrink(x, delta);

        for(i = 0; i < x->row; i++){
          for(j = 0; j < x->col-1; j++){
            x_tmp->data[j] = x->data[i][j];
          }
          x->data[i][ x->col-1] = func(x_tmp);
        }

      }
    }
    else if(x_r_res >= x->data[x->row-1][x->col-1]){
      for(i = 0; i < c->size; i++){
        x_ic->data[i] = c->data[i] - gamma*(x_r->data[i]-c->data[i]);
        x_tmp->data[i] = x_ic->data[i];
      }
      x_ic_res = func(x_tmp);
      //printf("Compute inside contraction: %f\n", x_ic_res);
      /* if f(ic) < f(n+1) */
      if(x_ic_res < x->data[x->row-1][x->col-1]){
        replace_xnp1(x, x_ic, x_ic_res);
      }
      else{
        shrink(x, delta);
        for(i = 0; i < x->row; i++){
          for(j = 0; j < x->col-1; j++){
            x_tmp->data[j] = x->data[i][j];
          }
          x->data[i][ x->col-1] = func(x_tmp);
        }
      }
    }
    MatrixSort(x, x->col-1);
    /*
    puts("Unsorted matrix");
    PrintMatrix(x);
    MatrixSort(x, x->col-1);

    puts("Sorted Matrix");
    PrintMatrix(x);
    puts("Centroid no WORST");
    PrintDVector(c);
    puts("Worst function value");
    printf("%f\n", x->data[x->row-1][x->col-1]);
    printf("------------------\n");
    sleep(2);
    */

    if(fabs(x->data[x->row-1][x->col-1]-x->data[0][x->col-1]) < xtol){
      break;
    }
    else{
      iter_ +=1;
    }
  }

  MatrixSort(x, x->col-1);

  /*puts("Final matrix");
  PrintMatrix(x);*/

  DVectorResize(best, x0->size);
  for(i = 0; i < x0->size; i++)
    best->data[i] = x->data[0][i];
  res = x->data[0][x0->size];

  DelDVector(&x_r);
  DelDVector(&x_e);
  DelDVector(&x_oc);
  DelDVector(&x_ic);
  DelDVector(&x_i);
  DelDVector(&c);
  DelDVector(&x_tmp);
  DelMatrix(&x);
  return res;
}

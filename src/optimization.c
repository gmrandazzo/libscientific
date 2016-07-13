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

typedef struct{
  dvector *x;
  double fval;
  double (*func)();
} simplex_arg;

/*
 *  NelderMeadSimplex:
 *  1) Perform k+1 steps with k the number of variable to optimise
 *  2) Rank the response from the best to the worst which is the last one x_wost
 *  3) Establish the new condition according to the function:
 *      x_new =  c - c - x_wost
 *      where c is the centroid calculate with all the best less the worst one.
 */

 void do_simplex(void *arg_)
 {
   simplex_arg *arg;
   arg = (simplex_arg*) arg_;
   //run the function and collect the results to argh->fval!
   arg->fval = arg->func(arg->x);
}

/*for maximization */
int cmpmax(const void *a, const void *b)
{
  const double *a_ = *(const double **)a;
  const double *b_ = *(const double **)b;
  int c = (int)(sizeof(a_)/sizeof(a_[0]))+1;
  //printf("cmpmax: %d\n", c);
  return b_[c] - a_[c];
}

/* for minimization */
int cmpmin(const void *a, const void *b)
{
  const double *a_ = *(const double **)a;
  const double *b_ = *(const double **)b;
  int c = (int)(sizeof(a_)/sizeof(a_[0]))+1;
  //printf("cmpmin: %d\n", c);
  return a_[c] - b_[c];
}

void gen_centroids(matrix *x, dvector *c)
{
  /* x is already ordered */
  size_t i, j;
  for(j = 0; j < x->col-1; j++){
    c->data[j] = 0.f;
    for(i = 0; i < x->row-1; i++){
      c->data[j] += x->data[i][j];
    }
    c->data[j] /= (double)(x->row-1.0);
  }
}

double NelderMeadSimplex(double (*func)(), dvector *x0, dvector *step, double xtol, size_t iter, dvector **best, enum OPT_TYPE otype)
{
  size_t i, j, k, iter_;
  matrix *x;
  dvector *c;
  double worst[2] = { 0.f, 0.f};
  double res = 0.f;
  simplex_arg arg;

  k = x0->size+1;
  if(otype == maximization)
    arg.fval = 0.f;
  else
    arg.fval = 9999.f;

  arg.func = func;
  NewDVector(&arg.x, x0->size);
  NewDVector(&c, x0->size);

  /*1) create k+1 steps */
  NewMatrix(&x, x0->size+1, x0->size+1);
  if(step != NULL){
    for(int i = 0; i < x0->size+1; i++){
      for(int j = 0; j < x0->size; j++){
        if(i-1 == j){
          x->data[i][j] = x0->data[j] + step->data[j];
        }
        else{
          x->data[i][j] = x0->data[j];
        }
      }
      if(otype == maximization)
        x->data[i][x0->size] = 0.f;
      else
        x->data[i][x0->size] = 9999.f; /* value of function */
    }
  }
  else{
    for(int i = 0; i < x0->size+1; i++){
      for(int j = 0; j < x0->size; j++){
        if(i-1 == j){
          x->data[i][j] = x0->data[j] + 0.5;
        }
        else{
          x->data[i][j] = x0->data[j];
        }
      }
      if(otype == maximization)
        x->data[i][x0->size] = 0.f;
      else
        x->data[i][x0->size] = 9999.f; /* value of function */
    }
  }

  /* initialization */
  for(i = 0; i < k; i++){
    if(otype == maximization)
      arg.fval = 0.f;
    else
      arg.fval = 9999.f;

    for(j = 0; j < x0->size; j++)
      arg.x->data[j] = x->data[i][j];
    arg.func = func;
    do_simplex(&arg);
    x->data[i][x0->size] = arg.fval;
  }
  if(otype == maximization)
    qsort(x->data, x->row, sizeof(double*), cmpmax);
  else
    qsort(x->data, x->row, sizeof(double*), cmpmin);

  /* the last is the worst one!! which is replaced with the new conditions */
  gen_centroids(x, c);
  worst[0] = worst[1] = x->data[k-1][k-1];

  for(j = 0; j < x0->size; j++){
    arg.x->data[j] = c->data[j] + 2.f*(c->data[j]-x->data[k-1][j]);
    x->data[k-1][j] = arg.x->data[j];
  }

  iter_ = 0;
  while(iter_ < iter)
  {
    printf("Worst: %.3f %.3f\n", worst[0], worst[1]);
    if(otype == maximization)
      arg.fval = 0.f;
    else
      arg.fval = 9999.f;
    for(j = 0; j < x0->size; j++)
      arg.x->data[j] = x->data[k-1][j];
    arg.func = func;
    do_simplex(&arg);
    x->data[k-1][k-1] = arg.fval;

    if(otype == maximization)
      qsort(x->data, x->row, sizeof(double*), cmpmax);
    else
      qsort(x->data, x->row, sizeof(double*), cmpmin);

    PrintMatrix(x);
    sleep(2);
    /* the last is the worst one!! which is replaced with the new conditions */
    gen_centroids(x, c);
    PrintDVector(c);

    /*check if the actual answer is better than all the other's */
    for(i = 0; i < x->row; i++){
      if(arg.fval > x->data[i][x->col-1]){
        continue;
      }
      else{
        break;
      }
    }

    /*if yes */
    if(i == x->row-1){
      x->data[x0->size][x0->size] = arg.fval;
      /* x = c + 2*(c-x) */
      puts("x = c + 2*(c-x)");
      for(j = 0; j < x0->size; j++){
        arg.x->data[j] = c->data[j] + 2.f*(c->data[j]-x->data[k-1][j]);
        x->data[k-1][j] = arg.x->data[j];
      }
    }
    else{
      if(arg.fval > worst[0] && arg.fval < worst[1]){
        puts("x = c + 0.5*(c-x)");
        /*x = c + 0.5*(c-x)*/
        for(j = 0; j < x0->size; j++){
          /* set the new condition and add to the last worst */
          arg.x->data[j] = c->data[j] + 0.5*(c->data[j] - x->data[k-1][j]);
          x->data[k-1][j] = arg.x->data[j];
        }
      }
      else if(arg.fval < worst[0] && arg.fval < worst[1]){
        puts("x = c - 0.5*(c-x)");
        /*x = c - 0.5*(c-x)*/
        for(j = 0; j < x0->size; j++){
          /* set the new condition and add to the last worst */
          arg.x->data[j] = c->data[j] - 0.5*(c->data[j] - x->data[k-1][j]);
          x->data[k-1][j] = arg.x->data[j];
        }
      }
      else{
        puts("x = c + c - x");
        /*x = c + c - x */
        for(j = 0; j < x0->size; j++){
          /* set the new condition and add to the last worst */
          arg.x->data[j] = c->data[j] + c->data[j] - x->data[k-1][j];
          x->data[k-1][j] = arg.x->data[j];
        }
      }
    }
    puts("suco");
    PrintMatrix(x);
    for(j = 0; j < x0->size; j++){
      /* set the new condition and add to the last worst */
      printf("%.3f - %.3f\n", c->data[j], x->data[k-1][j]);
      arg.x->data[j] = (c->data[j] + c->data[j]) - x->data[k-1][j];
      x->data[k-1][j] = arg.x->data[j];
    }

    double tmp = worst[0];
    worst[0] = x->data[k-1][k-1];
    worst[1] = tmp;

    if(fabs(x->data[0][k-1]-x->data[1][k-1]) > xtol){
      iter_ +=1;
    }
    else{
      break;
    }
  }

  if(otype == maximization)
    qsort(x->data, x->row, sizeof(double*), cmpmax);
  else
    qsort(x->data, x->row, sizeof(double*), cmpmin);


  DVectorResize(best, x0->size);
  for(i = 0; i < x0->size; i++)
    (*best)->data[i] = x->data[0][i];
  res = x->data[0][x0->size];

  DelDVector(&c);
  DelDVector(&arg.x);
  DelMatrix(&x);
  return res;
}

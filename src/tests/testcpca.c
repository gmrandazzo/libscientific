/* Unit tests for the cpca module.
 * Copyright (C) 2022-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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
#include <math.h>

#include "cpca.h"
#include "pca.h"
#include "tensor.h"
#include "numeric.h"

void PrepareMatrix4PCA(tensor *t, CPCAMODEL *model, matrix *X)
{
  size_t i, j, k, c, row, col;
  col = 0;
  row = 0;
  for(k = 0; k < t->order; k++){
    row = t->m[k]->row;
    col += t->m[k]->col;
  }

  ResizeMatrix(X, row, col);

  c = 0;
  for(k = 0; k < t->order; k++){
    double m = model->scaling_factor->data[k];
    for(j = 0; j < t->m[k]->col; j++){
      for(i = 0; i < t->m[k]->row; i++){
        X->data[i][c] = (t->m[k]->data[i][j] - model->colaverage->d[k]->data[j]);
        X->data[i][c] /= model->colscaling->d[k]->data[j];
        X->data[i][c] /= m; /* block scaling N.B.: this is applied to reproduce the pca score */
      }
      c += 1;
    }
  }
}

void test4()
{
  puts("Test 4: CPCA on 100x500x3 tensor ");
  size_t i, j, k;
  size_t nobj = 100;
  size_t nvars = 500;
  size_t norder = 3;
  tensor *t; // Definition of the input matrix
  CPCAMODEL *model; // Definition of the PCA model
  NewTensor(&t, norder);
  for(k = 0; k < norder; k++){
      NewTensorMatrix(t, k, nobj, nvars);
      srand_((uint32_t)time(0));
      for(i = 0; i < nobj; i++){
          for(j = 0; j < nvars; j++){
              t->m[k]->data[i][j] = randDouble(0,100);
          }
      }
  }
  NewCPCAModel(&model);
  CPCA(t, 1, 5, model);
  puts("CPCA on 100x500x3 tensor: OK");
  DelCPCAModel(&model);
  DelTensor(&t);
}

void test3()
{
  puts("Test 3: Projection test");
  size_t i, j, k;
  size_t rowsize=20;
  uivector *colsizes;
  tensor *t;
  matrix *p_super_scores;
  tensor *p_block_scores;

  initUIVector(&colsizes);
  UIVectorAppend(colsizes, 10);
  UIVectorAppend(colsizes, 7);
  UIVectorAppend(colsizes, 13);

  NewTensor(&t, colsizes->size);
  for(k = 0; k < colsizes->size; k++){
    NewTensorMatrix(t, k, rowsize, colsizes->data[k]);
    srand((unsigned int)(rowsize+colsizes->data[k]));
    for(i = 0; i < rowsize; i++){
      for(j = 0; j < colsizes->data[k]; j++){
        t->m[k]->data[i][j] = rand() % 100;
      }
    }
  }

  CPCAMODEL *model;
  NewCPCAModel(&model);

  CPCA(t, 1, 5, model);
  //PrintCPCA(model);

  initMatrix(&p_super_scores);
  initTensor(&p_block_scores);
  CPCAScorePredictor(t,
                     model,
                     5,
                     p_super_scores,
                     p_block_scores);

  for(i = 0; i < model->super_scores->row; i++){
    for(j = 0; j < model->super_scores->col; j++){
      if(FLOAT_EQ(model->super_scores->data[i][j], p_super_scores->data[i][j], 1E-4)){
        continue;
      }
      else{
        abort();
      }
    }
  }

  for(k = 0; k < model->block_scores->order; k++){
    for(i = 0; i < model->block_scores->m[k]->row; i++){
      for(j = 0; j < model->block_scores->m[k]->col; j++){
        if(FLOAT_EQ(model->block_scores->m[k]->data[i][j], p_block_scores->m[k]->data[i][j], 1E-4)){
          continue;
        }
        else{
          abort();
        }
      }
    }
  }

  puts("Projection test: OK");
  DelMatrix(&p_super_scores);
  DelTensor(&p_block_scores);
  DelCPCAModel(&model);
  DelTensor(&t);
  DelUIVector(&colsizes);
}

void test2()
{
  puts("Test 2: CPCA Variance test!");
 /*
  * we consider four blocks of five variables where one block
  * contains a strong direction that is not available in the other blocks.
  * X1 = [d1 d1 d1 d1 d1 ]
  * X2 = [d2 randn(4)]
  * X3 = [d2 randn(4)],
  * X4 = [d2 randn(4)]
  * All blocks have 50 observations
  * randn = random data
  *
  */
  size_t i, j, k;
  size_t rowsize=10;
  size_t colsize=5;
  size_t nblocks=4;
  tensor *t;

  /* d1 and d2 contains the cosine direction of the principal component */
  double d1 = -0.44;
  double d2 = +0.44;

  srand_((uint32_t)(rowsize+colsize+nblocks));
  NewTensor(&t, nblocks);
  for(k = 0; k < nblocks; k++){
    NewTensorMatrix(t, k, rowsize, colsize);
    for(i = 0; i < rowsize; i++){
      for(j = 0; j < colsize; j++){
        t->m[k]->data[i][j] = randDouble(-1, 1);
      }
    }
  }

  for(i = 0; i < rowsize; i++){
    for(j = 0; j < colsize; j++){
      t->m[0]->data[i][j] = (t->m[0]->data[i][j]*d1*d1);
    }

    t->m[1]->data[i][0] = (t->m[1]->data[i][0]*d2*d2);
    t->m[2]->data[i][0] = (t->m[2]->data[i][0]*d2*d2);
    t->m[3]->data[i][0] = (t->m[3]->data[i][0]*d2*d2);
  }


  CPCAMODEL *model;
  NewCPCAModel(&model);

  CPCA(t, 1, 5, model);
  //PrintCPCA(model);

  DelCPCAModel(&model);
  DelTensor(&t);
}

void test1(){
  puts("Test 1: Algorithm Check. CPCA super scores = PCA scores");
  size_t i, j, k;
  size_t rowsize=20;
  uivector *colsizes;
  tensor *t;
  size_t npc = 5;
  initUIVector(&colsizes);
  UIVectorAppend(colsizes, 10);
  UIVectorAppend(colsizes, 7);
  UIVectorAppend(colsizes, 13);

  NewTensor(&t, colsizes->size);
  for(k = 0; k < colsizes->size; k++){
    NewTensorMatrix(t, k, rowsize, colsizes->data[k]);
    srand((unsigned int)(rowsize+colsizes->data[k]));
    for(i = 0; i < rowsize; i++){
      for(j = 0; j < colsizes->data[k]; j++){
        t->m[k]->data[i][j] = rand() % 100;
      }
    }
  }

  CPCAMODEL *model;
  NewCPCAModel(&model);

  CPCA(t, 1, npc, model);
  //PrintCPCA(model);

  PCAMODEL *pca;
  matrix *X;
  initMatrix(&X);
  PrepareMatrix4PCA(t, model, X);

  NewPCAModel(&pca);
  PCA(X, -1, npc, pca, NULL);

  /* COMPARE THE SUPERSCORE WITH THE PCA RESULTS  */
  for(j = 0; j < model->super_scores->col; j++){
    for(i = 0; i < model->super_scores->row; i++){
      if(FLOAT_EQ(fabs(model->super_scores->data[i][j]), fabs(pca->scores->data[i][j]), 1E-1)){
        continue;
      }
      else{
        printf("%f =! %f\n", model->super_scores->data[i][j], pca->scores->data[i][j]);
        abort();
      }
    }
  }

  DelMatrix(&X);
  DelPCAModel(&pca);
  DelCPCAModel(&model);
  DelTensor(&t);
  DelUIVector(&colsizes);
}

/*
 * Sanity check: decomposition and recomposition must give back the data.
 *
 * Two blocks (3 and 4 columns) are generated from the same 3 latent
 * factors, so the concatenated data have rank 3. With npc = 3 (the maximum
 * allowed by the smallest block) the deflation removes everything, and each
 * block must be recomposed exactly from the super scores and its block
 * loadings:
 *
 *   X_b = (sum_k t_T[:,k] * p_b[:,k]') * colscaling_b + colaverage_b
 *
 * Every component is also checked to explain some variance and the total
 * explained variance must be 100%.
 */
void test5()
{
  printf("Test 5: CPCA decomposition and recomposition of every block = ");
  const size_t nobj = 30, nlat = 3;
  const size_t ncols[2] = {3, 4};
  size_t i, j, k, b;
  double T[30][3];
  tensor *t;
  CPCAMODEL *model;

  for(i = 0; i < nobj; i++){
    for(k = 0; k < nlat; k++)
      T[i][k] = sin(0.37*(double)((i+1)*(k+2))) * (double)(nlat - k);
  }

  NewTensor(&t, 2);
  for(b = 0; b < 2; b++){
    NewTensorMatrix(t, b, nobj, ncols[b]);
    for(i = 0; i < nobj; i++){
      for(j = 0; j < ncols[b]; j++){
        double x = 10.0*(double)(b+1); /* offset: checks the centering */
        for(k = 0; k < nlat; k++)
          x += T[i][k] * cos(1.3*(double)((j+1)*(k+1)) + (double)b);
        t->m[b]->data[i][j] = x;
      }
    }
  }

  NewCPCAModel(&model);
  CPCA(t, 1, nlat, model);

  double totvar = 0.0;
  for(k = 0; k < nlat; k++){
    if(model->total_expvar->data[k] <= 1e-6){
      printf("FAIL: PC%zu explains %g%% of variance\n", k+1, model->total_expvar->data[k]);
      abort();
    }
    totvar += model->total_expvar->data[k];
  }
  if(!FLOAT_EQ(totvar, 100.0, 1e-4)){
    printf("FAIL: total explained variance %g%%\n", totvar);
    abort();
  }

  for(b = 0; b < 2; b++){
    for(i = 0; i < nobj; i++){
      for(j = 0; j < ncols[b]; j++){
        double x = 0.0;
        for(k = 0; k < nlat; k++)
          x += model->super_scores->data[i][k] * model->block_loadings->m[b]->data[j][k];
        x = x * model->colscaling->d[b]->data[j] + model->colaverage->d[b]->data[j];
        if(!FLOAT_EQ(x, t->m[b]->data[i][j], 1e-6)){
          printf("FAIL: block %zu (%zu,%zu) recomposed %g != %g\n", b, i, j, x, t->m[b]->data[i][j]);
          abort();
        }
      }
    }
  }

  printf("OK.\n");
  DelCPCAModel(&model);
  DelTensor(&t);
}


static void fill_refit_data(matrix *m, double offset, size_t seed)
{
  for(size_t i = 0; i < m->row; i++)
    for(size_t j = 0; j < m->col; j++)
      m->data[i][j] = offset + sin(0.7*(double)((i+1)*(j+seed)) + (double)seed);
}

/*
 * Regression test: fitting the same model twice must give the same result as
 * a new model (block results used to be appended and the column averages of
 * the first fit reused).
 */
void test6()
{
  printf("Test 6: CPCA refit a model on new data = ");
  tensor *a, *b;
  CPCAMODEL *refit, *fresh;
  NewTensor(&a, 2); NewTensor(&b, 2);
  for(size_t k = 0; k < 2; k++){
    NewTensorMatrix(a, k, 12, 3 + k); fill_refit_data(a->m[k], 0.0, 1 + k);
    NewTensorMatrix(b, k, 12, 3 + k); fill_refit_data(b->m[k], 5.0, 3 + k);
  }
  NewCPCAModel(&refit); NewCPCAModel(&fresh);
  CPCA(a, 1, 2, refit);
  CPCA(b, 1, 2, refit);
  CPCA(b, 1, 2, fresh);
  if(refit->total_expvar->size != fresh->total_expvar->size ||
     refit->block_scores->order != fresh->block_scores->order ||
     refit->colaverage->size != fresh->colaverage->size){
    printf("FAIL: results were appended instead of replaced\n");
    abort();
  }
  for(size_t k = 0; k < 2; k++)
    for(size_t j = 0; j < b->m[k]->col; j++)
      if(!FLOAT_EQ(refit->colaverage->d[k]->data[j], fresh->colaverage->d[k]->data[j], 1e-12)){
        printf("FAIL: block %zu column %zu average differs after refit\n", k, j);
        abort();
      }
  for(size_t i = 0; i < 12; i++)
    for(size_t c = 0; c < 2; c++)
      if(!FLOAT_EQ(refit->super_scores->data[i][c], fresh->super_scores->data[i][c], 1e-9)){
        printf("FAIL: super score (%zu,%zu) differs after refit\n", i, c);
        abort();
      }
  printf("OK.\n");
  DelCPCAModel(&refit); DelCPCAModel(&fresh);
  DelTensor(&a); DelTensor(&b);
}

int main(void){
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  return 0;
}

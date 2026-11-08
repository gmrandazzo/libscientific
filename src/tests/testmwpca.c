/* Unit tests for the mwpca module.
 * Copyright (C) 2016-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "tensor.h"
#include "upca.h"
#include "numeric.h"

void test5()
{
  size_t i, j;
  tensor *t;
  matrix *ps;
  UPCAMODEL *m;

  puts(">>>>>>> Test 5: Predicted scores of the training tensor = model scores");

  NewTensor(&t, 2);

  NewTensorMatrix(t, 0, 3, 2);
  NewTensorMatrix(t, 1, 3, 2);

  setTensorValue(t, 0, 0, 0, 0.424264);  setTensorValue(t, 0, 0, 1, 0.565685);
  setTensorValue(t, 0, 1, 0, 0.565685);  setTensorValue(t, 0, 1, 1, 0.424264);
  setTensorValue(t, 0, 2, 0, 0.707101);  setTensorValue(t, 0, 2, 1, 0.707101);

  setTensorValue(t, 1, 0, 0, 0.565685);  setTensorValue(t, 1, 0, 1, 0.424264);
  setTensorValue(t, 1, 1, 0, 0.424264);  setTensorValue(t, 1, 1, 1, 0.565685);
  setTensorValue(t, 1, 2, 0, 0.707101);  setTensorValue(t, 1, 2, 1, 0.707101);

  puts("Tensor");
  PrintTensor(t);

  NewUPCAModel(&m);

  UPCA(t, 2, 1, m, NULL);
  initMatrix(&ps);
  UPCAScorePredictor(t, m,2, ps);

  for(i = 0; i < m->scores->row; i++){
    for(j = 0; j < m->scores->col; j++){
      if(FLOAT_EQ(m->scores->data[i][j], ps->data[i][j], 1E-6)){
        continue;
      }
      else{
        abort();
      }
    }
  }

  DelUPCAModel(&m);
  DelMatrix(&ps);
  DelTensor(&t);
}

void test4()
{
  tensor *t;
  UPCAMODEL *m;

  puts(">>>>>>> Test 4: Compute Multi Way PCA");

  NewTensor(&t, 2);

  NewTensorMatrix(t, 0, 4, 2);
  NewTensorMatrix(t, 1, 4, 2);

  setTensorValue(t, 0, 0, 0, 7);  setTensorValue(t, 0, 0, 1, 9);
  setTensorValue(t, 0, 1, 0, 5);  setTensorValue(t, 0, 1, 1, 4);
  setTensorValue(t, 0, 2, 0, 2);  setTensorValue(t, 0, 2, 1, 10);
  setTensorValue(t, 0, 3, 0, 4);  setTensorValue(t, 0, 3, 1, 16);

  setTensorValue(t, 1, 0, 0, 7);  setTensorValue(t, 1, 0, 1, 9);
  setTensorValue(t, 1, 1, 0, 5);  setTensorValue(t, 1, 1, 1, 4);
  setTensorValue(t, 1, 2, 0, 4);  setTensorValue(t, 1, 2, 1, 10);
  setTensorValue(t, 1, 3, 0, 6);  setTensorValue(t, 1, 3, 1, 16);

  puts("Tensor");
  PrintTensor(t);

  NewUPCAModel(&m);

  UPCA(t, 2, 1, m, NULL);

  DelUPCAModel(&m);
  DelTensor(&t);
}

/*Test from:
 *
 * PCA Nipals algorithm for Multi-Way Principal Component Analysis
 * Journal of Chemometrics, Vol 1, 41-56 (1987)
 * Svante Woold, Paul Geladi and Kim Esbensen
 */
void test3()
{
  tensor *t;

  UPCAMODEL *m;

  tensor *tpred;
  matrix *predscore;

  puts(">>>>>>> Test 3: Compute Multi Way PCA and Prediction");

  NewTensor(&t, 2);

  NewTensorMatrix(t, 0, 3, 2);
  NewTensorMatrix(t, 1, 3, 2);

  setTensorValue(t, 0, 0, 0, 0.424264);  setTensorValue(t, 0, 0, 1, 0.565685);
  setTensorValue(t, 0, 1, 0, 0.565685);  setTensorValue(t, 0, 1, 1, 0.424264);
  setTensorValue(t, 0, 2, 0, 0.707101);  setTensorValue(t, 0, 2, 1, 0.707101);

  setTensorValue(t, 1, 0, 0, 0.565685);  setTensorValue(t, 1, 0, 1, 0.424264);
  setTensorValue(t, 1, 1, 0, 0.424264);  setTensorValue(t, 1, 1, 1, 0.565685);
  setTensorValue(t, 1, 2, 0, 0.707101);  setTensorValue(t, 1, 2, 1, 0.707101);


  NewTensor(&tpred, 2);
  NewTensorMatrix(tpred, 0, 1, 2);
  NewTensorMatrix(tpred, 1, 1, 2);

  setTensorValue(tpred, 0, 0, 0, 0.5);  setTensorValue(tpred, 0, 0, 1, 0.6);

  setTensorValue(tpred, 1, 0, 0, 0.6);  setTensorValue(tpred, 1, 0, 1, 0.4);


  puts("Tensor");
  PrintTensor(t);


  NewUPCAModel(&m);

  UPCA(t, 2, 1, m, NULL);

  PrintUPCAModel(m);


  puts("Tensor to pred");
  PrintTensor(tpred);

  initMatrix(&predscore);
  UPCAScorePredictor(tpred, m, 12345, predscore);

  puts("Predicted Scores");
  PrintMatrix(predscore);

  DelMatrix(&predscore);
  DelTensor(&tpred);

  DelUPCAModel(&m);
  DelTensor(&t);
}


/*
 * Test with Random numbers
 */
void test2()
{
  size_t i, j, k;
  tensor *t;

  UPCAMODEL *m;

  puts(">>>>>>> Test 2: Compute Multi Way PCA with random variables");

  NewTensor(&t, 5);

  for(i = 0; i < t->order; i++){
    NewTensorMatrix(t, i, 30, 10);
  }

  srand((unsigned int)time(0));
  for(i = 0; i < t->order; i++){
    for(j = 0; j < t->m[i]->row; j++){
      for(k = 0; k < t->m[i]->col; k++)
          setTensorValue(t, i, j, k, rand()/((double)(RAND_MAX)+1));
    }
  }

  /*
  puts("Tensor");
  PrintTensor(t);
  */

  NewUPCAModel(&m);
  scisignal run = SIGSCIENTIFICRUN;
  UPCA(t, 3, 1, m, &run);

  PrintUPCAModel(m);

  DelUPCAModel(&m);
  DelTensor(&t);
}

/*Test from:
 *
 * PCA Nipals algorithm for Multi-Way Principal Component Analysis
 * Journal of Chemometrics, Vol 1, 41-56 (1987)
 * Svante Woold, Paul Geladi and Kim Esbensen
 */
void test1()
{
  tensor *t;
  UPCAMODEL *m;

  puts(">>>>>>> Test 1: Compute Multi Way PCA");

  NewTensor(&t, 2);

  NewTensorMatrix(t, 0, 3, 2);
  NewTensorMatrix(t, 1, 3, 2);

  setTensorValue(t, 0, 0, 0, 0.424264);  setTensorValue(t, 0, 0, 1, 0.565685);
  setTensorValue(t, 0, 1, 0, 0.565685);  setTensorValue(t, 0, 1, 1, 0.424264);
  setTensorValue(t, 0, 2, 0, 0.707101);  setTensorValue(t, 0, 2, 1, 0.707101);

  setTensorValue(t, 1, 0, 0, 0.565685);  setTensorValue(t, 1, 0, 1, 0.424264);
  setTensorValue(t, 1, 1, 0, 0.424264);  setTensorValue(t, 1, 1, 1, 0.565685);
  setTensorValue(t, 1, 2, 0, 0.707101);  setTensorValue(t, 1, 2, 1, 0.707101);

  puts("Tensor");
  PrintTensor(t);

  NewUPCAModel(&m);

  UPCA(t, 2, 1, m, NULL);

  PrintUPCAModel(m);

  DelUPCAModel(&m);
  DelTensor(&t);
}

/*
 * Sanity checks:
 * - a rank-1 trilinear tensor (X[k][i][j] = a_i * b_j * c_k plus column
 *   offsets) is described exactly by one UPCA component, so decomposing
 *   and recomposing it must give back the input, with and without scaling
 *   (the scaling used to be reverted as x + scaling + average);
 * - refitting a model must replace it (scores and loadings used to be
 *   appended and the old column averages reused).
 */
void test6()
{
  puts(">>>>>>> Test 6: UPCA decomposition/recomposition and refit");
  tensor *X, *Xr;
  UPCAMODEL *m;
  NewTensor(&X, 3);
  for(size_t k = 0; k < 3; k++){
    NewTensorMatrix(X, k, 15, 4);
    for(size_t i = 0; i < 15; i++)
      for(size_t j = 0; j < 4; j++)
        X->m[k]->data[i][j] = 3.0*(double)(j+1) + sin(0.9*(double)(i+1)) * (1.0+0.5*(double)j) * (1.0+0.7*(double)k);
  }
  for(size_t scaling = 0; scaling <= 1; scaling++){
    NewUPCAModel(&m);
    UPCA(X, 1, scaling, m, NULL);
    UPCA(X, 1, scaling, m, NULL); /* refit */
    if(m->scores->col != 1 || m->loadings->order != 1 || m->colaverage->size != 3){
      printf("FAIL: after refit scores %zux%zu, %zu loadings, %zu averages\n",
             m->scores->row, m->scores->col, m->loadings->order, m->colaverage->size);
      abort();
    }
    if(!FLOAT_EQ(m->varexp->data[0], 100.0, 1e-6)){
      printf("FAIL: one component explains %g%%, expected 100%%\n", m->varexp->data[0]);
      abort();
    }
    initTensor(&Xr);
    UPCAIndVarPredictor(m->scores, m->loadings, m->colaverage, m->colscaling, 1, Xr);
    for(size_t k = 0; k < 3; k++)
      for(size_t i = 0; i < 15; i++)
        for(size_t j = 0; j < 4; j++)
          if(!FLOAT_EQ(Xr->m[k]->data[i][j], X->m[k]->data[i][j], 1e-8)){
            printf("FAIL: scaling %zu: X[%zu][%zu][%zu] recomposed %g != %g\n",
                   scaling, k, i, j, Xr->m[k]->data[i][j], X->m[k]->data[i][j]);
            abort();
          }
    DelTensor(&Xr);
    DelUPCAModel(&m);
  }
  DelTensor(&X);
  puts("OK.");
}

int main(void)
{
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  return 0;
}

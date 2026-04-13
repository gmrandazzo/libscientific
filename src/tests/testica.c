/* Unit tests for the ica module.
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "numeric.h"
#include "ica.h"
#include "datasets.h"
#include "scientificinfo.h"
#include "io.h"

void test3()
{
  puts("Test ICA 3: ICA Prediction");
  size_t i, j;
  matrix *m, *_;
  matrix *psignals;
  ICAMODEL *model;

  initMatrix(&m);
  initMatrix(&_);
  iris(m, _);

  NewICAModel(&model);

  ICA(m, 1, 2, model);
  initMatrix(&psignals);
  ICASignalPredictor(m, model, psignals);
  for(i = 0; i < model->S->row; i++){
    for(j = 0; j < model->S->col; j++){
      if(FLOAT_EQ(model->S->data[i][j], psignals->data[i][j], 1e-12)){
        continue;
      }
      else{
        abort();
      }
    }
  }
  puts("OK.");
  DelICAModel(&model);
  DelMatrix(&m);
  DelMatrix(&psignals);
  DelMatrix(&_);
}


void test2()
{
  puts("Test ICA 2: ICA on iris dataset");
  matrix *m, *_;
  ICAMODEL *model;

  initMatrix(&m);
  initMatrix(&_);
  iris(m, _);

  NewICAModel(&model);

  ICA(m, 1, 2, model);
  PrintICA(model);
  puts("OK.");
  DelICAModel(&model);
  DelMatrix(&m);
  DelMatrix(&_);
}


void test1()
{
  puts("Test ICA 1: ICA on a random 10x3 dataset");
  matrix *m;
  ICAMODEL *model;

  NewMatrix(&m, 10, 3);
  m->data[0][0] = -1.; m->data[0][1] = -1.; m->data[0][2] = -2.;
  m->data[1][0] = 2.75643348; m->data[1][1] = 3.26710563; m->data[1][2] = 4.02353911;
  m->data[2][0] = -0.84668509; m->data[2][1] = -1.64556477; m->data[2][2] = -0.49224986;
  m->data[3][0] = 0.52000394; m->data[3][1] = 1.92666864; m->data[3][2] = 0.44667258;
  m->data[4][0] = -0.15235893; m->data[4][1] = -1.52062391;  m->data[4][2] = 0.32701716;
  m->data[5][0] = 1.39949457; m->data[5][1] = 2.14419173; m->data[5][2] = 1.5436863;
  m->data[6][0] = -2.27972909; m->data[6][1] = -2.80653121; m->data[6][2] = -3.0862603;
  m->data[7][0] = -1.67717986; m->data[7][1] = -2.61636771; m->data[7][2] = -2.29354756;
  m->data[8][0] = 1.2186072;  m->data[8][1] = 1.72041471; m->data[8][2] = 0.93902191;
  m->data[9][0] = -2.28790332; m->data[9][1] = -3.14395166; m->data[9][2] = -3.43185497;

  NewICAModel(&model);
  PrintMatrix(m);
  ICA(m, 0, 2, model);
  PrintICA(model);
  puts("OK.");
  DelICAModel(&model);
  DelMatrix(&m);
}

void test4()
{
  puts("Test ICA 4: Signal separation and CSV output");
  size_t n_samples = 200;
  size_t n_sources = 2;
  matrix *S_orig; /* Original sources */
  matrix *X;      /* Mixed signals */
  matrix *A;      /* Mixing matrix */
  ICAMODEL *model;

  NewMatrix(&S_orig, n_samples, n_sources);
  for(size_t i = 0; i < n_samples; i++){
    double t = (double)i / 20.0;
    S_orig->data[i][0] = sin(t); /* Sine wave */
    S_orig->data[i][1] = (sin(2.0*t) > 0) ? 1.0 : -1.0; /* Square wave */
  }

  /* Mixing matrix */
  NewMatrix(&A, n_sources, n_sources);
  A->data[0][0] = 0.5; A->data[0][1] = 0.5;
  A->data[1][0] = 0.2; A->data[1][1] = 0.8;

  /* Mix signals: X = S_orig * A^T (since S_orig is n_samples x n_sources) */
  NewMatrix(&X, n_samples, n_sources);
  for(size_t i = 0; i < n_samples; i++){
    for(size_t j = 0; j < n_sources; j++){
      for(size_t k = 0; k < n_sources; k++){
        X->data[i][j] += S_orig->data[i][k] * A->data[j][k];
      }
    }
  }

  NewICAModel(&model);
  ICA_ext(X, 1, n_sources, 1.0, 1e-8, 5000, model);

  /* Output original, mixed and separated signals to CSV */
  matrix *out;
  NewMatrix(&out, n_samples, n_sources * 3);
  for(size_t i = 0; i < n_samples; i++){
    out->data[i][0] = S_orig->data[i][0];
    out->data[i][1] = S_orig->data[i][1];
    out->data[i][2] = X->data[i][0];
    out->data[i][3] = X->data[i][1];
    out->data[i][4] = model->S->data[i][0];
    out->data[i][5] = model->S->data[i][1];
  }

  WriteMatrixCSV("ica_signals.csv", out);
  puts("Results saved to ica_signals.csv. OK.");

  DelMatrix(&out);
  DelICAModel(&model);
  DelMatrix(&X);
  DelMatrix(&A);
  DelMatrix(&S_orig);
}


int main(void)
{
  test1();
  /*test2();
  test3();*/
  test4();
  return 0;
}
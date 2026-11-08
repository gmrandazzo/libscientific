/* Unit tests for the pca module.
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

#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include "numeric.h"
#include "pca.h"
#include "pca.h"
#include "datasets.h"
#include "scientificinfo.h"

void test11()
{
  printf("Test 11 - PCA T squared Contributions = ");
  matrix *m; /* Data matrix */
  matrix *contributions;
  dvector *spe;
  dvector *spe_ground_truth;
  PCAMODEL *model;
  int run = SIGSCIENTIFICRUN;

  NewMatrix(&m, 14, 7);

  m->data[0][0] = 4.0000;  m->data[0][1] = 4.0000;  m->data[0][2] = 1.0000;  m->data[0][3] = 84.1400;  m->data[0][4] = 1.0500;  m->data[0][5] = 235.1500;  m->data[0][6] = 357.1500;
  m->data[1][0] = 5.0000;  m->data[1][1] = 5.0000;  m->data[1][2] = 1.0000;  m->data[1][3] = 79.1000;  m->data[1][4] = 0.9780;  m->data[1][5] = 1.5090;  m->data[1][6] = 231.0000;
  m->data[2][0] = 4.0000;  m->data[2][1] = 5.0000;  m->data[2][2] = 1.0000;  m->data[2][3] = 67.0900;  m->data[2][4] = 0.9700;  m->data[2][5] = 249.0000;  m->data[2][6] = 403.0000;
  m->data[3][0] = 4.0000;  m->data[3][1] = 4.0000;  m->data[3][2] = 1.0000;  m->data[3][3] = 68.0700;  m->data[3][4] = 0.9360;  m->data[3][5] = 187.3500;  m->data[3][6] = 304.5500;
  m->data[4][0] = 3.0000;  m->data[4][1] = 4.0000;  m->data[4][2] = 2.0000;  m->data[4][3] = 68.0800;  m->data[4][4] = 1.0300;  m->data[4][5] = 363.0000;  m->data[4][6] = 529.0000;
  m->data[5][0] = 9.0000;  m->data[5][1] = 7.0000;  m->data[5][2] = 1.0000;  m->data[5][3] = 129.1600;  m->data[5][4] = 1.0900;  m->data[5][5] = 258.0000;  m->data[5][6] = 510.0000;
  m->data[6][0] = 10.0000;  m->data[6][1] = 8.0000;  m->data[6][2] = 0.0000;  m->data[6][3] = 128.1600;  m->data[6][4] = 1.1500;  m->data[6][5] = 352.0000;  m->data[6][6] = 491.0000;
  m->data[7][0] = 6.0000;  m->data[7][1] = 6.0000;  m->data[7][2] = 0.0000;  m->data[7][3] = 78.1118;  m->data[7][4] = 0.8765;  m->data[7][5] = 278.6400;  m->data[7][6] = 353.3000;
  m->data[8][0] = 16.0000;  m->data[8][1] = 10.0000;  m->data[8][2] = 0.0000;  m->data[8][3] = 202.2550;  m->data[8][4] = 1.2710;  m->data[8][5] = 429.1500;  m->data[8][6] = 666.6500;
  m->data[9][0] = 6.0000;  m->data[9][1] = 12.0000;  m->data[9][2] = 0.0000;  m->data[9][3] = 84.1600;  m->data[9][4] = 0.7800;  m->data[9][5] = 279.0000;  m->data[9][6] = 354.0000;
  m->data[10][0] = 4.0000;  m->data[10][1] = 8.0000;  m->data[10][2] = 1.0000;  m->data[10][3] = 72.1100;  m->data[10][4] = 0.8900;  m->data[10][5] = 164.5000;  m->data[10][6] = 339.0000;
  m->data[11][0] = 4.0000;  m->data[11][1] = 9.0000;  m->data[11][2] = 1.0000;  m->data[11][3] = 71.1100;  m->data[11][4] = 0.8660;  m->data[11][5] = 210.0000;  m->data[11][6] = 360.0000;
  m->data[12][0] = 5.0000;  m->data[12][1] = 11.0000;  m->data[12][2] = 1.0000;  m->data[12][3] = 85.1500;  m->data[12][4] = 0.8620;  m->data[12][5] = 266.0000;  m->data[12][6] = 379.0000;
  m->data[13][0] = 5.0000;  m->data[13][1] = 10.0000;  m->data[13][2] = 1.0000;  m->data[13][3] = 86.1300;  m->data[13][4] = 0.8800;  m->data[13][5] = 228.0000;  m->data[13][6] = 361.0000;

  NewDVector(&spe_ground_truth, 14);
  spe_ground_truth->data[0] = 0.0492;
  spe_ground_truth->data[1] = 0.0059;
  spe_ground_truth->data[2] = 0.0407;
  spe_ground_truth->data[3] = 0.0555;
  spe_ground_truth->data[4] = 0.0025;
  spe_ground_truth->data[5] = 0.0049;
  spe_ground_truth->data[6] = 0.0312;
  spe_ground_truth->data[7] = 0.0017;
  spe_ground_truth->data[8] = 0.0070;
  spe_ground_truth->data[9] = 0.0030;
  spe_ground_truth->data[10] = 0.0242;
  spe_ground_truth->data[11] = 0.0115;
  spe_ground_truth->data[12] = 0.0236;
  spe_ground_truth->data[13] = 0.0166;

  NewPCAModel(&model);

  PrintMatrix(m);
  PCA(m, 1, 5, model, &run);
  initDVector(&spe);
  initMatrix(&contributions);
  PCATsqContributions(m, model, 5, spe, contributions);
  PrintDVector(spe);
  PrintMatrix(contributions);
  for(size_t i = 0; i < spe->size; i++){
    if(FLOAT_EQ(spe->data[i], spe_ground_truth->data[i], 1e-4)){
      continue;
    }
    else{
      printf("%f %f\n", spe->data[i], spe_ground_truth->data[i]);
      abort();
    }
  }
  printf("OK.\n");

  DelDVector(&spe_ground_truth);
  DelDVector(&spe);
  DelMatrix(&contributions);
  DelPCAModel(&model);
  DelMatrix(&m);
}

void test10(){
  printf("Test PCA 10 - Check PCAIndVarPredictor  = ");
  size_t i, j;
  matrix *m, *_, *rm;
  PCAMODEL *model;

  initMatrix(&m);
  initMatrix(&_);
  iris(m, _);

  NewPCAModel(&model);
  PCA(m, 1, 4, model, NULL);

  initMatrix(&rm);
  PCAIndVarPredictor(model->scores, model->loadings, model->colaverage, model->colscaling, 4, rm);

  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(FLOAT_EQ(m->data[i][j], rm->data[i][j], 1E-6)){
        continue;
      }
      else{
        abort();
      }
    }
  }
  printf("OK.\n");
  DelPCAModel(&model);
  DelMatrix(&rm);
  DelMatrix(&m);
  DelMatrix(&_);
}

void test9()
{
  printf("Test PCA 9 - Check algorithm consistency = ");
  size_t i;
  size_t id;
  matrix *m1;
  matrix *m2;
  uivector *ids;
  PCAMODEL *mod1;
  PCAMODEL *mod2;
  NewMatrix(&m1, 20, 3);
  for(i = 0; i < 20; i++){
    m1->data[i][0] = (double)(i+1); m1->data[i][1] = (double)(i+1); m1->data[i][2] = (double)(i+1);
  }
  srand_((uint32_t)time(NULL));
  NewMatrix(&m2, m1->row, m1->col);
  initUIVector(&ids);
  i = 0;
  while( i < 20){
    id = (size_t)randInt(0, 20);
    if(UIVectorHasValue(ids, id) == 0){
      continue;
    }
    else{
      UIVectorAppend(ids, id);
      m2->data[i][0] = m1->data[id][0]; m2->data[i][1] = m1->data[id][1]; m2->data[i][2] = m1->data[id][2];
      i++;
    }
  }
  NewPCAModel(&mod1);
  NewPCAModel(&mod2);
  PCA(m1, 1, 2, mod1, NULL);
  PCA(m2, 1, 2, mod2, NULL);
  for(i = 0; i < mod1->scores->row; i++){
    id = ids->data[i];
    if(FLOAT_EQ(mod1->scores->data[id][0], mod2->scores->data[i][0], 1e-2) &&
       FLOAT_EQ(mod1->scores->data[id][1], mod2->scores->data[i][1], 1e-2)){
      continue;
    }
    else{
      printf("%f %f\n", mod1->scores->data[id][0], mod2->scores->data[i][0]);
      printf("%f %f\n", mod1->scores->data[id][1], mod2->scores->data[i][1]);
      abort();
    }
  }
  printf("OK.\n");
  DelPCAModel(&mod1);
  DelPCAModel(&mod2);
  DelUIVector(&ids);
  DelMatrix(&m1);
  DelMatrix(&m2);
}

void test8()
{
  printf("Test PCA 8 - Score prediction test = ");
  size_t i, j;
  matrix *m, *_, *p;
  PCAMODEL *model;

  initMatrix(&m);
  initMatrix(&_);
  iris(m, _);

  NewPCAModel(&model);

  PCA(m, 1, 2, model, NULL);
  initMatrix(&p);
  PCAScorePredictor(m, model, 2, p);


  for(i = 0; i < model->scores->row; i++){
    for(j = 0; j < model->scores->col; j++){
      if(FLOAT_EQ(model->scores->data[i][j], p->data[i][j], 1E-6)){
        continue;
      }
      else{
        abort();
      }
    }
  }
  printf("OK.\n");
  DelPCAModel(&model);
  DelMatrix(&p);
  DelMatrix(&m);
  DelMatrix(&_);
}

void test7()
{
  printf("Test PCA 7 - PCA on iris dataset = ");
  matrix *m, *_;
  PCAMODEL *model;

  initMatrix(&m);
  initMatrix(&_);
  iris(m, _);

  NewPCAModel(&model);

  PCA(m, 1, 2, model, NULL);
  PrintPCA(model);
  DelPCAModel(&model);
  DelMatrix(&m);
  DelMatrix(&_);
}

void test6()
{
  puts("Test PCA 6 = PCA on a simple dataset = ");
  matrix *m;
  PCAMODEL *model;

  NewMatrix(&m, 3,3);

  m->data[0][0] = 7; m->data[0][1] = 4; m->data[0][2] = 4;
  m->data[1][0] = 4; m->data[1][1] = 4; m->data[1][2] = 4;
  m->data[2][0] = 1; m->data[2][1] = 4; m->data[2][2] = 7;

  PrintMatrix(m);

  NewPCAModel(&model);

  PCA(m, 0, 5, model, NULL);
  printf("OK.\n");
  PrintPCA(model);

  DelPCAModel(&model);
  DelMatrix(&m);
}

void test5()
{
  printf("Test PCA 5 - PCA on a simple dataset = ");
  matrix *m; /* Data matrix */
  PCAMODEL *model;

  NewMatrix(&m, 3, 4);

  setMatrixValue(m, 0, 0, 0.424264);  setMatrixValue(m, 0, 1, 0.565685);  setMatrixValue(m, 0, 2, 0.565685);  setMatrixValue(m, 0, 3, 0.424264);
  setMatrixValue(m, 1, 0, 0.565685);  setMatrixValue(m, 1, 1, 0.424264);  setMatrixValue(m, 1, 2, 0.424264);  setMatrixValue(m, 1, 3, 0.565685);
  setMatrixValue(m, 2, 0, 0.707101);  setMatrixValue(m, 2, 1, 0.707101);  setMatrixValue(m, 2, 2, 0.707101);  setMatrixValue(m, 2, 3, 0.707101);


  NewPCAModel(&model);
  PCA(m, 1, 5, model, NULL);
  printf("OK.\n");
  PrintPCA(model);

  DelPCAModel(&model);
  DelMatrix(&m);
}

void test4()
{
  printf("Test PCA 4 - PCA on a simple dataset = ");
  matrix *m; /* Data matrix */
  PCAMODEL *model;
  int run = SIGSCIENTIFICRUN;
  NewMatrix(&m, 3, 2);

  setMatrixValue(m, 0, 0, 0.424264);  setMatrixValue(m, 0, 1, 0.565685);
  setMatrixValue(m, 1, 0, 0.565685);  setMatrixValue(m, 1, 1, 0.424264);
  setMatrixValue(m, 2, 0, 0.707101);  setMatrixValue(m, 2, 1, 0.707101);

  NewPCAModel(&model);

  PCA(m, 1, 5, model, &run);
  printf("OK.\n");
  PrintPCA(model);
  DelPCAModel(&model);
  DelMatrix(&m);
}

void test3()
{
  printf("Test PCA 3 - PCA on a simple dataset = ");
  matrix *m; /* Data matrix */
  PCAMODEL *model;
  int run = SIGSCIENTIFICRUN;
  NewMatrix(&m, 10, 2);
  m->data[0][0] = 2.2;  m->data[0][1] = 3.3;
  m->data[1][0] = 2.9;  m->data[1][1] = 4.5;
  m->data[2][0] = 4.0;  m->data[2][1] = 3.5;
  m->data[3][0] = 4.500;  m->data[3][1] = 4.7;
  m->data[4][0] = 5.1;  m->data[4][1] = 6.0;
  m->data[5][0] = 5.5;  m->data[5][1] = 4.5;
  m->data[6][0] = 6.5000;  m->data[6][1] = 6.4;
  m->data[7][0] = 7.2;  m->data[7][1] = 6.2;
  m->data[8][0] = 8.8;  m->data[8][1] = 3.9;
  m->data[9][0] = 10;  m->data[9][1] = 6.2;


  NewPCAModel(&model);
  PCA(m, 1, 5, model, &run);
  printf("OK.\n");
  PrintPCA(model);

  DelPCAModel(&model);
  DelMatrix(&m);
}


void test2()
{
  printf("Test PCA 2 - PCARankValidation test = ");
  size_t i;
  matrix *m; /* Data matrix */
  dvector *R;
  int run = SIGSCIENTIFICSTOP;
  NewMatrix(&m, 14, 7);

  m->data[0][0] = 4.0000;  m->data[0][1] = 4.0000;  m->data[0][2] = 1.0000;  m->data[0][3] = 84.1400;  m->data[0][4] = 1.0500;  m->data[0][5] = 235.1500;  m->data[0][6] = 357.1500;
  m->data[1][0] = 5.0000;  m->data[1][1] = 5.0000;  m->data[1][2] = 1.0000;  m->data[1][3] = 79.1000;  m->data[1][4] = 0.9780;  m->data[1][5] = 1.5090;  m->data[1][6] = 231.0000;
  m->data[2][0] = 4.0000;  m->data[2][1] = 5.0000;  m->data[2][2] = 1.0000;  m->data[2][3] = 67.0900;  m->data[2][4] = 0.9700;  m->data[2][5] = 249.0000;  m->data[2][6] = 403.0000;
  m->data[3][0] = 4.0000;  m->data[3][1] = 4.0000;  m->data[3][2] = 1.0000;  m->data[3][3] = 68.0700;  m->data[3][4] = 0.9360;  m->data[3][5] = 187.3500;  m->data[3][6] = 304.5500;
  m->data[4][0] = 3.0000;  m->data[4][1] = 4.0000;  m->data[4][2] = 2.0000;  m->data[4][3] = 68.0800;  m->data[4][4] = 1.0300;  m->data[4][5] = 363.0000;  m->data[4][6] = 529.0000;
  m->data[5][0] = 9.0000;  m->data[5][1] = 7.0000;  m->data[5][2] = 1.0000;  m->data[5][3] = 129.1600;  m->data[5][4] = 1.0900;  m->data[5][5] = 258.0000;  m->data[5][6] = 510.0000;
  m->data[6][0] = 10.0000;  m->data[6][1] = 8.0000;  m->data[6][2] = 0.0000;  m->data[6][3] = 128.1600;  m->data[6][4] = 1.1500;  m->data[6][5] = 352.0000;  m->data[6][6] = 491.0000;
  m->data[7][0] = 6.0000;  m->data[7][1] = 6.0000;  m->data[7][2] = 0.0000;  m->data[7][3] = 78.1118;  m->data[7][4] = 0.8765;  m->data[7][5] = 278.6400;  m->data[7][6] = 353.3000;
  m->data[8][0] = 16.0000;  m->data[8][1] = 10.0000;  m->data[8][2] = 0.0000;  m->data[8][3] = 202.2550;  m->data[8][4] = 1.2710;  m->data[8][5] = 429.1500;  m->data[8][6] = 666.6500;
  m->data[9][0] = 6.0000;  m->data[9][1] = 12.0000;  m->data[9][2] = 0.0000;  m->data[9][3] = 84.1600;  m->data[9][4] = 0.7800;  m->data[9][5] = 279.0000;  m->data[9][6] = 354.0000;
  m->data[10][0] = 4.0000;  m->data[10][1] = 8.0000;  m->data[10][2] = 1.0000;  m->data[10][3] = 72.1100;  m->data[10][4] = 0.8900;  m->data[10][5] = 164.5000;  m->data[10][6] = 339.0000;
  m->data[11][0] = 4.0000;  m->data[11][1] = 9.0000;  m->data[11][2] = 1.0000;  m->data[11][3] = 71.1100;  m->data[11][4] = 0.8660;  m->data[11][5] = 210.0000;  m->data[11][6] = 360.0000;
  m->data[12][0] = 5.0000;  m->data[12][1] = 11.0000;  m->data[12][2] = 1.0000;  m->data[12][3] = 85.1500;  m->data[12][4] = 0.8620;  m->data[12][5] = 266.0000;  m->data[12][6] = 379.0000;
  m->data[13][0] = 5.0000;  m->data[13][1] = 10.0000;  m->data[13][2] = 1.0000;  m->data[13][3] = 86.1300;  m->data[13][4] = 0.8800;  m->data[13][5] = 228.0000;  m->data[13][6] = 361.0000;

  initDVector(&R);

  PCARankValidation(m, 5, 1, 3, 20, R, &run);
  printf("OK.\n");
  puts("-----------------");
  puts("PC\t R^2");
  for(i = 0; i < R->size; i++)
    printf("%zu\t%f\n", i+1, getDVectorValue(R, i));
  puts("-----------------");

  DelDVector(&R);
  DelMatrix(&m);
}

void test1()
{
  printf("Test 1 - PCA on a simple dataset = ");
  matrix *m; /* Data matrix */
  PCAMODEL *model;
  int run = SIGSCIENTIFICRUN;

  NewMatrix(&m, 14, 7);


  m->data[0][0] = 4.0000;  m->data[0][1] = 4.0000;  m->data[0][2] = 1.0000;  m->data[0][3] = 84.1400;  m->data[0][4] = 1.0500;  m->data[0][5] = 235.1500;  m->data[0][6] = 357.1500;
  m->data[1][0] = 5.0000;  m->data[1][1] = 5.0000;  m->data[1][2] = 1.0000;  m->data[1][3] = 79.1000;  m->data[1][4] = 0.9780;  m->data[1][5] = 1.5090;  m->data[1][6] = 231.0000;
  m->data[2][0] = 4.0000;  m->data[2][1] = 5.0000;  m->data[2][2] = 1.0000;  m->data[2][3] = 67.0900;  m->data[2][4] = 0.9700;  m->data[2][5] = 249.0000;  m->data[2][6] = 403.0000;
  m->data[3][0] = 4.0000;  m->data[3][1] = 4.0000;  m->data[3][2] = 1.0000;  m->data[3][3] = 68.0700;  m->data[3][4] = 0.9360;  m->data[3][5] = 187.3500;  m->data[3][6] = 304.5500;
  m->data[4][0] = 3.0000;  m->data[4][1] = 4.0000;  m->data[4][2] = 2.0000;  m->data[4][3] = 68.0800;  m->data[4][4] = 1.0300;  m->data[4][5] = 363.0000;  m->data[4][6] = 529.0000;
  m->data[5][0] = 9.0000;  m->data[5][1] = 7.0000;  m->data[5][2] = 1.0000;  m->data[5][3] = 129.1600;  m->data[5][4] = 1.0900;  m->data[5][5] = 258.0000;  m->data[5][6] = 510.0000;
  m->data[6][0] = 10.0000;  m->data[6][1] = 8.0000;  m->data[6][2] = 0.0000;  m->data[6][3] = 128.1600;  m->data[6][4] = 1.1500;  m->data[6][5] = 352.0000;  m->data[6][6] = 491.0000;
  m->data[7][0] = 6.0000;  m->data[7][1] = 6.0000;  m->data[7][2] = 0.0000;  m->data[7][3] = 78.1118;  m->data[7][4] = 0.8765;  m->data[7][5] = 278.6400;  m->data[7][6] = 353.3000;
  m->data[8][0] = 16.0000;  m->data[8][1] = 10.0000;  m->data[8][2] = 0.0000;  m->data[8][3] = 202.2550;  m->data[8][4] = 1.2710;  m->data[8][5] = 429.1500;  m->data[8][6] = 666.6500;
  m->data[9][0] = 6.0000;  m->data[9][1] = 12.0000;  m->data[9][2] = 0.0000;  m->data[9][3] = 84.1600;  m->data[9][4] = 0.7800;  m->data[9][5] = 279.0000;  m->data[9][6] = 354.0000;
  m->data[10][0] = 4.0000;  m->data[10][1] = 8.0000;  m->data[10][2] = 1.0000;  m->data[10][3] = 72.1100;  m->data[10][4] = 0.8900;  m->data[10][5] = 164.5000;  m->data[10][6] = 339.0000;
  m->data[11][0] = 4.0000;  m->data[11][1] = 9.0000;  m->data[11][2] = 1.0000;  m->data[11][3] = 71.1100;  m->data[11][4] = 0.8660;  m->data[11][5] = 210.0000;  m->data[11][6] = 360.0000;
  m->data[12][0] = 5.0000;  m->data[12][1] = 11.0000;  m->data[12][2] = 1.0000;  m->data[12][3] = 85.1500;  m->data[12][4] = 0.8620;  m->data[12][5] = 266.0000;  m->data[12][6] = 379.0000;
  m->data[13][0] = 5.0000;  m->data[13][1] = 10.0000;  m->data[13][2] = 1.0000;  m->data[13][3] = 86.1300;  m->data[13][4] = 0.8800;  m->data[13][5] = 228.0000;  m->data[13][6] = 361.0000;


  NewPCAModel(&model);

  PCA(m, 1, 5, model, &run);
  printf("OK.\n");
  PrintPCA(model);

  DelPCAModel(&model);
  DelMatrix(&m);
}

/*
 * Regression test: every component must get its own NIPALS iteration
 * budget (max_iter), and a component that exhausts it must still be
 * stored and deflated.
 *
 * X = C * diag(s) * R where C has orthogonal zero-mean cosine columns and
 * R is a Householder reflection, so the exact PCA is known:
 * loadings = rows of R, explained variance proportional to s^2.
 * s_k = 0.99^k gives small gaps between consecutive eigenvalues: each
 * component needs ~20-200 NIPALS iterations, well within the per-component
 * budget, but more than 1000 in total. When the iteration counter was not
 * reset between components (bug), from PC12 on every component ran a single
 * iteration, was never stored and came back as zero scores/loadings
 * with 0% explained variance.
 */
void test12()
{
  printf("Test PCA 12 - Per-component NIPALS iteration budget  = ");
  enum { NROW = 128, NCOL = 20 };
  const double pi = 3.14159265358979323846;
  double sv[NCOL], R[NCOL][NCOL];
  double vv = 0.0, ss = 0.0, totvar = 0.0;
  size_t i, j, k, a, b;
  matrix *m, *rm;
  PCAMODEL *model;

  for(k = 0; k < NCOL; k++){
    sv[k] = pow(0.99, (double)k);
    ss += sv[k]*sv[k];
  }

  /* Householder reflection R = I - 2 v v' / v'v with v = (1, 2, ..., NCOL) */
  for(i = 0; i < NCOL; i++)
    vv += (double)((i+1)*(i+1));
  for(i = 0; i < NCOL; i++)
    for(j = 0; j < NCOL; j++)
      R[i][j] = (i == j ? 1.0 : 0.0) - 2.0*(double)((i+1)*(j+1))/vv;

  NewMatrix(&m, NROW, NCOL);
  for(i = 0; i < NROW; i++){
    for(j = 0; j < NCOL; j++){
      double x = 0.0;
      for(k = 0; k < NCOL; k++)
        x += cos(2.0*pi*(double)((k+1)*i)/(double)NROW) * sv[k] * R[k][j];
      m->data[i][j] = x;
    }
  }

  NewPCAModel(&model);
  PCA(m, 0, NCOL, model, NULL);

  for(a = 0; a < NCOL; a++){
    /* 1. explained variance of every component matches the exact one */
    double expected = sv[a]*sv[a]/ss*100.0;
    if(!FLOAT_EQ(model->varexp->data[a], expected, 1e-3)){
      printf("FAIL: PC%zu explains %g%% of variance, expected %g%%\n",
             a+1, model->varexp->data[a], expected);
      abort();
    }
    totvar += model->varexp->data[a];

    /* 2. loadings are orthonormal */
    for(b = a; b < NCOL; b++){
      double d = 0.0;
      for(j = 0; j < NCOL; j++)
        d += model->loadings->data[j][a]*model->loadings->data[j][b];
      if(!FLOAT_EQ(d, (a == b ? 1.0 : 0.0), 1e-6)){
        printf("FAIL: loadings PC%zu.PC%zu = %g\n", a+1, b+1, d);
        abort();
      }
    }

    /* 3. loadings match the exact solution (up to the sign) */
    double c = 0.0;
    for(j = 0; j < NCOL; j++)
      c += model->loadings->data[j][a]*R[a][j];
    if(!FLOAT_EQ(fabs(c), 1.0, 1e-4)){
      printf("FAIL: PC%zu |cos| with exact loading = %g\n", a+1, fabs(c));
      abort();
    }
  }

  /* 4. all the components together explain the whole matrix */
  if(!FLOAT_EQ(totvar, 100.0, 1e-3)){
    printf("FAIL: total explained variance %g%%\n", totvar);
    abort();
  }

  /* 5. scores x loadings' reconstruct the input */
  initMatrix(&rm);
  PCAIndVarPredictor(model->scores, model->loadings, model->colaverage, model->colscaling, NCOL, rm);
  for(i = 0; i < NROW; i++){
    for(j = 0; j < NCOL; j++){
      if(!FLOAT_EQ(m->data[i][j], rm->data[i][j], 1e-6)){
        printf("FAIL: reconstruction (%zu,%zu) %g != %g\n", i, j, rm->data[i][j], m->data[i][j]);
        abort();
      }
    }
  }

  printf("OK.\n");
  DelMatrix(&rm);
  DelPCAModel(&model);
  DelMatrix(&m);
}

int main(void)
{
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
  test8();
  test9();
  test10();
  test11();
  test12();
  return 0;
}

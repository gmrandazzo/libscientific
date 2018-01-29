/* testepls.c
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

#include "epls.h"
#include "modelvalidation.h"
#include "matrix.h"
#include "vector.h"
#include "numeric.h"
#include <time.h>

void TestEPLS3()
{
  matrix *x, *y, *y_validation, *y_validation_residuals; /* Data matrix */
  EPLSMODEL *m;

  NewMatrix(&x, 14, 6);
  NewMatrix(&y, 14, 1);

  x->data[0][0] = 4.0000;  x->data[0][1] = 4.0000;  x->data[0][2] = 1.0000;  x->data[0][3] = 84.1400;  x->data[0][4] = 1.0500;  x->data[0][5] = 235.1500;
  x->data[1][0] = 5.0000;  x->data[1][1] = 5.0000;  x->data[1][2] = 1.0000;  x->data[1][3] = 79.1000;  x->data[1][4] = 0.9780;  x->data[1][5] = 231;
  x->data[2][0] = 4.0000;  x->data[2][1] = 5.0000;  x->data[2][2] = 1.0000;  x->data[2][3] = 67.0900;  x->data[2][4] = 0.9700;  x->data[2][5] = 249.0000;
  x->data[3][0] = 4.0000;  x->data[3][1] = 4.0000;  x->data[3][2] = 1.0000;  x->data[3][3] = 68.0700;  x->data[3][4] = 0.9360;  x->data[3][5] = 187.3500;
  x->data[4][0] = 3.0000;  x->data[4][1] = 4.0000;  x->data[4][2] = 2.0000;  x->data[4][3] = 68.0800;  x->data[4][4] = 1.0300;  x->data[4][5] = 363.0000;
  x->data[5][0] = 9.0000;  x->data[5][1] = 7.0000;  x->data[5][2] = 1.0000;  x->data[5][3] = 129.1600;  x->data[5][4] = 1.0900;  x->data[5][5] = 258.0000;
  x->data[6][0] = 10.0000;  x->data[6][1] = 8.0000;  x->data[6][2] = 0.0000;  x->data[6][3] = 128.1600;  x->data[6][4] = 1.1500;  x->data[6][5] = 352.0000;
  x->data[7][0] = 6.0000;  x->data[7][1] = 6.0000;  x->data[7][2] = 0.0000;  x->data[7][3] = 78.1118;  x->data[7][4] = 0.8765;  x->data[7][5] = 278.6400;
  x->data[8][0] = 16.0000;  x->data[8][1] = 10.0000;  x->data[8][2] = 0.0000;  x->data[8][3] = 202.2550;  x->data[8][4] = 1.2710;  x->data[8][5] = 429.1500;
  x->data[9][0] = 6.0000;  x->data[9][1] = 12.0000;  x->data[9][2] = 0.0000;  x->data[9][3] = 84.1600;  x->data[9][4] = 0.7800;  x->data[9][5] = 279.0000;
  x->data[10][0] = 4.0000;  x->data[10][1] = 8.0000;  x->data[10][2] = 1.0000;  x->data[10][3] = 72.1100;  x->data[10][4] = 0.8900;  x->data[10][5] = 164.5000;
  x->data[11][0] = 4.0000;  x->data[11][1] = 9.0000;  x->data[11][2] = 1.0000;  x->data[11][3] = 71.1100;  x->data[11][4] = 0.8660;  x->data[11][5] = 210.0000;
  x->data[12][0] = 5.0000;  x->data[12][1] = 11.0000;  x->data[12][2] = 1.0000;  x->data[12][3] = 85.1500;  x->data[12][4] = 0.8620;  x->data[12][5] = 266.0000;
  x->data[13][0] = 5.0000;  x->data[13][1] = 10.0000;  x->data[13][2] = 1.0000;  x->data[13][3] = 86.1300;  x->data[13][4] = 0.8800;  x->data[13][5] = 228.0000;


  y->data[0][0] = 357.1500;
  y->data[1][0] = 388.0000;
  y->data[2][0] = 403.0000;
  y->data[3][0] = 304.5500;
  y->data[4][0] = 529.0000;
  y->data[5][0] = 510.0000;
  y->data[6][0] = 491.0000;
  y->data[7][0] = 353.3000;
  y->data[8][0] = 666.6500;
  y->data[9][0] = 354.0000;
  y->data[10][0] = 339.0000;
  y->data[11][0] = 360.0000;
  y->data[12][0] = 379.0000;
  y->data[13][0] = 361.0000;

  printf("Test EPLS 3: Dynamic Random Subspace Method\n");

  /*Allocate the final output*/
  NewEPLSModel(&m);

  ELearningParameters eparm;
  eparm.algorithm = BaggingRandomSubspaceMethod;
  eparm.n_models = 1000;
  eparm.r_fix = 4;
  eparm.trainsize = 0.7;
  size_t nlv = 5; /*This will ignored and set to max 4*/
  EPLS(x, y, nlv, 1, 0, m, eparm, NULL);

  //PrintEPLSModel(m);

  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = 4;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;
  size_t nthreads = 1;
  initMatrix(&y_validation);
  initMatrix(&y_validation_residuals);
  BootstrapRandomGroupsCV(&minpt, 5, 20, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Median);
  //LeaveOneOut(&minpt, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Median);
  //PrintMatrix(y);
  //PrintMatrix(y_validation);
  matrix *q2;
  matrix *sdep;
  matrix *bias;
  initMatrix(&q2);
  initMatrix(&sdep);
  initMatrix(&bias);
  EPLSRegressionStatistics(y, y_validation, &q2, &sdep, &bias);
  puts("Q2 Dynamic Random Subspace Method EPLS");
  PrintMatrix(q2);
  puts("Standard Deviation Prediction Error Dynamic Random Subspace Method EPLS");
  PrintMatrix(sdep);
  puts("BIAS Dynamic Random Subspace Method EPLS");
  PrintMatrix(bias);
  DelMatrix(&bias);
  DelMatrix(&sdep);
  DelMatrix(&q2);
  DelMatrix(&y_validation);
  DelMatrix(&y_validation_residuals);
  DelEPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);
}

void TestEPLS2()
{
  matrix *x, *y, *y_validation, *y_validation_residuals; /* Data matrix */
  EPLSMODEL *m;

  NewMatrix(&x, 14, 6);
  NewMatrix(&y, 14, 1);

  x->data[0][0] = 4.0000;  x->data[0][1] = 4.0000;  x->data[0][2] = 1.0000;  x->data[0][3] = 84.1400;  x->data[0][4] = 1.0500;  x->data[0][5] = 235.1500;
  x->data[1][0] = 5.0000;  x->data[1][1] = 5.0000;  x->data[1][2] = 1.0000;  x->data[1][3] = 79.1000;  x->data[1][4] = 0.9780;  x->data[1][5] = 231;
  x->data[2][0] = 4.0000;  x->data[2][1] = 5.0000;  x->data[2][2] = 1.0000;  x->data[2][3] = 67.0900;  x->data[2][4] = 0.9700;  x->data[2][5] = 249.0000;
  x->data[3][0] = 4.0000;  x->data[3][1] = 4.0000;  x->data[3][2] = 1.0000;  x->data[3][3] = 68.0700;  x->data[3][4] = 0.9360;  x->data[3][5] = 187.3500;
  x->data[4][0] = 3.0000;  x->data[4][1] = 4.0000;  x->data[4][2] = 2.0000;  x->data[4][3] = 68.0800;  x->data[4][4] = 1.0300;  x->data[4][5] = 363.0000;
  x->data[5][0] = 9.0000;  x->data[5][1] = 7.0000;  x->data[5][2] = 1.0000;  x->data[5][3] = 129.1600;  x->data[5][4] = 1.0900;  x->data[5][5] = 258.0000;
  x->data[6][0] = 10.0000;  x->data[6][1] = 8.0000;  x->data[6][2] = 0.0000;  x->data[6][3] = 128.1600;  x->data[6][4] = 1.1500;  x->data[6][5] = 352.0000;
  x->data[7][0] = 6.0000;  x->data[7][1] = 6.0000;  x->data[7][2] = 0.0000;  x->data[7][3] = 78.1118;  x->data[7][4] = 0.8765;  x->data[7][5] = 278.6400;
  x->data[8][0] = 16.0000;  x->data[8][1] = 10.0000;  x->data[8][2] = 0.0000;  x->data[8][3] = 202.2550;  x->data[8][4] = 1.2710;  x->data[8][5] = 429.1500;
  x->data[9][0] = 6.0000;  x->data[9][1] = 12.0000;  x->data[9][2] = 0.0000;  x->data[9][3] = 84.1600;  x->data[9][4] = 0.7800;  x->data[9][5] = 279.0000;
  x->data[10][0] = 4.0000;  x->data[10][1] = 8.0000;  x->data[10][2] = 1.0000;  x->data[10][3] = 72.1100;  x->data[10][4] = 0.8900;  x->data[10][5] = 164.5000;
  x->data[11][0] = 4.0000;  x->data[11][1] = 9.0000;  x->data[11][2] = 1.0000;  x->data[11][3] = 71.1100;  x->data[11][4] = 0.8660;  x->data[11][5] = 210.0000;
  x->data[12][0] = 5.0000;  x->data[12][1] = 11.0000;  x->data[12][2] = 1.0000;  x->data[12][3] = 85.1500;  x->data[12][4] = 0.8620;  x->data[12][5] = 266.0000;
  x->data[13][0] = 5.0000;  x->data[13][1] = 10.0000;  x->data[13][2] = 1.0000;  x->data[13][3] = 86.1300;  x->data[13][4] = 0.8800;  x->data[13][5] = 228.0000;


  y->data[0][0] = 357.1500;
  y->data[1][0] = 388.0000;
  y->data[2][0] = 403.0000;
  y->data[3][0] = 304.5500;
  y->data[4][0] = 529.0000;
  y->data[5][0] = 510.0000;
  y->data[6][0] = 491.0000;
  y->data[7][0] = 353.3000;
  y->data[8][0] = 666.6500;
  y->data[9][0] = 354.0000;
  y->data[10][0] = 339.0000;
  y->data[11][0] = 360.0000;
  y->data[12][0] = 379.0000;
  y->data[13][0] = 361.0000;

  printf("Test EPLS 2: Fixed Random Subspace Method\n");

  /*Allocate the final output*/
  NewEPLSModel(&m);

  ELearningParameters eparm;
  eparm.algorithm = FixedRandomSubspaceMethod;
  eparm.n_models = 1000;
  eparm.r_fix = 3;
  size_t nlv = 6;
  EPLS(x, y, nlv, 1, 0, m, eparm, NULL);

  //PrintEPLSModel(m);

  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = nlv;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;
  size_t nthreads = 4;
  initMatrix(&y_validation);
  initMatrix(&y_validation_residuals);
  //BootstrapRandomGroupsCV(&minpt, 5, 20, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  LeaveOneOut(&minpt, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  //PrintMatrix(y);
  //PrintMatrix(y_validation);
  matrix *q2;
  matrix *sdep;
  matrix *bias;
  initMatrix(&q2);
  initMatrix(&sdep);
  initMatrix(&bias);
  EPLSRegressionStatistics(y, y_validation, &q2, &sdep, &bias);
  puts("Q2 Fixed Random Subspace Method EPLS");
  PrintMatrix(q2);
  puts("Standard Deviation Prediction Error Fixed Random Subspace Method EPLS");
  PrintMatrix(sdep);
  puts("BIAS Fixed Random Subspace Method EPLS");
  PrintMatrix(bias);
  DelMatrix(&bias);
  DelMatrix(&sdep);
  DelMatrix(&q2);
  DelMatrix(&y_validation);
  DelMatrix(&y_validation_residuals);
  DelEPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);
}

void TestEPLS1()
{
  matrix *x, *y, *y_validation, *y_validation_residuals; /* Data matrix */
  EPLSMODEL *m;

  NewMatrix(&x, 14, 6);
  NewMatrix(&y, 14, 1);

  x->data[0][0] = 4.0000;  x->data[0][1] = 4.0000;  x->data[0][2] = 1.0000;  x->data[0][3] = 84.1400;  x->data[0][4] = 1.0500;  x->data[0][5] = 235.1500;
  x->data[1][0] = 5.0000;  x->data[1][1] = 5.0000;  x->data[1][2] = 1.0000;  x->data[1][3] = 79.1000;  x->data[1][4] = 0.9780;  x->data[1][5] = 231;
  x->data[2][0] = 4.0000;  x->data[2][1] = 5.0000;  x->data[2][2] = 1.0000;  x->data[2][3] = 67.0900;  x->data[2][4] = 0.9700;  x->data[2][5] = 249.0000;
  x->data[3][0] = 4.0000;  x->data[3][1] = 4.0000;  x->data[3][2] = 1.0000;  x->data[3][3] = 68.0700;  x->data[3][4] = 0.9360;  x->data[3][5] = 187.3500;
  x->data[4][0] = 3.0000;  x->data[4][1] = 4.0000;  x->data[4][2] = 2.0000;  x->data[4][3] = 68.0800;  x->data[4][4] = 1.0300;  x->data[4][5] = 363.0000;
  x->data[5][0] = 9.0000;  x->data[5][1] = 7.0000;  x->data[5][2] = 1.0000;  x->data[5][3] = 129.1600;  x->data[5][4] = 1.0900;  x->data[5][5] = 258.0000;
  x->data[6][0] = 10.0000;  x->data[6][1] = 8.0000;  x->data[6][2] = 0.0000;  x->data[6][3] = 128.1600;  x->data[6][4] = 1.1500;  x->data[6][5] = 352.0000;
  x->data[7][0] = 6.0000;  x->data[7][1] = 6.0000;  x->data[7][2] = 0.0000;  x->data[7][3] = 78.1118;  x->data[7][4] = 0.8765;  x->data[7][5] = 278.6400;
  x->data[8][0] = 16.0000;  x->data[8][1] = 10.0000;  x->data[8][2] = 0.0000;  x->data[8][3] = 202.2550;  x->data[8][4] = 1.2710;  x->data[8][5] = 429.1500;
  x->data[9][0] = 6.0000;  x->data[9][1] = 12.0000;  x->data[9][2] = 0.0000;  x->data[9][3] = 84.1600;  x->data[9][4] = 0.7800;  x->data[9][5] = 279.0000;
  x->data[10][0] = 4.0000;  x->data[10][1] = 8.0000;  x->data[10][2] = 1.0000;  x->data[10][3] = 72.1100;  x->data[10][4] = 0.8900;  x->data[10][5] = 164.5000;
  x->data[11][0] = 4.0000;  x->data[11][1] = 9.0000;  x->data[11][2] = 1.0000;  x->data[11][3] = 71.1100;  x->data[11][4] = 0.8660;  x->data[11][5] = 210.0000;
  x->data[12][0] = 5.0000;  x->data[12][1] = 11.0000;  x->data[12][2] = 1.0000;  x->data[12][3] = 85.1500;  x->data[12][4] = 0.8620;  x->data[12][5] = 266.0000;
  x->data[13][0] = 5.0000;  x->data[13][1] = 10.0000;  x->data[13][2] = 1.0000;  x->data[13][3] = 86.1300;  x->data[13][4] = 0.8800;  x->data[13][5] = 228.0000;


  y->data[0][0] = 357.1500;
  y->data[1][0] = 388.0000;
  y->data[2][0] = 403.0000;
  y->data[3][0] = 304.5500;
  y->data[4][0] = 529.0000;
  y->data[5][0] = 510.0000;
  y->data[6][0] = 491.0000;
  y->data[7][0] = 353.3000;
  y->data[8][0] = 666.6500;
  y->data[9][0] = 354.0000;
  y->data[10][0] = 339.0000;
  y->data[11][0] = 360.0000;
  y->data[12][0] = 379.0000;
  y->data[13][0] = 361.0000;

  printf("Test EPLS 1: Bagging PLS\n");

  /*Allocate the final output*/
  NewEPLSModel(&m);

  ELearningParameters eparm;
  eparm.algorithm = Bagging;
  eparm.n_models = 2;
  eparm.trainsize = 0.7;
  size_t nlv = 5;
  EPLS(x, y, nlv, 1, 0, m, eparm, NULL);

  //PrintEPLSModel(m);

  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = nlv;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;
  size_t nthreads = 4;
  initMatrix(&y_validation);
  initMatrix(&y_validation_residuals);
  //BootstrapRandomGroupsCV(&minpt, 5, 20, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  LeaveOneOut(&minpt, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  /*PrintMatrix(y);
  PrintMatrix(y_validation);*/
  matrix *q2;
  matrix *sdep;
  matrix *bias;
  initMatrix(&q2);
  initMatrix(&sdep);
  initMatrix(&bias);
  EPLSRegressionStatistics(y, y_validation, &q2, &sdep, &bias);
  puts("Q2 Bagging EPLS");
  PrintMatrix(q2);
  puts("Standard Deviation Prediction Error Bagging EPLS");
  PrintMatrix(sdep);
  puts("BIAS Bagging EPLS");
  PrintMatrix(bias);
  DelMatrix(&bias);
  DelMatrix(&sdep);
  DelMatrix(&q2);
  DelMatrix(&y_validation);
  DelMatrix(&y_validation_residuals);
  DelEPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);
}

int main(void)
{
  /*test 1- 5*/
  //TestEPLS1();
  //TestEPLS2();
  TestEPLS3();
}

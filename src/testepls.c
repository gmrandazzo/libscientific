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

void TestEPLS5()
{
  matrix *x, *y, *y_validation, *y_validation_residuals; /* Data matrix */
  uivector *groups;
  EPLSMODEL *m;

  NewMatrix(&x, 14, 6);
  NewMatrix(&y, 14, 1);
  NewUIVector(&groups, 14);

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

  groups->data[0] = 1;
  groups->data[1] = 1;
  groups->data[2] = 0;
  groups->data[3] = 1;
  groups->data[4] = 0;
  groups->data[5] = 0;
  groups->data[6] = 0;
  groups->data[7] = 1;
  groups->data[8] = 0;
  groups->data[9] = 1;
  groups->data[10] = 1;
  groups->data[11] = 1;
  groups->data[12] = 1;
  groups->data[13] = 1;

  printf("Test EPLS 5: Dynamic Random Subspace Method and KFold Cross Validation\n");

  /*Allocate the final output*/
  NewEPLSModel(&m);

  ELearningParameters eparm;
  eparm.algorithm = BaggingRandomSubspaceMethod;
  eparm.n_models = 1000;
  eparm.r_fix = 4;
  eparm.trainsize = 0.7;
  size_t nlv = 4; /*This will ignored and set to max 4*/
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
  //BootstrapRandomGroupsCV(&minpt, 5, 20, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  //LeaveOneOut(&minpt, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  KFoldCV(&minpt, groups, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  PrintMatrix(y);
  PrintMatrix(y_validation);
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
  DelUIVector(&groups);
  DelMatrix(&x);
  DelMatrix(&y);
}

void TestEPLS4()
{
  matrix *x, *y, *y_validation, *y_validation_residuals; /* Data matrix */
  matrix *x_test, *y_test, *y_test_predicted;
  EPLSMODEL *m;

  NewMatrix(&x, 60, 2);
  NewMatrix(&y, 60, 1);

  x->data[0][0] = 0.309205; x->data[0][1] = -0.127041; y->data[0][0] = 0.000000;
  x->data[1][0] = -0.486914; x->data[1][1] = 0.572599; y->data[1][0] = 1.000000;
  x->data[2][0] = 0.509549; x->data[2][1] = -1.664272; y->data[2][0] = 1.000000;
  x->data[3][0] = 1.352429; x->data[3][1] = 0.608595; y->data[3][0] = 1.000000;
  x->data[4][0] = 1.896001; x->data[4][1] = -0.426535; y->data[4][0] = 1.000000;
  x->data[5][0] = 0.747795; x->data[5][1] = 0.549574; y->data[5][0] = 0.000000;
  x->data[6][0] = -1.198059; x->data[6][1] = 1.400476; y->data[6][0] = 0.000000;
  x->data[7][0] = -1.406588; x->data[7][1] = -1.036872; y->data[7][0] = 0.000000;
  x->data[8][0] = -0.330754; x->data[8][1] = -0.553255; y->data[8][0] = 1.000000;
  x->data[9][0] = 0.817738; x->data[9][1] = 1.643003; y->data[9][0] = 0.000000;
  x->data[10][0] = -1.306676; x->data[10][1] = -0.607112; y->data[10][0] = 0.000000;
  x->data[11][0] = -0.626639; x->data[11][1] = -0.450971; y->data[11][0] = 1.000000;
  x->data[12][0] = -0.221201; x->data[12][1] = 1.008688; y->data[12][0] = 0.000000;
  x->data[13][0] = -0.851516; x->data[13][1] = -0.060227; y->data[13][0] = 1.000000;
  x->data[14][0] = -0.597399; x->data[14][1] = 0.278299; y->data[14][0] = 1.000000;
  x->data[15][0] = 0.155818; x->data[15][1] = -1.145737; y->data[15][0] = 1.000000;
  x->data[16][0] = 0.078135; x->data[16][1] = 0.184332; y->data[16][0] = 0.000000;
  x->data[17][0] = 0.241302; x->data[17][1] = -1.743719; y->data[17][0] = 1.000000;
  x->data[18][0] = -0.424436; x->data[18][1] = 1.308545; y->data[18][0] = 0.000000;
  x->data[19][0] = 0.008145; x->data[19][1] = -0.586563; y->data[19][0] = 1.000000;
  x->data[20][0] = 0.031382; x->data[20][1] = 1.278677; y->data[20][0] = 0.000000;
  x->data[21][0] = -2.027163; x->data[21][1] = -0.075527; y->data[21][0] = 0.000000;
  x->data[22][0] = -0.202164; x->data[22][1] = 0.098295; y->data[22][0] = 0.000000;
  x->data[23][0] = -1.770933; x->data[23][1] = -0.156025; y->data[23][0] = 0.000000;
  x->data[24][0] = 1.019755; x->data[24][1] = -0.399649; y->data[24][0] = 1.000000;
  x->data[25][0] = 1.224674; x->data[25][1] = -1.455377; y->data[25][0] = 1.000000;
  x->data[26][0] = 0.198951; x->data[26][1] = 2.276942; y->data[26][0] = 0.000000;
  x->data[27][0] = 0.275245; x->data[27][1] = 0.589769; y->data[27][0] = 0.000000;
  x->data[28][0] = -0.984345; x->data[28][1] = 1.295465; y->data[28][0] = 0.000000;
  x->data[29][0] = 0.510383; x->data[29][1] = -1.099258; y->data[29][0] = 1.000000;
  x->data[30][0] = 0.887964; x->data[30][1] = -1.220074; y->data[30][0] = 1.000000;
  x->data[31][0] = 0.159400; x->data[31][1] = -1.055063; y->data[31][0] = 1.000000;
  x->data[32][0] = -0.312520; x->data[32][1] = 1.853613; y->data[32][0] = 0.000000;
  x->data[33][0] = -1.720182; x->data[33][1] = -0.113235; y->data[33][0] = 0.000000;
  x->data[34][0] = -2.051619; x->data[34][1] = 0.329559; y->data[34][0] = 0.000000;
  x->data[35][0] = -0.814450; x->data[35][1] = 0.038820; y->data[35][0] = 0.000000;
  x->data[36][0] = -0.163612; x->data[36][1] = -1.537821; y->data[36][0] = 1.000000;
  x->data[37][0] = 1.539163; x->data[37][1] = -0.406466; y->data[37][0] = 1.000000;
  x->data[38][0] = 0.755665; x->data[38][1] = 0.626711; y->data[38][0] = 0.000000;
  x->data[39][0] = -0.685180; x->data[39][1] = 0.564629; y->data[39][0] = 1.000000;
  x->data[40][0] = -0.420818; x->data[40][1] = -1.031264; y->data[40][0] = 1.000000;
  x->data[41][0] = -1.611818; x->data[41][1] = -1.161224; y->data[41][0] = 0.000000;
  x->data[42][0] = -0.360045; x->data[42][1] = 0.031320; y->data[42][0] = 1.000000;
  x->data[43][0] = 1.108777; x->data[43][1] = -0.565967; y->data[43][0] = 1.000000;
  x->data[44][0] = 1.821248; x->data[44][1] = -0.466077; y->data[44][0] = 1.000000;
  x->data[45][0] = 0.249423; x->data[45][1] = 0.370420; y->data[45][0] = 0.000000;
  x->data[46][0] = -0.041567; x->data[46][1] = -0.718505; y->data[46][0] = 1.000000;
  x->data[47][0] = 0.430949; x->data[47][1] = -0.516105; y->data[47][0] = 0.000000;
  x->data[48][0] = 0.230930; x->data[48][1] = -1.165639; y->data[48][0] = 1.000000;
  x->data[49][0] = 0.262854; x->data[49][1] = 0.089705; y->data[49][0] = 0.000000;
  x->data[50][0] = -1.491288; x->data[50][1] = 0.682659; y->data[50][0] = 0.000000;
  x->data[51][0] = 1.855848; x->data[51][1] = 0.397113; y->data[51][0] = 1.000000;
  x->data[52][0] = -0.341948; x->data[52][1] = -0.502986; y->data[52][0] = 1.000000;
  x->data[53][0] = 0.890636; x->data[53][1] = -0.738440; y->data[53][0] = 1.000000;
  x->data[54][0] = -1.237843; x->data[54][1] = 0.726437; y->data[54][0] = 0.000000;
  x->data[55][0] = 0.153508; x->data[55][1] = 0.891228; y->data[55][0] = 0.000000;
  x->data[56][0] = -1.636016; x->data[56][1] = 1.144093; y->data[56][0] = 0.000000;
  x->data[57][0] = -0.037348; x->data[57][1] = 0.539623; y->data[57][0] = 0.000000;
  x->data[58][0] = -0.224132; x->data[58][1] = 0.944773; y->data[58][0] = 0.000000;
  x->data[59][0] = -0.994714; x->data[59][1] = -0.785859; y->data[59][0] = 1.000000;

  NewMatrix(&x_test, 40, 2);
  NewMatrix(&y_test, 40, 1);
  x_test->data[0][0] = -1.551825; x_test->data[0][1] = 0.469108; y_test->data[0][0] = 0.000000;
  x_test->data[1][0] = -1.143596; x_test->data[1][1] = 0.762807; y_test->data[1][0] = 0.000000;
  x_test->data[2][0] = -0.287244; x_test->data[2][1] = 1.052591; y_test->data[2][0] = 0.000000;
  x_test->data[3][0] = 0.486019; x_test->data[3][1] = 1.960419; y_test->data[3][0] = 0.000000;
  x_test->data[4][0] = 0.593482; x_test->data[4][1] = 0.252789; y_test->data[4][0] = 0.000000;
  x_test->data[5][0] = 0.592601; x_test->data[5][1] = -1.448428; y_test->data[5][0] = 1.000000;
  x_test->data[6][0] = 0.207643; x_test->data[6][1] = 0.888564; y_test->data[6][0] = 0.000000;
  x_test->data[7][0] = -0.719785; x_test->data[7][1] = 1.723347; y_test->data[7][0] = 0.000000;
  x_test->data[8][0] = -0.878941; x_test->data[8][1] = 0.653493; y_test->data[8][0] = 1.000000;
  x_test->data[9][0] = -0.528803; x_test->data[9][1] = 1.285421; y_test->data[9][0] = 0.000000;
  x_test->data[10][0] = -0.831901; x_test->data[10][1] = 0.932441; y_test->data[10][0] = 0.000000;
  x_test->data[11][0] = -0.029746; x_test->data[11][1] = -0.634758; y_test->data[11][0] = 1.000000;
  x_test->data[12][0] = 1.717335; x_test->data[12][1] = -0.178755; y_test->data[12][0] = 1.000000;
  x_test->data[13][0] = -0.091303; x_test->data[13][1] = -0.938243; y_test->data[13][0] = 1.000000;
  x_test->data[14][0] = 1.316503; x_test->data[14][1] = -1.107888; y_test->data[14][0] = 1.000000;
  x_test->data[15][0] = 0.371231; x_test->data[15][1] = -1.301925; y_test->data[15][0] = 1.000000;
  x_test->data[16][0] = -0.173442; x_test->data[16][1] = 1.503620; y_test->data[16][0] = 0.000000;
  x_test->data[17][0] = 1.397303; x_test->data[17][1] = 0.390400; y_test->data[17][0] = 1.000000;
  x_test->data[18][0] = 1.127568; x_test->data[18][1] = -1.652896; y_test->data[18][0] = 1.000000;
  x_test->data[19][0] = 0.626822; x_test->data[19][1] = -2.143372; y_test->data[19][0] = 1.000000;
  x_test->data[20][0] = -0.344084; x_test->data[20][1] = -0.115017; y_test->data[20][0] = 1.000000;
  x_test->data[21][0] = 1.333808; x_test->data[21][1] = -0.730181; y_test->data[21][0] = 1.000000;
  x_test->data[22][0] = 0.118921; x_test->data[22][1] = 0.142101; y_test->data[22][0] = 0.000000;
  x_test->data[23][0] = -1.552685; x_test->data[23][1] = -0.308749; y_test->data[23][0] = 0.000000;
  x_test->data[24][0] = 1.894578; x_test->data[24][1] = -1.573924; y_test->data[24][0] = 1.000000;
  x_test->data[25][0] = -0.977863; x_test->data[25][1] = 1.510417; y_test->data[25][0] = 0.000000;
  x_test->data[26][0] = 1.150042; x_test->data[26][1] = -0.464233; y_test->data[26][0] = 0.000000;
  x_test->data[27][0] = -0.317847; x_test->data[27][1] = -1.493817; y_test->data[27][0] = 1.000000;
  x_test->data[28][0] = 1.370862; x_test->data[28][1] = 0.223750; y_test->data[28][0] = 1.000000;
  x_test->data[29][0] = -1.559623; x_test->data[29][1] = -0.264803; y_test->data[29][0] = 0.000000;
  x_test->data[30][0] = 1.695311; x_test->data[30][1] = -0.911746; y_test->data[30][0] = 1.000000;
  x_test->data[31][0] = -0.980486; x_test->data[31][1] = 0.847862; y_test->data[31][0] = 0.000000;
  x_test->data[32][0] = -0.536313; x_test->data[32][1] = -0.205491; y_test->data[32][0] = 1.000000;
  x_test->data[33][0] = 1.600973; x_test->data[33][1] = -0.095316; y_test->data[33][0] = 1.000000;
  x_test->data[34][0] = 0.988402; x_test->data[34][1] = -1.759330; y_test->data[34][0] = 1.000000;
  x_test->data[35][0] = 1.818721; x_test->data[35][1] = 0.734654; y_test->data[35][0] = 1.000000;
  x_test->data[36][0] = -0.731536; x_test->data[36][1] = 1.242361; y_test->data[36][0] = 0.000000;
  x_test->data[37][0] = 0.670509; x_test->data[37][1] = -1.012216; y_test->data[37][0] = 1.000000;
  x_test->data[38][0] = -0.759938; x_test->data[38][1] = 1.311125; y_test->data[38][0] = 0.000000;
  x_test->data[39][0] = -0.224656; x_test->data[39][1] = 1.702718; y_test->data[39][0] = 0.000000;

  printf("Test EPLS 4: Dynamic Random Subspace Method Discriminant Analysis\n");

  /*Allocate the final output*/
  NewEPLSModel(&m);

  ELearningParameters eparm;
  eparm.algorithm = Bagging;
  eparm.n_models = 100;
  eparm.trainsize = 0.3;

  size_t nlv = 5; /*This will ignored and set to max 4*/
  EPLS(x, y, nlv, 1, 0, m, eparm, NULL);

  //PrintEPLSModel(m);

  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = 4;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;
  size_t nthreads = 4;
  initMatrix(&y_validation);
  initMatrix(&y_validation_residuals);
  //BootstrapRandomGroupsCV(&minpt, 5, 20, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  LeaveOneOut(&minpt, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  /*PrintMatrix(y);
  PrintMatrix(y_validation);*/
  tensor *roc;
  matrix *roc_auc;
  tensor *pr;
  matrix *pr_ap;
  initTensor(&roc);
  initMatrix(&roc_auc);
  initTensor(&pr);
  initMatrix(&pr_ap);
  EPLSDiscriminantAnalysisStatistics(y, y_validation, &roc, &roc_auc, &pr, &pr_ap);

  puts("ROC AUCs Boosting EPLS-DA");
  PrintMatrix(roc_auc);
  puts("Precision-Recall Boosting EPLS-DA");
  PrintMatrix(pr_ap);


  initMatrix(&y_test_predicted);
  EPLSYPRedictorAllLV(x_test, m, Averaging, NULL, &y_test_predicted);
  matrix *roc_auc_test, *pr_ap_test;
  initMatrix(&roc_auc_test);
  initMatrix(&pr_ap_test);

  EPLSDiscriminantAnalysisStatistics(y, y_validation, NULL, &roc_auc_test, NULL, &pr_ap_test);

  puts("ROC AUCs External TEST Boosting EPLS-DA");
  PrintMatrix(roc_auc_test);
  puts("Precision-Recall External TEST Boosting EPLS-DA");
  PrintMatrix(pr_ap_test);

  DelMatrix(&pr_ap_test);
  DelMatrix(&roc_auc_test);
  DelMatrix(&y_test_predicted);
  DelMatrix(&roc_auc);
  DelMatrix(&pr_ap);
  DelTensor(&roc);
  DelTensor(&pr);
  DelMatrix(&y_validation);
  DelMatrix(&y_validation_residuals);
  DelEPLSModel(&m);
  DelMatrix(&x_test);
  DelMatrix(&y_test);
  DelMatrix(&x);
  DelMatrix(&y);
}

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
  eparm.n_models = 10000;
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
  size_t nthreads = 4;
  initMatrix(&y_validation);
  initMatrix(&y_validation_residuals);
  BootstrapRandomGroupsCV(&minpt, 5, 20, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  //LeaveOneOut(&minpt, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  PrintMatrix(y);
  PrintMatrix(y_validation);
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
  size_t nthreads = 1;
  initMatrix(&y_validation);
  initMatrix(&y_validation_residuals);
  BootstrapRandomGroupsCV(&minpt, 5, 20, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
  //LeaveOneOut(&minpt, _EPLS_, &y_validation, &y_validation_residuals, nthreads, NULL, 2, eparm, Averaging);
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
  //TestEPLS3();
  //TestEPLS4();
  TestEPLS5();
}

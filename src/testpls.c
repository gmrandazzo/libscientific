/* testpls.c
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

#include "array.h"
#include "pls.h"
#include "modelvalidation.h"
#include "matrix.h"
#include "vector.h"
#include "numeric.h"
#include <time.h>


void TestPLS15()
{
  matrix *x, *y;

  PLSMODEL *m;

  array *roc;
  array *precision_recall;
  matrix *roc_auc;
  matrix *precision_recall_ap;

  matrix *xpred;
  matrix *xpredscores;

  matrix *ypred;

  size_t nobj = 100;
  size_t nobj_to_pred = 50;
  size_t nfeats = 250;
  size_t ntarg = 1;
  NewMatrix(&x, nobj, nfeats);
  NewMatrix(&xpred, nobj_to_pred, nfeats);
  NewMatrix(&y, nobj, ntarg);

  srand(nobj+nfeats+ntarg);
  for(size_t i = 0; i < nobj; i++){
    for(size_t j = 0; j < nfeats; j++){
      setMatrixValue(x, i, j, randDouble(-100,100));
    }

    for(size_t j = 0; j < ntarg; j++){
      setMatrixValue(y, i, j, randInt(0,2));
    }
  }

  for(size_t i = 0; i < nobj_to_pred; i++){
    for(size_t j = 0; j < nfeats; j++){
      setMatrixValue(xpred, i, j, randDouble(-100,100));
    }
  }


  /*NewMatrix(&x, 4, 2);
  NewMatrix(&y, 4, 1);

  setMatrixValue(x, 0, 0, 5); setMatrixValue(x, 0, 1, 7);
  setMatrixValue(x, 1, 0, 3); setMatrixValue(x, 1, 1, 12);
  setMatrixValue(x, 2, 0, 2); setMatrixValue(x, 2, 1, 1);
  setMatrixValue(x, 3, 0, 4); setMatrixValue(x, 3, 1, 8);

  setMatrixValue(y, 0, 0, 0);
  setMatrixValue(y, 1, 0, 0);
  setMatrixValue(y, 2, 0, 1);
  setMatrixValue(y, 3, 0, 1);


  NewMatrix(&xpred, 3, 2);
  setMatrixValue(xpred, 0, 0, 5); setMatrixValue(xpred, 0, 1, 5);  // y = 7,33
  setMatrixValue(xpred, 1, 0, 12); setMatrixValue(xpred, 1, 1, 4); // y = 13,99
  setMatrixValue(xpred, 2, 0, 21); setMatrixValue(xpred, 2, 1, 15); // y = 26,66
  */

  puts("Test 15 - PLS DA Model and prediction");
  /*Allocate the final output*/
  NewPLSModel(&m);

  /*puts("X:");
  PrintMatrix(x);

  puts("Y:");
  PrintMatrix(y);*/

  ssignal run = SIGSCIENTIFICRUN;
  PLS(x, y, 4, 1, 0, m, &run);



  /*VALIDATE THE MODEL */
  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = 4;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, &m->predicted_y, &m->pred_residuals, 4, NULL, 0);
  PrintMatrix(m->predicted_y);

  initArray(&roc);
  initArray(&precision_recall);
  initMatrix(&roc_auc);
  initMatrix(&precision_recall_ap);

  PLSDiscriminantAnalysisStatistics(y, m->predicted_y, &roc, &roc_auc, &precision_recall, &precision_recall_ap);

  PrintPLSModel(m);

  puts("ROC AUC's for Trainig set");
  PrintMatrix(roc_auc);
  puts("ROC Curves for each LV");
  PrintArray(roc);

  puts("Precision-Recall or Trainig set");
  PrintMatrix(precision_recall_ap);
  puts("Precision-Recall Curves for each LV");
  PrintArray(precision_recall);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, &xpredscores);

  PLSYPredictor(xpredscores, m, 1, &ypred);

  /*puts("\nPrediction scores...");
  puts("x");
  PrintMatrix(xpredscores);

  puts("Extimated Y:");
  PrintMatrix(ypred);
  */
  DelMatrix(&roc_auc);
  DelMatrix(&precision_recall_ap);
  DelArray(&roc);
  DelArray(&precision_recall);

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);

  DelMatrix(&xpredscores);

  DelMatrix(&xpred);
  DelMatrix(&ypred);
}

void TestPLS14()
{
  matrix *x, *y;
  PLSMODEL *m;

  NewMatrix(&x, 14, 6);
  NewMatrix(&y, 14, 2);

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

  y->data[0][1] = 0.1500;
  y->data[1][1] = 0.0000;
  y->data[2][1] = 1.0000;
  y->data[3][1] = 0.5500;
  y->data[4][1] = 2.0000;
  y->data[5][1] = 2.0000;
  y->data[6][1] = 1.0000;
  y->data[7][1] = 0.3000;
  y->data[8][1] = 3.6500;
  y->data[9][1] = 0.0000;
  y->data[10][1] = 0.0000;
  y->data[11][1] = 0.0000;
  y->data[12][1] = 0.0000;
  y->data[13][1] = 1.0000;

  printf("Test PLS 14\n");

  NewPLSModel(&m);
  PLS(x, y, 100, 1, 0, m, NULL);

  /*VALIDATE THE MODEL */
  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = 999;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  //BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, &m->predicted_y, &m->pred_residuals, 4, NULL, 0);
  LeaveOneOut(&minpt, _PLS_, &m->predicted_y, &m->pred_residuals, 4, NULL, 0);
  PrintMatrix(m->predicted_y);

  PLSRegressionStatistics(y, m->predicted_y, &m->q2y, &m->sdep, &m->bias);
  PrintMatrix(m->predicted_y);
  puts("Q2 Cross Validation");
  PrintMatrix(m->q2y);
  puts("SDEP Cross Validation");
  PrintMatrix(m->sdep);
  puts("BIAS Cross Validation");
  PrintMatrix(m->bias);

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);
}


void TestPLS13()
{
  size_t nobj = 800;
  matrix *mx, *my;
  matrix *q2, *sdep, *bias;
  matrix *predicted_y;
  matrix *predicted_residuals;
  puts("Test13: Simple Calculation PLS Model with RGCV with random data");


  NewMatrix(&mx, nobj, 128);
  NewMatrix(&my, nobj, 1);

  srand(nobj);
  for(size_t i = 0; i < nobj; i++){
    for(size_t j = 0; j < 128; j++){
      setMatrixValue(mx, i, j, randDouble(0,20));
    }

    setMatrixValue(my, i, 0, randDouble(0,1));
  }

  initMatrix(&q2);
  initMatrix(&sdep);
  initMatrix(&bias);
  initMatrix(&predicted_y);
  initMatrix(&predicted_residuals);


  MODELINPUT minpt;
  minpt.mx = &mx;
  minpt.my = &my;
  minpt.nlv = 10;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  //LeaveOneOut(&minpt, _PLS_, &predicted_y, &predicted_residuals, 4, NULL, 0);
  BootstrapRandomGroupsCV(&minpt, 5, 20, _PLS_, &predicted_y, &predicted_residuals, 4, NULL, 0);
  PLSRegressionStatistics(my, predicted_y, &q2, &sdep, &bias);

  puts("Q2 Cross Validation");
  PrintMatrix(q2);
  puts("SDEP Cross Validation");
  PrintMatrix(sdep);
  puts("BIAS Cross Validation");
  PrintMatrix(bias);

  DelMatrix(&predicted_y);
  DelMatrix(&predicted_residuals);
  DelMatrix(&q2);
  DelMatrix(&sdep);
  DelMatrix(&bias);
  DelMatrix(&mx);
  DelMatrix(&my);
}

void TestPLS12()
{
  size_t i, j, k, nvar;
  matrix *x, *x_, *y, *map; /* Data matrix */
  uivector *varselected, *vardistribution;
  PLSMODEL *m;

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

  printf("Test PLS 12\n");

  initUIVector(&varselected);
  initUIVector(&vardistribution);
  initMatrix(&map);
  /*LOO and Y Scrambling
  PLSSpearmannVariableSelection(x, y, NULL, NULL,
                                   1, 0, 999, 1, 0, 0,
                                   0.1,
                                   &varselected, &map, &vardistribution, 4, NULL);
*/
  puts("Variable Selected");
  PrintUIVector(varselected);

  puts("Models Map");
  PrintMatrix(map);

  puts("Variable Distribution");
  PrintUIVector(vardistribution);

  nvar = 0;
  for(i = 0; i < varselected->size; i++){
    if(getUIVectorValue(varselected, i) == 1){
      nvar++;
    }
    else{
      continue;
    }
  }

  NewMatrix(&x_, x->row, nvar);

  for(i = 0; i < x->row; i++){
    k = 0;
    for(j = 0; j < x->col; j++){
      if(getUIVectorValue(varselected, j) == 1){
        setMatrixValue(x_, i, k, getMatrixValue(x, i, j));
        k++;
      }
      else{
        continue;
      }
    }
  }

  PrintMatrix(x_);
  PrintMatrix(y);
  NewPLSModel(&m);
  PLS(x_, y, 100, 1, 0, m, NULL);

  MODELINPUT minpt;
  minpt.mx = &x_;
  minpt.my = &y;
  minpt.nlv = 10;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  LeaveOneOut(&minpt, _PLS_, &m->predicted_y, &m->pred_residuals, 4, NULL, 0);
  PLSRegressionStatistics(y, m->predicted_y, &m->q2y, &m->sdep, &m->bias);

  /*PrintPLSModel(m);*/
  puts("Q^2");
  PrintMatrix(m->q2y);

  puts("SDEP");
  PrintMatrix(m->sdep);

  puts("BIAS");
  PrintMatrix(m->bias);

  puts("LOO PREDICTED Y");
  PrintMatrix(m->predicted_y);

  puts("RealY");
  PrintMatrix(y);

  DelPLSModel(&m);

  DelMatrix(&map);
  DelUIVector(&varselected);
  DelUIVector(&vardistribution);
  DelMatrix(&x_);
  DelMatrix(&x);
  DelMatrix(&y);

}

void TestPLS11()
{
  matrix *x, *y;
  PLSMODEL *m;

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

  printf("Test PLS 11\n");

  NewPLSModel(&m);
  PLS(x, y, 100, 1, 0, m, NULL);

  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = 100;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  LeaveOneOut(&minpt, _PLS_, &m->predicted_y, &m->pred_residuals, 4, NULL, 0);
  PLSRegressionStatistics(y, m->predicted_y, &m->q2y, &m->sdep, &m->bias);

  ValidationArg varg;
  varg.vtype = LOO;
  YScrambling(&minpt, _PLS_, varg, 100, &m->r2q2scrambling, 4, NULL);

  /*PrintPLSModel(m);*/
  puts("Q^2");
  PrintMatrix(m->q2y);

  puts("SDEP");
  PrintMatrix(m->sdep);

  puts("BIAS");
  PrintMatrix(m->bias);

  puts("PREDICTED Y");
  PrintMatrix(m->predicted_y);

  puts("PREDICTED RESIDUALS");
  PrintMatrix(m->pred_residuals);

  puts("YSCRAMBLING RESULTS");
  PrintMatrix(m->r2q2scrambling);

  puts("RealY");
  PrintMatrix(y);

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);
}

void TestPLS10()
{
  size_t i, j, k, nvar;
  matrix *x, *x_, *y, *map; /* Data matrix */
  uivector *varselected, *vardistribution;
  PLSMODEL *m;

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

  printf("Test PLS 10\n");

  initUIVector(&varselected);
  initUIVector(&vardistribution);
  initMatrix(&map);
  /*LOO
  PSLPSOVariableSelection(x, y, NULL, NULL,
                       1, 0, 999, 1, 0, 0,
                       400, 1.0,
                       &varselected, &map, &vardistribution,  4, NULL);*/





  /*RG CV
  PSLPSOVariableSelection(x, y, NULL, NULL,
                       1, 0, 999, 2, 2, 20, 0,
                       30, 0.4,
                       &varselected, &map, &vardistribution, 4, NULL);*/




  /* GENETIC ALGORITM VARIABLE SELECTION  LOO
  PLSGAVariableSelection(x, y, NULL, NULL,
                       1, 0, 100, 1, 0, 0, 0,
                       40, 0.3, 0.5, 0, 0, 0.9,
                       &varselected, &map, &vardistribution, 4, NULL);*/

  /* GENETIC ALGORITM VARIABLE SELECTION  RG: 2 groups, 20 iterations..
  PLSGAVariableSelection(x, y, NULL, NULL,
                       1, 0, 1, 1, 7, 20,
                       50, 0.5, 0.5, 0.0, 0.0, 0.9,
                       &varselected, &map, &vardistribution, 4, NULL);
*/



  puts("Variable Selected");
  PrintUIVector(varselected);

  puts("Average of Variables Distribution");
  PrintUIVector(vardistribution);

  puts("Models Map");
  PrintMatrix(map);

  nvar = 0;
  for(i = 0; i < varselected->size; i++){
    if(getUIVectorValue(varselected, i) == 1){
      nvar++;
    }
    else{
      continue;
    }
  }

  NewMatrix(&x_, x->row, nvar);

  for(i = 0; i < x->row; i++){
    k = 0;
    for(j = 0; j < x->col; j++){
      if(getUIVectorValue(varselected, j) == 1){
        setMatrixValue(x_, i, k, getMatrixValue(x, i, j));
        k++;
      }
      else{
        continue;
      }
    }
  }

  NewPLSModel(&m);
  PLS(x_, y, 100, 1, 0, m, NULL);

  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = 100;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  LeaveOneOut(&minpt, _PLS_, &m->predicted_y, &m->pred_residuals, 4, NULL, 0);
  PLSRegressionStatistics(y, m->predicted_y, &m->q2y, &m->sdep, &m->bias);

  /*PrintPLSModel(m);*/
  puts("Q^2");
  PrintMatrix(m->q2y);

  puts("SDEP");
  PrintMatrix(m->sdep);

  puts("BIAS");
  PrintMatrix(m->bias);

  puts("LOO PREDICTED Y");
  PrintMatrix(m->predicted_y);

  puts("RealY");
  PrintMatrix(y);

  DelPLSModel(&m);

  DelUIVector(&vardistribution);
  DelMatrix(&map);
  DelUIVector(&varselected);
  DelMatrix(&x_);
  DelMatrix(&x);
  DelMatrix(&y);

}

/*
 * points are based on this relaction:
 *
 * y = x1 + x2/3 + 0.6666
 *
 * This test it's working for prediction
 *
 * The estimated value are equal to the exact value:

   y1 = 7,33
   y2 = 13,99
   y3 = 26,66
 *
 */
void TestPLS9()
{
  matrix *x, *y;

  PLSMODEL *m;

  matrix *r2y;
  matrix *sdec;
  matrix *bias;

  matrix *xpred;
  matrix *xpredscores;

  matrix *ypred;

  NewMatrix(&x, 4, 2);
  NewMatrix(&y, 4, 1);

  setMatrixValue(x, 0, 0, 5); setMatrixValue(x, 0, 1, 7);
  setMatrixValue(x, 1, 0, 3); setMatrixValue(x, 1, 1, 12);
  setMatrixValue(x, 2, 0, 2); setMatrixValue(x, 2, 1, 1);
  setMatrixValue(x, 3, 0, 4); setMatrixValue(x, 3, 1, 8);

  setMatrixValue(y, 0, 0, 8);
  setMatrixValue(y, 1, 0, 7.66);
  setMatrixValue(y, 2, 0, 3.0);
  setMatrixValue(y, 3, 0, 7.33);


  NewMatrix(&xpred, 3, 2);
  setMatrixValue(xpred, 0, 0, 5); setMatrixValue(xpred, 0, 1, 5); /* y = 7,33 */
  setMatrixValue(xpred, 1, 0, 12); setMatrixValue(xpred, 1, 1, 4); /* y = 13,99 */
  setMatrixValue(xpred, 2, 0, 21); setMatrixValue(xpred, 2, 1, 15); /* y = 26,66 */

  puts("Test 9 - PLS Model and prediction");
  /*Allocate the final output*/
  NewPLSModel(&m);

  puts("X:");
  PrintMatrix(x);

  puts("Y:");
  PrintMatrix(y);

  ssignal run = SIGSCIENTIFICRUN;
  PLS(x, y, 4, 1, 0, m, &run);

  initMatrix(&r2y);
  initMatrix(&sdec);
  initMatrix(&bias);

  PLSYPredictorAllLV(x, m, NULL, &m->recalculated_y);
  PLSRegressionStatistics(y, m->recalculated_y, &r2y, &sdec, &bias);

  PrintPLSModel(m);

  puts("R^2 for Y Scores");
  PrintMatrix(r2y);

  puts("SDEC");
  PrintMatrix(sdec);

  puts("BIAS");
  PrintMatrix(bias);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, &xpredscores);

  PLSYPredictor(xpredscores, m, 1, &ypred);

  puts("\nPrediction scores...");
  puts("x");
  PrintMatrix(xpredscores);

  puts("Extimated Y:");
  PrintMatrix(ypred);
  DelMatrix(&bias);
  DelMatrix(&sdec);
  DelMatrix(&r2y);

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);

  DelMatrix(&xpredscores);

  DelMatrix(&xpred);
  DelMatrix(&ypred);
}

void TestPLS8()
{
  matrix *x, *y;

  PLSMODEL *m;

  matrix *xpred;
  matrix *xpredscores;

  matrix *ypred;

  NewMatrix(&x, 3, 2);
  NewMatrix(&y, 3, 1);

  setMatrixValue(x, 0, 0, 4); setMatrixValue(x, 0, 1, 3);
  setMatrixValue(x, 1, 0, 2); setMatrixValue(x, 1, 1, 2);
  setMatrixValue(x, 2, 0, 5); setMatrixValue(x, 2, 1, 2);

  setMatrixValue(y, 0, 0, 50);
  setMatrixValue(y, 1, 0, 86);
  setMatrixValue(y, 2, 0, 20);


  NewMatrix(&xpred, 1, 2);
  setMatrixValue(xpred, 0, 0, 62); setMatrixValue(xpred, 0, 1, 1);

  puts("Test 8");
  /*Allocate the final output*/
  NewPLSModel(&m);



  puts("X:");
  PrintMatrix(x);

  puts("Y:");
  PrintMatrix(y);

  PLS(x, y, 4, 1, 0, m, NULL);

  PrintPLSModel(m);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, &xpredscores);

  PLSYPredictor(xpredscores, m, 2, &ypred);


  puts("\nPrediction scores...");
  puts("x");
  PrintMatrix(xpredscores);

  puts("Extimated Y:");
  PrintMatrix(ypred);

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);

  DelMatrix(&xpredscores);

  DelMatrix(&xpred);
  DelMatrix(&ypred);
}

void TestPLS7()
{
  matrix *x, *y;

  PLSMODEL *m;

  matrix *xpred;
  matrix *xpredscores;

  matrix *ypred;

  NewMatrix(&x, 3, 2);
  NewMatrix(&y, 3, 1);

  setMatrixValue(x, 0, 0, 37); setMatrixValue(x, 0, 1, 12);
  setMatrixValue(x, 1, 0, 62); setMatrixValue(x, 1, 1, 40);
  setMatrixValue(x, 2, 0, 13); setMatrixValue(x, 2, 1, 2);

  setMatrixValue(y, 0, 0, 50);
  setMatrixValue(y, 1, 0, 86);
  setMatrixValue(y, 2, 0, 20);


  NewMatrix(&xpred, 1, 2);
  setMatrixValue(xpred, 0, 0, 62); setMatrixValue(xpred, 0, 1, 1);

  puts("Test 7");
  /*Allocate the final output*/
  NewPLSModel(&m);

  puts("X:");
  PrintMatrix(x);

  puts("Y:");
  PrintMatrix(y);

  PLS(x, y, 4, 1, 0, m, NULL);

  PrintPLSModel(m);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, &xpredscores);

  PLSYPredictor(xpredscores, m, 2, &ypred);


  puts("\nPrediction scores...");
  puts("x");
  PrintMatrix(xpredscores);

  puts("Extimated Y");
  PrintMatrix(ypred);

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);

  DelMatrix(&xpredscores);

  DelMatrix(&xpred);
  DelMatrix(&ypred);
}

void TestPLS6()
{
  matrix *x, *y;

  PLSMODEL *m;
  matrix *xpred;
  matrix *xpredscores;

  matrix *ypred;

  NewMatrix(&x, 3, 4);
  NewMatrix(&y, 3, 1);

  setMatrixValue(x, 0, 0, 37); setMatrixValue(x, 0, 1, 12); setMatrixValue(x, 0, 2, 4); setMatrixValue(x, 0, 3, 3);
  setMatrixValue(x, 1, 0, 62); setMatrixValue(x, 1, 1, 40); setMatrixValue(x, 1, 2, 2); setMatrixValue(x, 1, 3, 2);
  setMatrixValue(x, 2, 0, 13); setMatrixValue(x, 2, 1, 2); setMatrixValue(x, 2, 2, 5); setMatrixValue(x, 2, 3, 2);

  setMatrixValue(y, 0, 0, 50);
  setMatrixValue(y, 1, 0, 86);
  setMatrixValue(y, 2, 0, 20);


  NewMatrix(&xpred, 1, 4);
  setMatrixValue(xpred, 0, 0, 62); setMatrixValue(xpred, 0, 1, 1); setMatrixValue(xpred, 0, 2, 62); setMatrixValue(xpred, 0, 3, 1);

  puts("Test 6");
  /*Allocate the final output*/
  NewPLSModel(&m);

  puts("X:");
  PrintMatrix(x);

  puts("Y:");
  PrintMatrix(y);

  PLS(x, y, 4, 1, 0, m, NULL);

  PrintPLSModel(m);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 4, &xpredscores);

  PLSYPredictor(xpredscores, m, 4, &ypred);

  puts("\nPrediction scores...");
  puts("x");
  PrintMatrix(xpredscores);

  puts("Extimated Y");
  PrintMatrix(ypred);

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);

  DelMatrix(&xpredscores);

  DelMatrix(&xpred);
  DelMatrix(&ypred);
}

void TestPLS5()
{

  matrix *x, *y, *predicted_y, *pred_residuals; /* Data matrix */
  matrix *q2y;
  matrix *sdep;

  NewMatrix(&x, 14, 6);
  NewMatrix(&y, 14, 1);

  x->data[0][0] = 4.000; x->data[0][1] = 4.000; x->data[0][2] = 1.000; x->data[0][3] = 84.140; x->data[0][4] = 1.050; x->data[0][5] = 235.150;
  x->data[1][0] = 5.000; x->data[1][1] = 5.000; x->data[1][2] = 1.000; x->data[1][3] = 79.100; x->data[1][4] = 0.978; x->data[1][5] = 231.000;
  x->data[2][0] = 4.000; x->data[2][1] = 5.000; x->data[2][2] = 1.000; x->data[2][3] = 67.090; x->data[2][4] = 0.970; x->data[2][5] = 249.000;
  x->data[3][0] = 4.000; x->data[3][1] = 4.000; x->data[3][2] = 1.000; x->data[3][3] = 68.070; x->data[3][4] = 0.936; x->data[3][5] = 187.350;
  x->data[4][0] = 3.000; x->data[4][1] = 4.000; x->data[4][2] = 2.000; x->data[4][3] = 68.080; x->data[4][4] = 1.030; x->data[4][5] = 363.000;
  x->data[5][0] = 9.000; x->data[5][1] = 7.000; x->data[5][2] = 1.000; x->data[5][3] = 129.160; x->data[5][4] = 1.090; x->data[5][5] = 258.000;
  x->data[6][0] = 10.000; x->data[6][1] = 8.000; x->data[6][2] = 0.000; x->data[6][3] = 128.160; x->data[6][4] = 1.150; x->data[6][5] = 352.00;
  x->data[7][0] = 6.000; x->data[7][1] = 6.000; x->data[7][2] = 0.000; x->data[7][3] = 78.110; x->data[7][4] = 0.880; x->data[7][5] = 278.64;
  x->data[8][0] = 16.000; x->data[8][1] = 10.000; x->data[8][2] = 0.000; x->data[8][3] = 202.260; x->data[8][4] = 1.271; x->data[8][5] = 429.15;
  x->data[9][0] = 6.000; x->data[9][1] = 12.000; x->data[9][2] = 0.000; x->data[9][3] = 84.160; x->data[9][4] = 0.780; x->data[9][5] = 279.00;
  x->data[10][0] = 4.000; x->data[10][1] = 8.000; x->data[10][2] = 1.000; x->data[10][3] = 72.110; x->data[10][4] = 0.890; x->data[10][5] = 164.50;
  x->data[11][0] = 4.000; x->data[11][1] = 9.000; x->data[11][2] = 1.000; x->data[11][3] = 71.110; x->data[11][4] = 0.866; x->data[11][5] = 210.00;
  x->data[12][0] = 5.000; x->data[12][1] = 11.000; x->data[12][2] = 1.000; x->data[12][3] = 85.150; x->data[12][4] = 0.862; x->data[12][5] = 266.00;
  x->data[13][0] = 5.000; x->data[13][1] = 10.000; x->data[13][2] = 1.000; x->data[13][3] = 86.130; x->data[13][4] = 0.880; x->data[13][5] = 228.00;

  y->data[0][0] = 357.150;
  y->data[1][0] = 388.000;
  y->data[2][0] = 403.000;
  y->data[3][0] = 304.550;
  y->data[4][0] = 529.000;
  y->data[5][0] = 510.000;
  y->data[6][0] = 491.000;
  y->data[7][0] = 353.300;
  y->data[8][0] = 666.650;
  y->data[9][0] = 354.000;
  y->data[10][0] = 339.000;
  y->data[11][0] = 360.000;
  y->data[12][0] = 379.000;
  y->data[13][0] = 361.000;


  printf("Test PLS 5\n");

  puts("Matrix X and Y");
  PrintMatrix(x);
  PrintMatrix(y);


  /*Allocate the final output*/
  initMatrix(&q2y);
  initMatrix(&sdep);

  ssignal run = SIGSCIENTIFICRUN;

  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = 999;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  initMatrix(&predicted_y);
  initMatrix(&pred_residuals);
  BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, &predicted_y, &pred_residuals, 4, &run, 0);
  PLSRegressionStatistics(y, predicted_y, &q2y, &sdep, NULL);

  puts("Q2 Cross Validation");
  PrintMatrix(q2y);

  puts("SDEP Cross Validation");
  PrintMatrix(sdep);

  DelMatrix(&sdep);
  DelMatrix(&q2y);
  DelMatrix(&x);
  DelMatrix(&y);
  DelMatrix(&predicted_y);
  DelMatrix(&pred_residuals);
}



void TestPLS4()
{

  matrix *x, *y; /* Data matrix */

  PLSMODEL *m;

  matrix *xpred/*, *ypred*/;
  matrix *xpredscores;

  matrix *ypred;

  NewMatrix(&x, 12, 6);
  NewMatrix(&y, 12, 1);

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



  NewMatrix(&xpred, 2, 6);
  /*NewMatrix(&ypred, 2, 1);*/

  xpred->data[0][0] = 5.0000;  xpred->data[0][1] = 10.0000;  xpred->data[0][2] = 1.0000;  xpred->data[0][3] = 86.1300;  xpred->data[0][4] = 0.8800;  xpred->data[0][5] = 228.0000;
  xpred->data[1][0] = 5.0000;  xpred->data[1][1] = 11.0000;  xpred->data[1][2] = 1.0000;  xpred->data[1][3] = 85.1500;  xpred->data[1][4] = 0.8620;  xpred->data[1][5] = 266.0000;

  /*ypred->data[0][0] = 361.0000;
  ypred->data[1][0] = 379.0000;*/

  printf("Test PLS 3\n");


    /*Allocate the final output*/
  NewPLSModel(&m);

  PLS(x, y, 7, 1, 0, m, NULL);

  PrintPLSModel(m);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, &xpredscores);

  PLSYPredictor(xpredscores, m, 2, &ypred);

  puts("\nPrediction scores...");
  puts("x");
  PrintMatrix(xpredscores);

  puts("Extimated Y");
  PrintMatrix(ypred);

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);

  DelMatrix(&xpredscores);
  DelMatrix(&xpred);
  DelMatrix(&ypred);
}

void TestPLS3()
{

  matrix *x, *y; /* Data matrix */

  PLSMODEL *m;

  matrix *xpred/*, *ypred*/;
  matrix *xpredscores;
  matrix *ypred;

  NewMatrix(&x, 13, 6);
  NewMatrix(&y, 13, 1);

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


  NewMatrix(&xpred, 1, 6);
  /*NewMatrix(&ypred, 1, 1);*/

  xpred->data[0][0] = 5.0000;  xpred->data[0][1] = 10.0000;  xpred->data[0][2] = 1.0000;  xpred->data[0][3] = 86.1300;  xpred->data[0][4] = 0.8800;  xpred->data[0][5] = 228.0000;

  /*ypred->data[0][0] = 361.0000;*/

  printf("Test PLS 3\n");

  /*Allocate the final output*/
  NewPLSModel(&m);

  PLS(x, y, 7, 1, 0, m, NULL);

  PrintPLSModel(m);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, &xpredscores);

  PLSYPredictor(xpredscores, m, 2, &ypred);

  puts("\nPrediction scores...");
  puts("x");
  PrintMatrix(xpredscores);

  puts("Extimated Y");
  PrintMatrix(ypred);

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);

  DelMatrix(&xpredscores);

  DelMatrix(&xpred);
  DelMatrix(&ypred);
}

void TestPLS2_2()
{
  matrix *x, *y;

  PLSMODEL *m;

  NewMatrix(&x, 5, 4);
  NewMatrix(&y, 5, 3);

  setMatrixValue(x, 0, 0, 7); setMatrixValue(x, 0, 1, 7); setMatrixValue(x, 0, 2, 13); setMatrixValue(x, 0, 3, 7);
  setMatrixValue(x, 1, 0, 4); setMatrixValue(x, 1, 1, 3); setMatrixValue(x, 1, 2, 14); setMatrixValue(x, 1, 3, 7);
  setMatrixValue(x, 2, 0, 10); setMatrixValue(x, 2, 1, 5); setMatrixValue(x, 2, 2, 12); setMatrixValue(x, 2, 3, 5);
  setMatrixValue(x, 3, 0, 16); setMatrixValue(x, 3, 1, 7); setMatrixValue(x, 3, 2, 11); setMatrixValue(x, 3, 3, 3);
  setMatrixValue(x, 4, 0, 13); setMatrixValue(x, 4, 1, 3); setMatrixValue(x, 4, 2, 10); setMatrixValue(x, 4, 3, 3);

  setMatrixValue(y, 0, 0, 14); setMatrixValue(y, 0, 1, 7); setMatrixValue(y, 0, 2, 8);
  setMatrixValue(y, 1, 0, 10); setMatrixValue(y, 1, 1, 7); setMatrixValue(y, 1, 2, 6);
  setMatrixValue(y, 2, 0, 8); setMatrixValue(y, 2, 1, 5); setMatrixValue(y, 2, 2, 5);
  setMatrixValue(y, 3, 0, 2); setMatrixValue(y, 3, 1, 4); setMatrixValue(y, 3, 2, 7);
  setMatrixValue(y, 4, 0, 6); setMatrixValue(y, 4, 1, 3); setMatrixValue(y, 4, 2, 4);



  printf("Test PLS 2\n");
  NewPLSModel(&m);
  ssignal run = SIGSCIENTIFICSTOP;
  PLS(x, y, 4, 1, 1, m, &run);

  PrintPLSModel(m);

  DelPLSModel(&m);
  DelMatrix(&y);
  DelMatrix(&x);
}

void TestPLS2()
{
  size_t row, col;
  matrix *x, *y;

  PLSMODEL *m;

  NewMatrix(&x, 5, 4);
  NewMatrix(&y, 5, 3);

  setMatrixValue(x, 0, 0, 7); setMatrixValue(x, 0, 1, 7); setMatrixValue(x, 0, 2, 13); setMatrixValue(x, 0, 3, 7);
  setMatrixValue(x, 1, 0, 4); setMatrixValue(x, 1, 1, 3); setMatrixValue(x, 1, 2, 14); setMatrixValue(x, 1, 3, 7);
  setMatrixValue(x, 2, 0, 10); setMatrixValue(x, 2, 1, 5); setMatrixValue(x, 2, 2, 12); setMatrixValue(x, 2, 3, 5);
  setMatrixValue(x, 3, 0, 16); setMatrixValue(x, 3, 1, 7); setMatrixValue(x, 3, 2, 11); setMatrixValue(x, 3, 3, 3);
  setMatrixValue(x, 4, 0, 13); setMatrixValue(x, 4, 1, 3); setMatrixValue(x, 4, 2, 10); setMatrixValue(x, 4, 3, 3);

  setMatrixValue(y, 0, 0, 14); setMatrixValue(y, 0, 1, 7); setMatrixValue(y, 0, 2, 8);
  setMatrixValue(y, 1, 0, 10); setMatrixValue(y, 1, 1, 7); setMatrixValue(y, 1, 2, 6);
  setMatrixValue(y, 2, 0, 8); setMatrixValue(y, 2, 1, 5); setMatrixValue(y, 2, 2, 5);
  setMatrixValue(y, 3, 0, 2); setMatrixValue(y, 3, 1, 4); setMatrixValue(y, 3, 2, 7);
  setMatrixValue(y, 4, 0, 6); setMatrixValue(y, 4, 1, 3); setMatrixValue(y, 4, 2, 4);

  printf("Test PLS 2\n");
  NewPLSModel(&m);
  ssignal run = SIGSCIENTIFICRUN;

  PLS(x, y, 3, 0, 0, m, &run);

  PrintPLSModel(m);

  /*VALIDATE THE MODEL */
  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = 999;
  minpt.xautoscaling = 0;
  minpt.yautoscaling = 0;

  //BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, &m->predicted_y, &m->pred_residuals, 4, NULL, 0);
  LeaveOneOut(&minpt, _PLS_, &m->predicted_y, &m->pred_residuals, 4, &run, 0);
  PLSRegressionStatistics(y, m->predicted_y, &m->q2y, &m->sdep, &m->bias);
  PrintMatrix(m->predicted_y);
  puts("Q2 Cross Validation");
  PrintMatrix(m->q2y);
  puts("SDEP Cross Validation");
  PrintMatrix(m->sdep);
  puts("BIAS Cross Validation");
  PrintMatrix(m->bias);

  puts("PREDICTED Y");
  PrintMatrix(m->predicted_y);

  puts("REAL Y");
  PrintMatrix(y);


  MatrixGetMaxValue(m->sdep, &row, &col);

  printf("cutoff : %lu  %lu\n", row, col);

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);
}

void TestPLS1()
{
  matrix *x, *y; /* Data matrix */
  dvector *betas;
  PLSMODEL *m;

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

  printf("Test PLS 1\n");

  PrintMatrix(x);
  PrintMatrix(y);

  /*Allocate the final output*/
  NewPLSModel(&m);

  PLS(x, y, 777, 0, 0, m, NULL);

  PrintPLSModel(m);

  PrintMatrix(y);


  /*VALIDATE THE MODEL */
  MODELINPUT minpt;
  minpt.mx = &x;
  minpt.my = &y;
  minpt.nlv = 999;
  minpt.xautoscaling = 0;
  minpt.yautoscaling = 0;

  BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, &m->predicted_y, &m->pred_residuals, 4, NULL, 0);
  //LeaveOneOut(&minpt, _PLS_, &m->predicted_y, &m->pred_residuals, 4, NULL, 0);
  PLSRegressionStatistics(y, m->predicted_y, &m->q2y, &m->sdep, &m->bias);
  PrintMatrix(m->predicted_y);
  puts("Q2 Cross Validation");
  PrintMatrix(m->q2y);
  puts("SDEP Cross Validation");
  PrintMatrix(m->sdep);
  puts("BIAS Cross Validation");
  PrintMatrix(m->bias);

  puts("Q^2");
  PrintMatrix(m->q2y);

  puts("SDEP");
  PrintMatrix(m->sdep);

  puts("BIAS");
  PrintMatrix(m->bias);

  puts("BETA COEFFICIENTS");
  initDVector(&betas);
  PLSBetasCoeff(m, GetLVCCutoff(m->q2y), &betas);
  PrintDVector(betas);

  puts("PREDICTED Y");
  PrintMatrix(m->predicted_y);

  puts("PREDICTED RESIDUALS");
  PrintMatrix(m->pred_residuals);

  puts("REAL Y");
  PrintMatrix(y);

  DelDVector(&betas);
  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);
}

int main(void)
{
  /*test 1- 5*/
  TestPLS1();
  /*TestPLS2();
  TestPLS3();
  TestPLS4();
  TestPLS5();*/

  /*test 6-9
  TestPLS6();
  TestPLS7();
  TestPLS8();
  TestPLS9();*/
  //TestPLS10();
TestPLS11();
  //TestPLS12();
  /*TestPLS13();
  TestPLS14();*/
//TestPLS15();
  return 0;
}

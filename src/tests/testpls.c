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

#include "tensor.h"
#include "pls.h"
#include "modelvalidation.h"
#include "matrix.h"
#include "vector.h"
#include "numeric.h"
#include "datasets.h"
#include <time.h>

void TestPLS15()
{
  printf("Test PLS 15 - check PLS convergence criteria: ");
  size_t i;
  size_t j;
  size_t id;
  matrix *x1;
  matrix *y1;
  matrix *x2;
  matrix *y2;
  uivector *ids;
  PLSMODEL *mod1;
  PLSMODEL *mod2;
  initMatrix(&x1);
  initMatrix(&y1);
  boston_house_price(x1, y1);
  NewMatrix(&x2, x1->row, x1->col);
  NewMatrix(&y2, y1->row, y1->col);
  srand_(time(NULL));
  initUIVector(&ids);
  i = 0;
  while( i < x1->row){
    id = randInt(0, x1->row);
    if(UIVectorHasValue(ids, id) == 0){
      continue;
    }
    else{
      UIVectorAppend(ids, id);
      for(j = 0; j < x1->col; j++){
        x2->data[i][j] = x1->data[id][j];
      }
      for(j = 0; j < y1->col; j++){
        y2->data[i][j] = y1->data[id][j];
      }
      i++;
    }
  }
  /*Allocate the final output*/
  NewPLSModel(&mod1);
  NewPLSModel(&mod2);
  PLS(x1, y1, 3, 1, 0, mod1, NULL);
  PLS(x2, y2, 3, 1, 0, mod2, NULL);

  for(i = 0; i < mod1->xscores->row; i++){
    id = ids->data[i];
    if(FLOAT_EQ(mod1->xscores->data[id][0], mod2->xscores->data[i][0], 1e-5) &&
       FLOAT_EQ(mod1->xscores->data[id][1], mod2->xscores->data[i][1], 1e-5)){
      continue;
    }
    else{
      printf("%f %f\n", mod1->xscores->data[id][0], mod2->xscores->data[i][0]);
      printf("%f %f\n", mod1->xscores->data[id][1], mod2->xscores->data[i][1]);
      abort();
    }
  }
  printf("OK.\n");
  
  DelUIVector(&ids);
  DelPLSModel(&mod1);
  DelPLSModel(&mod2);
  DelMatrix(&x1);
  DelMatrix(&y1);
  DelMatrix(&x2);
  DelMatrix(&y2);
}

void TestPLS14()
{
  printf("Test PLS 14 - check PLS Score prediction: ");
  size_t i, j;
  matrix *x, *y, *ps; /* Data matrix */
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

  /*Allocate the final output*/
  NewPLSModel(&m);

  PLS(x, y, 3, 0, 0, m, NULL);
  initMatrix(&ps);

  PLSScorePredictor(x, m, 3, ps);
  for(i = 0; i < m->xscores->row; i++){
    for(j = 0; j < m->xscores->col; j++){
      if(FLOAT_EQ(m->xscores->data[i][j], ps->data[i][j], 1E-6)){
        continue;
      }
      else{
        abort();
      }
    }
  }
  printf("OK.\n");
  DelPLSModel(&m);
  DelMatrix(&ps);
  DelMatrix(&x);
  DelMatrix(&y);
}

void TestPLS13()
{
  printf("Test 13 - PLS Regression Model and validation via KFoldCV: ");
  matrix *x, *y;
  uivector *groups;

  PLSMODEL *m;

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
  NewUIVector(&groups, nobj);

  srand_(nobj+nfeats+ntarg);
  for(size_t i = 0; i < nobj; i++){
    for(size_t j = 0; j < nfeats; j++){
      setMatrixValue(x, i, j, randDouble(-100,100));
    }

    double y_ = 0.f;
    for(size_t j = 0; j < nfeats; j++){
      y_ += x->data[i][j];
    }
    y_ /= (double)x->row;

    for(size_t j = 0; j < ntarg; j++){
      setMatrixValue(y, i, j, y_);
      //setMatrixValue(y, i, j, randDouble(0, 2));
    }
    groups->data[i] = randInt(0, 6);
    //groups->data[i] = randInt(0, 6);
  }

  for(size_t i = 0; i < nobj_to_pred; i++){
    for(size_t j = 0; j < nfeats; j++){
      setMatrixValue(xpred, i, j, randDouble(-100,100));
    }
  }

  /*Allocate the final output*/
  NewPLSModel(&m);

  /*puts("X:");
  PrintMatrix(x);

  puts("Y:");
  PrintMatrix(y);*/

  ssignal run = SIGSCIENTIFICRUN;
  PLS(x, y, 4, 1, 0, m, &run);


  /*VALIDATE THE MODEL */
  MODELINPUT minpt = initModelInput();
  minpt.mx = x;
  minpt.my = y;
  minpt.nlv = 4;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  KFoldCV(&minpt, groups, _PLS_, m->predicted_y, m->pred_residuals, 4, NULL, 0);

  //PrintMatrix(m->predicted_y);

  PLSRegressionStatistics(y, m->predicted_y, m->q2y, m->sdep, m->bias);

  printf("OK.\n");
  /*PrintMatrix(y);
  PrintMatrix(m->predicted_y);*/
  puts("Q2 Cross Validation");
  PrintMatrix(m->q2y);
  puts("SDEP Cross Validation");
  PrintMatrix(m->sdep);
  puts("BIAS Cross Validation");
  PrintMatrix(m->bias);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, xpredscores);

  PLSYPredictor(xpredscores, m, 1, ypred);

  /*puts("\nPrediction scores...");
  puts("x");
  PrintMatrix(xpredscores);

  puts("Extimated Y:");
  PrintMatrix(ypred);
  */

  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);
  DelUIVector(&groups);
  DelMatrix(&xpredscores);

  DelMatrix(&xpred);
  DelMatrix(&ypred);

}

void TestPLS12()
{
  printf("Test PLS 12 - PLS DA Model and prediction using random data\n");
  matrix *x, *y;

  PLSMODEL *m;

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

  srand_(nobj+nfeats+ntarg);
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

  /*Allocate the final output*/
  NewPLSModel(&m);

  /*puts("X:");
  PrintMatrix(x);

  puts("Y:");
  PrintMatrix(y);*/

  ssignal run = SIGSCIENTIFICRUN;
  PLS(x, y, 4, 1, 0, m, &run);

  /*VALIDATE THE MODEL */
  MODELINPUT minpt = initModelInput();
  minpt.mx = x;
  minpt.my = y;
  minpt.nlv = 4;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, m->predicted_y, m->pred_residuals, 4, NULL, 0);
  PrintMatrix(m->predicted_y);

  PLSDiscriminantAnalysisStatistics(y,
                                    m->predicted_y,
                                    m->roc_validation,
                                    m->roc_auc_validation,
                                    m->precision_recall_validation,
                                    m->precision_recall_ap_validation);
  
  PrintPLSModel(m);

  puts("ROC AUC's for Trainig set");
  PrintMatrix(m->roc_auc_validation);
  puts("ROC Curves for each LV");
  PrintTensor(m->roc_validation);

  puts("Precision-Recall or Trainig set");
  PrintMatrix(m->precision_recall_ap_validation);
  puts("Precision-Recall Curves for each LV");
  PrintTensor(m->precision_recall_validation);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, xpredscores);

  PLSYPredictor(xpredscores, m, 1, ypred);



  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);

  DelMatrix(&xpredscores);

  DelMatrix(&xpred);
  DelMatrix(&ypred);
}

void TestPLS11()
{
  size_t nobj = 800;
  matrix *mx, *my;
  matrix *q2, *sdep, *bias;
  matrix *predicted_y;
  matrix *predicted_residuals;
  puts("Test PLS 11: Simple Calculation PLS Model with RGCV using random data");

  NewMatrix(&mx, nobj, 128);
  NewMatrix(&my, nobj, 1);

  srand_(nobj);
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


  MODELINPUT minpt = initModelInput();
  minpt.mx = mx;
  minpt.my = my;
  minpt.nlv = 10;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  //LeaveOneOut(&minpt, _PLS_, predicted_y, predicted_residuals, 4, NULL, 0);
  BootstrapRandomGroupsCV(&minpt, 5, 20, _PLS_, predicted_y, predicted_residuals, 4, NULL, 0);
  PLSRegressionStatistics(my, predicted_y, q2, sdep, bias);

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

void TestPLS10()
{
  printf("Test PLS 10: YScrambling test\n");
  matrix *x, *y;
  PLSMODEL *m;

  initMatrix(&x);
  initMatrix(&y);
  residential_building(x, y);

  NewPLSModel(&m);
  PLS(x, y, 5, 1, 1, m, NULL);

  MODELINPUT minpt = initModelInput();
  minpt.mx = x;
  minpt.my = y;
  minpt.nlv = 5;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 1;

  LeaveOneOut(&minpt, _PLS_, m->predicted_y, m->pred_residuals, 1, NULL, 0);
  PLSRegressionStatistics(y, m->predicted_y, m->q2y, m->sdep, m->bias);

  ValidationArg varg = initValidationArg();
  varg.vtype = BootstrapRGCV;
  YScrambling(&minpt, _PLS_, varg, 2, m->yscrambling, 1, NULL);

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
  PrintMatrix(m->yscrambling);

  puts("RealY");
  PrintMatrix(y);

  DelPLSModel(&m);
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
  puts("Test PLS 9: another simple prediction test");
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

  PLSYPredictorAllLV(x, m, NULL, m->recalculated_y);
  PLSRegressionStatistics(y, m->recalculated_y, r2y, sdec, bias);

  PrintPLSModel(m);

  puts("R^2 for Y Scores");
  PrintMatrix(r2y);

  puts("SDEC");
  PrintMatrix(sdec);

  puts("BIAS");
  PrintMatrix(bias);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, xpredscores);

  PLSYPredictor(xpredscores, m, 1, ypred);

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
  puts("Test PLS 8: another simple prediction test");
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

  PLSScorePredictor(xpred, m, 2, xpredscores);

  PLSYPredictor(xpredscores, m, 2, ypred);


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
  puts("Test PLS 7: another simple prediction test");
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

  PLSScorePredictor(xpred, m, 2, xpredscores);

  PLSYPredictor(xpredscores, m, 2, ypred);


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
  puts("Test PLS 6: a simple prediction test with 3 task");
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

  PLSScorePredictor(xpred, m, 4, xpredscores);

  PLSYPredictor(xpredscores, m, 4, ypred);

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
  printf("Test PLS 5: test multitask partial least squares  \n");
  matrix *x, *y, *predicted_y, *pred_residuals; /* Data matrix */
  matrix *q2y;
  matrix *sdep;

  initMatrix(&x);
  initMatrix(&y);
  residential_building(x, y);

  puts("Matrix X and Y");
  PrintMatrix(x);
  PrintMatrix(y);

  /*Allocate the final output*/
  initMatrix(&q2y);
  initMatrix(&sdep);

  ssignal run = SIGSCIENTIFICRUN;

  MODELINPUT minpt = initModelInput();
  minpt.mx = x;
  minpt.my = y;
  minpt.nlv = 5;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 1;

  initMatrix(&predicted_y);
  initMatrix(&pred_residuals);
  BootstrapRandomGroupsCV(&minpt, 5, 20, _PLS_, predicted_y, pred_residuals, 4, &run, 0);
  PLSRegressionStatistics(y, predicted_y, q2y, sdep, NULL);

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
  printf("Test PLS 4: Prediction of two instance \n");
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


    /*Allocate the final output*/
  NewPLSModel(&m);

  PLS(x, y, 7, 1, 0, m, NULL);

  PrintPLSModel(m);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, xpredscores);

  PLSYPredictor(xpredscores, m, 2, ypred);

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
  printf("Test PLS 3: Prediction of one instance \n");
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

  xpred->data[0][0] = 5.0000;  xpred->data[0][1] = 10.0000;  xpred->data[0][2] = 1.0000;  xpred->data[0][3] = 86.1300;  xpred->data[0][4] = 0.8800;  xpred->data[0][5] = 228.0000;

  /*Allocate the final output*/
  NewPLSModel(&m);

  PLS(x, y, 7, 1, 0, m, NULL);

  PrintPLSModel(m);

  initMatrix(&xpredscores);
  initMatrix(&ypred);

  PLSScorePredictor(xpred, m, 2, xpredscores);

  PLSYPredictor(xpredscores, m, 2, ypred);

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

void TestPLS2()
{
  printf("Test PLS 2: boston house pricing\n");
  matrix *x, *y;

  PLSMODEL *m;

  initMatrix(&x);
  initMatrix(&y);
  boston_house_price(x,  y);

  NewPLSModel(&m);
  ssignal run = SIGSCIENTIFICRUN;

  PLS(x, y, 5, 0, 0, m, &run);

  /*VALIDATE THE MODEL */
  MODELINPUT minpt = initModelInput();
  minpt.mx = x;
  minpt.my = y;
  minpt.nlv = 5;
  minpt.xautoscaling = 1;
  minpt.yautoscaling = 0;

  //BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, m->predicted_y, m->pred_residuals, 4, NULL, 0);
  LeaveOneOut(&minpt, _PLS_, m->predicted_y, m->pred_residuals, 4, &run, 0);
  PLSRegressionStatistics(y, m->predicted_y, m->q2y, m->sdep, m->bias);
  PrintMatrix(m->predicted_y);
  puts("Q2 Cross Validation");
  PrintMatrix(m->q2y);
  puts("SDEP Cross Validation");
  PrintMatrix(m->sdep);
  puts("BIAS Cross Validation");
  PrintMatrix(m->bias);

  puts("PREDICTED Y");
  PrintMatrix(m->predicted_y);

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

  PLS(x, y, 3, 1, 0, m, NULL);

  PrintPLSModel(m);
  
  /*VALIDATE THE MODEL */
  MODELINPUT minpt = initModelInput();
  minpt.mx = x;
  minpt.my = y;
  minpt.nlv = 3;
  minpt.xautoscaling = 0;
  minpt.yautoscaling = 0;

  BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, m->predicted_y, m->pred_residuals, 1, NULL, 0);
  //LeaveOneOut(&minpt, _PLS_, m->predicted_y, m->pred_residuals, 4, NULL, 0);
  PLSRegressionStatistics(y, m->predicted_y, m->q2y, m->sdep, m->bias);
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

  PLSBetasCoeff(m, GetLVCCutoff(m->q2y)+1, betas);
  PrintDVector(betas);

  puts("PREDICTED Y");
  PrintMatrix(m->predicted_y);

  puts("PREDICTED RESIDUALS");
  PrintMatrix(m->pred_residuals);

  DelDVector(&betas);
  DelPLSModel(&m);
  DelMatrix(&x);
  DelMatrix(&y);
}

int main(void)
{
  /*test 1- 5*/
  TestPLS1();
  TestPLS2();
  TestPLS3();
  TestPLS4();
  TestPLS5();

  /*test 6-13*/
  TestPLS6();
  TestPLS7();
  TestPLS8();
  TestPLS9();
  //TestPLS10();
  TestPLS11();
  TestPLS12();
  TestPLS13();
  TestPLS14();
  TestPLS15();
  return 0;
}

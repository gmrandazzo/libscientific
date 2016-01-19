/* testmwpls.c
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "array.h"
#include "upls.h"
#include "upca.h"

void test9()
{
  array *ax;
  array *ay;

  dvector *q2x;
  array *q2y;
  array *sdep;
  array *yscrabl_q2, *yscrabl_sdep;

  UPLSMODEL *m;

  puts(">>>>>>> Test 9: Compute Multi Way LOO Cross Validation");

  NewArray(&ax, 2);

  NewArrayMatrix(&ax, 0, 10, 3);
  NewArrayMatrix(&ax, 1, 10, 3);

  setArrayValue(ax, 0, 0, 0, 37);  setArrayValue(ax, 0, 0, 1, 12);  setArrayValue(ax, 0, 0, 2, -1);
  setArrayValue(ax, 0, 1, 0, 62);  setArrayValue(ax, 0, 1, 1, 40);  setArrayValue(ax, 0, 1, 2, -0.331);
  setArrayValue(ax, 0, 2, 0, 13);  setArrayValue(ax, 0, 2, 1, 2);   setArrayValue(ax, 0, 2, 2, -0.731);
  setArrayValue(ax, 0, 3, 0, 62);  setArrayValue(ax, 0, 3, 1, 62);  setArrayValue(ax, 0, 3, 2, -0.893);
  setArrayValue(ax, 0, 4, 0, 28);  setArrayValue(ax, 0, 4, 1, 46);  setArrayValue(ax, 0, 4, 2, 0.283);
  setArrayValue(ax, 0, 5, 0, 15);  setArrayValue(ax, 0, 5, 1, 53);  setArrayValue(ax, 0, 5, 2, 0.940);
  setArrayValue(ax, 0, 6, 0, 49);  setArrayValue(ax, 0, 6, 1, 99);  setArrayValue(ax, 0, 6, 2, 0.231);
  setArrayValue(ax, 0, 7, 0, 88);  setArrayValue(ax, 0, 7, 1, 48);  setArrayValue(ax, 0, 7, 2, -1);
  setArrayValue(ax, 0, 8, 0, 15);  setArrayValue(ax, 0, 8, 1, 29);  setArrayValue(ax, 0, 8, 2, 0.88);
  setArrayValue(ax, 0, 9, 0, 63);  setArrayValue(ax, 0, 9, 1, 78);  setArrayValue(ax, 0, 9, 2, 0.321);

  setArrayValue(ax, 1, 0, 0, 4);   setArrayValue(ax, 1, 0, 1, 3);   setArrayValue(ax, 1, 0, 2, 0.30);
  setArrayValue(ax, 1, 1, 0, 2);   setArrayValue(ax, 1, 1, 1, 2);   setArrayValue(ax, 1, 1, 2, 0.12);
  setArrayValue(ax, 1, 2, 0, 5);   setArrayValue(ax, 1, 2, 1, 2);   setArrayValue(ax, 1, 2, 2, 0.53);
  setArrayValue(ax, 1, 3, 0, 1);   setArrayValue(ax, 1, 3, 1, 1);   setArrayValue(ax, 1, 3, 2, -1);
  setArrayValue(ax, 1, 4, 0, 7);   setArrayValue(ax, 1, 4, 1, 3);   setArrayValue(ax, 1, 4, 2, -13);
  setArrayValue(ax, 1, 5, 0, 5);   setArrayValue(ax, 1, 5, 1, 4);   setArrayValue(ax, 1, 5, 2, -56.9);
  setArrayValue(ax, 1, 6, 0, 12);  setArrayValue(ax, 1, 6, 1, 3);   setArrayValue(ax, 1, 6, 2, 0.36);
  setArrayValue(ax, 1, 7, 0, 3);   setArrayValue(ax, 1, 7, 1, 3);   setArrayValue(ax, 1, 7, 2, 0.4);
  setArrayValue(ax, 1, 8, 0, 26);  setArrayValue(ax, 1, 8, 1, 7);   setArrayValue(ax, 1, 8, 2, 0.1);
  setArrayValue(ax, 1, 9, 0, 5);   setArrayValue(ax, 1, 9, 1, 9);   setArrayValue(ax, 1, 9, 2, 0.37);


  NewArray(&ay, 2);

  NewArrayMatrix(&ay, 0, 10, 1);
  NewArrayMatrix(&ay, 1, 10, 1);

  setArrayValue(ay, 0, 0, 0, 50);
  setArrayValue(ay, 0, 1, 0, 86);
  setArrayValue(ay, 0, 2, 0, 20);
  setArrayValue(ay, 0, 3, 0, 95);
  setArrayValue(ay, 0, 4, 0, 61);
  setArrayValue(ay, 0, 5, 0, 50.5);
  setArrayValue(ay, 0, 6, 0, 113.5);
  setArrayValue(ay, 0, 7, 0, 119);
  setArrayValue(ay, 0, 8, 0, 62);
  setArrayValue(ay, 0, 9, 0, 116);

  setArrayValue(ay, 1, 0, 0, 50);
  setArrayValue(ay, 1, 1, 0, 86);
  setArrayValue(ay, 1, 2, 0, 20);
  setArrayValue(ay, 1, 3, 0, 95);
  setArrayValue(ay, 1, 4, 0, 61);
  setArrayValue(ay, 1, 5, 0, 50.5);
  setArrayValue(ay, 1, 6, 0, 113.5);
  setArrayValue(ay, 1, 7, 0, 119);
  setArrayValue(ay, 1, 8, 0, 62);
  setArrayValue(ay, 1, 9, 0, 116);

  initDVector(&q2x);
  initArray(&q2y);
  initArray(&sdep);

  NewUPLSModel(&m);

  ssignal run = SIGSCIENTIFICRUN;

  UPLS(ax, ay, 99, 1, 0, m, &run);

  /* we divide all the dataset in 5 groups */
  UPLSRandomGroupsCV(ax, ay, 1, 0, 3, 5, 20, &q2x, &q2y, &sdep, &m->predicted_y, &m->pred_residuals, &run);

  /*UPLSLOOCV(ax, ay, 1, 1, 3,&q2x, &q2y, &sdep, &m->predicted_y, &m->pred_residuals, &run);*/


  initArray(&yscrabl_q2);
  initArray(&yscrabl_sdep);
  UPLSYScrambling(ax, ay, 1, 1, 3, 3, 0, 0, 0, &yscrabl_q2, &yscrabl_sdep, NULL);

  puts("YSCRAMBL Q2");
  PrintArray(yscrabl_q2);
  puts("YSCRAMBL SDEP");
  PrintArray(yscrabl_sdep);

  puts("Q2Y");
  PrintArray(q2y);

  puts("SDEP");
  PrintArray(m->sdep);

  puts("Predicted Y");
  PrintArray(m->predicted_y);

  puts("RealY");
  PrintArray(ay);

  puts("Predicted Residuals");
  PrintArray(m->pred_residuals);

  if(q2y->order > 0){
    printf("CUTOFF: %u\n", (unsigned int)UPLSGetPCModelCutOff(q2y));
  }

  DelArray(&yscrabl_q2);
  DelArray(&yscrabl_sdep);
  DelUPLSModel(&m);
  DelArray(&sdep);
  DelDVector(&q2x);
  DelArray(&q2y);
  DelArray(&ay);
  DelArray(&ax);
}

void test8()
{
  array *ax;
  array *ay;

  dvector *q2x;
  array *q2y;
  array *sdep;
  array *predicted_y;

  puts(">>>>>>> Test 8: Compute Multi Way Random Group Cross Validation");

  NewArray(&ax, 2);

  NewArrayMatrix(&ax, 0, 10, 3);
  NewArrayMatrix(&ax, 1, 10, 3);

  setArrayValue(ax, 0, 0, 0, 37);  setArrayValue(ax, 0, 0, 1, 12);  setArrayValue(ax, 0, 0, 2, -1);
  setArrayValue(ax, 0, 1, 0, 62);  setArrayValue(ax, 0, 1, 1, 40);  setArrayValue(ax, 0, 1, 2, -0.331);
  setArrayValue(ax, 0, 2, 0, 13);  setArrayValue(ax, 0, 2, 1, 2);   setArrayValue(ax, 0, 2, 2, -0.731);
  setArrayValue(ax, 0, 3, 0, 62);  setArrayValue(ax, 0, 3, 1, 62);  setArrayValue(ax, 0, 3, 2, -0.893);
  setArrayValue(ax, 0, 4, 0, 28);  setArrayValue(ax, 0, 4, 1, 46);  setArrayValue(ax, 0, 4, 2, 0.283);
  setArrayValue(ax, 0, 5, 0, 15);  setArrayValue(ax, 0, 5, 1, 53);  setArrayValue(ax, 0, 5, 2, 0.940);
  setArrayValue(ax, 0, 6, 0, 49);  setArrayValue(ax, 0, 6, 1, 99);  setArrayValue(ax, 0, 6, 2, 0.231);
  setArrayValue(ax, 0, 7, 0, 88);  setArrayValue(ax, 0, 7, 1, 48);  setArrayValue(ax, 0, 7, 2, -1);
  setArrayValue(ax, 0, 8, 0, 15);  setArrayValue(ax, 0, 8, 1, 29);  setArrayValue(ax, 0, 8, 2, 0.88);
  setArrayValue(ax, 0, 9, 0, 63);  setArrayValue(ax, 0, 9, 1, 78);  setArrayValue(ax, 0, 9, 2, 0.321);

  setArrayValue(ax, 1, 0, 0, 4);   setArrayValue(ax, 1, 0, 1, 3);   setArrayValue(ax, 1, 0, 2, 0.30);
  setArrayValue(ax, 1, 1, 0, 2);   setArrayValue(ax, 1, 1, 1, 2);   setArrayValue(ax, 1, 1, 2, 0.12);
  setArrayValue(ax, 1, 2, 0, 5);   setArrayValue(ax, 1, 2, 1, 2);   setArrayValue(ax, 1, 2, 2, 0.53);
  setArrayValue(ax, 1, 3, 0, 1);   setArrayValue(ax, 1, 3, 1, 1);   setArrayValue(ax, 1, 3, 2, -1);
  setArrayValue(ax, 1, 4, 0, 7);   setArrayValue(ax, 1, 4, 1, 3);   setArrayValue(ax, 1, 4, 2, -13);
  setArrayValue(ax, 1, 5, 0, 5);   setArrayValue(ax, 1, 5, 1, 4);   setArrayValue(ax, 1, 5, 2, -56.9);
  setArrayValue(ax, 1, 6, 0, 12);  setArrayValue(ax, 1, 6, 1, 3);   setArrayValue(ax, 1, 6, 2, 0.36);
  setArrayValue(ax, 1, 7, 0, 3);   setArrayValue(ax, 1, 7, 1, 3);   setArrayValue(ax, 1, 7, 2, 0.4);
  setArrayValue(ax, 1, 8, 0, 26);  setArrayValue(ax, 1, 8, 1, 7);   setArrayValue(ax, 1, 8, 2, 0.1);
  setArrayValue(ax, 1, 9, 0, 5);   setArrayValue(ax, 1, 9, 1, 9);   setArrayValue(ax, 1, 9, 2, 0.37);


  NewArray(&ay, 2);

  NewArrayMatrix(&ay, 0, 10, 1);
  NewArrayMatrix(&ay, 1, 10, 1);

  setArrayValue(ay, 0, 0, 0, 50);
  setArrayValue(ay, 0, 1, 0, 86);
  setArrayValue(ay, 0, 2, 0, 20);
  setArrayValue(ay, 0, 3, 0, 95);
  setArrayValue(ay, 0, 4, 0, 61);
  setArrayValue(ay, 0, 5, 0, 50.5);
  setArrayValue(ay, 0, 6, 0, 113.5);
  setArrayValue(ay, 0, 7, 0, 119);
  setArrayValue(ay, 0, 8, 0, 62);
  setArrayValue(ay, 0, 9, 0, 116);

  setArrayValue(ay, 1, 0, 0, 50);
  setArrayValue(ay, 1, 1, 0, 86);
  setArrayValue(ay, 1, 2, 0, 20);
  setArrayValue(ay, 1, 3, 0, 95);
  setArrayValue(ay, 1, 4, 0, 61);
  setArrayValue(ay, 1, 5, 0, 50.5);
  setArrayValue(ay, 1, 6, 0, 113.5);
  setArrayValue(ay, 1, 7, 0, 119);
  setArrayValue(ay, 1, 8, 0, 62);
  setArrayValue(ay, 1, 9, 0, 116);

  initDVector(&q2x);
  initArray(&q2y);
  initArray(&sdep);
  initArray(&predicted_y);

  ssignal run = SIGSCIENTIFICRUN;
  /* we divide all the dataset in 5 groups */
  UPLSRandomGroupsCV(ax, ay, 1, 0, 3, 3, 20, &q2x, &q2y, &sdep, &predicted_y, NULL, &run);

  /*
  puts("Q2X");
  PrintDVector(q2x);
  */

  puts("Q2Y");
  PrintArray(q2y);

  puts("SDEP");
  PrintArray(sdep);

  puts("Predicted Y");
  PrintArray(predicted_y);

  puts("RealY");
  PrintArray(ay);

  if(q2y->order > 0){
    printf("CUTOFF: %u\n", (unsigned int)UPLSGetPCModelCutOff(q2y));
  }

  DelArray(&predicted_y);
  DelArray(&sdep);
  DelDVector(&q2x);
  DelArray(&q2y);
  DelArray(&ay);
  DelArray(&ax);
}

/*
 * Test from Fabio Broccatelli Calculation R(y)^2
 */
void test7()
{
  array *ax;
  array *ay;

  dvector *q2x;
  array *q2y;
  array *sdep;

  puts(">>>>>>> Test 7: Compute Multi Way Random Group Cross Validation");

  NewArray(&ax, 2);

  NewArrayMatrix(&ax, 0, 10, 2);
  NewArrayMatrix(&ax, 1, 10, 2);

  setArrayValue(ax, 0, 0, 0, 37);  setArrayValue(ax, 0, 0, 1, 12);
  setArrayValue(ax, 0, 1, 0, 62);  setArrayValue(ax, 0, 1, 1, 40);
  setArrayValue(ax, 0, 2, 0, 13);  setArrayValue(ax, 0, 2, 1, 2);
  setArrayValue(ax, 0, 3, 0, 62);  setArrayValue(ax, 0, 3, 1, 62);
  setArrayValue(ax, 0, 4, 0, 28);  setArrayValue(ax, 0, 4, 1, 46);
  setArrayValue(ax, 0, 5, 0, 15);  setArrayValue(ax, 0, 5, 1, 53);
  setArrayValue(ax, 0, 6, 0, 49);  setArrayValue(ax, 0, 6, 1, 99);
  setArrayValue(ax, 0, 7, 0, 88);  setArrayValue(ax, 0, 7, 1, 48);
  setArrayValue(ax, 0, 8, 0, 15);  setArrayValue(ax, 0, 8, 1, 29);
  setArrayValue(ax, 0, 9, 0, 63);  setArrayValue(ax, 0, 9, 1, 78);

  setArrayValue(ax, 1, 0, 0, 4);  setArrayValue(ax, 1, 0, 1, 3);
  setArrayValue(ax, 1, 1, 0, 2);  setArrayValue(ax, 1, 1, 1, 2);
  setArrayValue(ax, 1, 2, 0, 5);  setArrayValue(ax, 1, 2, 1, 2);
  setArrayValue(ax, 1, 3, 0, 1);  setArrayValue(ax, 1, 3, 1, 1);
  setArrayValue(ax, 1, 4, 0, 7);  setArrayValue(ax, 1, 4, 1, 3);
  setArrayValue(ax, 1, 5, 0, 5);  setArrayValue(ax, 1, 5, 1, 4);
  setArrayValue(ax, 1, 6, 0, 12);  setArrayValue(ax, 1, 6, 1, 3);
  setArrayValue(ax, 1, 7, 0, 3);  setArrayValue(ax, 1, 7, 1, 3);
  setArrayValue(ax, 1, 8, 0, 26);  setArrayValue(ax, 1, 8, 1, 7);
  setArrayValue(ax, 1, 9, 0, 5);  setArrayValue(ax, 1, 9, 1, 9);

  NewArray(&ay, 1);

  NewArrayMatrix(&ay, 0, 10, 1);

  setArrayValue(ay, 0, 0, 0, 50);
  setArrayValue(ay, 0, 1, 0, 86);
  setArrayValue(ay, 0, 2, 0, 20);
  setArrayValue(ay, 0, 3, 0, 95);
  setArrayValue(ay, 0, 4, 0, 61);
  setArrayValue(ay, 0, 5, 0, 50.5);
  setArrayValue(ay, 0, 6, 0, 113.5);
  setArrayValue(ay, 0, 7, 0, 119);
  setArrayValue(ay, 0, 8, 0, 62);
  setArrayValue(ay, 0, 9, 0, 116);



  puts("X:");
  PrintArray(ax);
  puts("Y:");
  PrintArray(ay);

  initDVector(&q2x);
  initArray(&q2y);
  initArray(&sdep);

  /* we divide all the dataset in 5 groups */
  UPLSRandomGroupsCV(ax, ay, 1, 0, 2, 5, 20, &q2x, &q2y, &sdep, NULL, NULL, NULL);

  puts("Q2X");
  PrintDVector(q2x);

  puts("Q2Y");
  PrintArray(q2y);

  puts("SDEP");
  PrintArray(sdep);

  DelArray(&sdep);
  DelDVector(&q2x);
  DelArray(&q2y);
  DelArray(&ay);
  DelArray(&ax);
}


/*
 * Test from Fabio Broccatelli Calculation R(y)^2
 */
void test6()
{
  array *ax;
  array *ay;

  array *axp;

  UPLSMODEL *m;

  array *py;
  matrix *pscores;

  dvector *r2x;
  array *r2y;
  array *sdec;


  puts(">>>>>>> Test 6: Compute Multi Way PLS and Prediction");

  NewArray(&ax, 2);

  NewArrayMatrix(&ax, 0, 3, 2);
  NewArrayMatrix(&ax, 1, 3, 2);

  setArrayValue(ax, 0, 0, 0, 37);  setArrayValue(ax, 0, 0, 1, 12);
  setArrayValue(ax, 0, 1, 0, 62);  setArrayValue(ax, 0, 1, 1, 40);
  setArrayValue(ax, 0, 2, 0, 13);  setArrayValue(ax, 0, 2, 1, 2);

  setArrayValue(ax, 1, 0, 0, 4);  setArrayValue(ax, 1, 0, 1, 3);
  setArrayValue(ax, 1, 1, 0, 2);  setArrayValue(ax, 1, 1, 1, 2);
  setArrayValue(ax, 1, 2, 0, 5);  setArrayValue(ax, 1, 2, 1, 2);


  NewArray(&ay, 2);

  NewArrayMatrix(&ay, 0, 3, 1);
  NewArrayMatrix(&ay, 1, 3, 1);

  setArrayValue(ay, 0, 0, 0, 50);
  setArrayValue(ay, 0, 1, 0, 86);
  setArrayValue(ay, 0, 2, 0, 20);


  setArrayValue(ay, 1, 0, 0, 50);
  setArrayValue(ay, 1, 1, 0, 86);
  setArrayValue(ay, 1, 2, 0, 20);


  NewArray(&axp, 2);
  NewArrayMatrix(&axp, 0, 1, 2);
  NewArrayMatrix(&axp, 1, 1, 2);

  setArrayValue(axp, 0, 0, 0, 62);  setArrayValue(axp, 0, 0, 1, 62);
  setArrayValue(axp, 1, 0, 0, 1);  setArrayValue(axp, 1, 0, 1, 1);


  NewUPLSModel(&m);

  /* Compute the model */
  ssignal run = SIGSCIENTIFICRUN;
  UPLS(ax, ay, 2, 1, 0, m, &run);


  initDVector(&r2x);
  initArray(&r2y);
  initArray(&sdec);

  /* calculating the r^2 for x and y model*/
  UPLSRSquared(ax, ay, m, 2, &r2x, &r2y, &sdec);


  puts("Data\nX:");
  PrintArray(ax);
  puts("Y:");
  PrintArray(ay);

  initArray(&py);
  initMatrix(&pscores);


  /*compute the xscore prediction*/
  UPLSScorePredictor(axp,m, 2, &pscores);

  /*compute the Y prediction*/
  UPLSYPredictor(pscores, m, 2, &py);


  PrintUPLSModel(m);

  puts("R^2 for Y");
  PrintArray(r2y);

  puts("R^2 for X");
  PrintDVector(r2x);

  puts("SDEC:");
  PrintArray(sdec);

  puts("PREDICTIONS");
  PrintArray(py);

  DelArray(&py);
  DelMatrix(&pscores);
  DelArray(&axp);

  DelArray(&sdec);
  DelDVector(&r2x);
  DelArray(&r2y);

  DelUPLSModel(&m);
  DelArray(&ay);
  DelArray(&ax);
}


/*
 * Test from Fabio Broccatelli
 */
void test5()
{
  array *ax;
  array *ay;

  array *axp;

  UPLSMODEL *m;

  array *py;
  matrix *pscores;

  puts(">>>>>>> Test 5: Compute Multi Way PLS and Prediction");

  NewArray(&ax, 2);

  NewArrayMatrix(&ax, 0, 3, 2);
  NewArrayMatrix(&ax, 1, 3, 2);

  setArrayValue(ax, 0, 0, 0, 37);  setArrayValue(ax, 0, 0, 1, 12);
  setArrayValue(ax, 0, 1, 0, 62);  setArrayValue(ax, 0, 1, 1, 40);
  setArrayValue(ax, 0, 2, 0, 13);  setArrayValue(ax, 0, 2, 1, 2);

  setArrayValue(ax, 1, 0, 0, 4);  setArrayValue(ax, 1, 0, 1, 3);
  setArrayValue(ax, 1, 1, 0, 2);  setArrayValue(ax, 1, 1, 1, 2);
  setArrayValue(ax, 1, 2, 0, 5);  setArrayValue(ax, 1, 2, 1, 2);


  NewArray(&ay, 2);

  NewArrayMatrix(&ay, 0, 3, 1);
  NewArrayMatrix(&ay, 1, 3, 1);

  setArrayValue(ay, 0, 0, 0, 50);
  setArrayValue(ay, 0, 1, 0, 86);
  setArrayValue(ay, 0, 2, 0, 20);


  setArrayValue(ay, 1, 0, 0, 50);
  setArrayValue(ay, 1, 1, 0, 86);
  setArrayValue(ay, 1, 2, 0, 20);


  NewArray(&axp, 2);
  NewArrayMatrix(&axp, 0, 1, 2);
  NewArrayMatrix(&axp, 1, 1, 2);

  setArrayValue(axp, 0, 0, 0, 62);  setArrayValue(axp, 0, 0, 1, 62);
  setArrayValue(axp, 1, 0, 0, 1);  setArrayValue(axp, 1, 0, 1, 1);


  NewUPLSModel(&m);

  UPLS(ax, ay, 2, 1, 0, m, NULL);

  puts("Data\nX:");
  PrintArray(ax);
  puts("Y:");
  PrintArray(ay);


  initArray(&py);
  initMatrix(&pscores);


  UPLSScorePredictor(axp, m, 2, &pscores);

  UPLSYPredictor(pscores, m, 2, &py);

  PrintUPLSModel(m);

  puts("Prediction Scores");
  PrintMatrix(pscores);

  puts("Prediction Y");
  PrintArray(py);

  DelArray(&py);
  DelMatrix(&pscores);
  DelArray(&axp);

  DelUPLSModel(&m);
  DelArray(&ay);
  DelArray(&ax);
}

void test4()
{
  array *ax;
  array *ay;

  array *axp;

  UPLSMODEL *m;

  array *py;
  matrix *pscores;

  puts(">>>>>>> Test 4: Compute Multi Way PLS and Prediction");

  NewArray(&ax, 2);

  NewArrayMatrix(&ax, 0, 3, 4);
  NewArrayMatrix(&ax, 1, 3, 4);

  setArrayValue(ax, 0, 0, 0, 0.424264);  setArrayValue(ax, 0, 0, 1, 0.565685);  setArrayValue(ax, 0, 0, 2, 0.4);  setArrayValue(ax, 0, 0, 3, 0.6);
  setArrayValue(ax, 0, 1, 0, 0.565685);  setArrayValue(ax, 0, 1, 1, 0.424264);  setArrayValue(ax, 0, 1, 2, 0.6);  setArrayValue(ax, 0, 1, 3, 0.4);
  setArrayValue(ax, 0, 2, 0, 0.707101);  setArrayValue(ax, 0, 2, 1, 0.707101);  setArrayValue(ax, 0, 2, 2, 0.1);  setArrayValue(ax, 0, 2, 3, 0.1);

  setArrayValue(ax, 1, 0, 0, 0.565685);  setArrayValue(ax, 1, 0, 1, 0.424264); setArrayValue(ax, 1, 0, 2, 0.6);  setArrayValue(ax, 1, 0, 3, 0.4);
  setArrayValue(ax, 1, 1, 0, 0.424264);  setArrayValue(ax, 1, 1, 1, 0.565685);  setArrayValue(ax, 1, 1, 2, 0.4);  setArrayValue(ax, 1, 1, 3, 0.6);
  setArrayValue(ax, 1, 2, 0, 0.707101);  setArrayValue(ax, 1, 2, 1, 0.707101);  setArrayValue(ax, 1, 2, 2, 0.1);  setArrayValue(ax, 1, 2, 3, 0.1);


  NewArray(&ay, 2);

  NewArrayMatrix(&ay, 0, 3, 1);
  NewArrayMatrix(&ay, 1, 3, 1);

  setArrayValue(ay, 0, 0, 0, 1.0);
  setArrayValue(ay, 0, 1, 0, 2.0);
  setArrayValue(ay, 0, 2, 0, 3.0);


  setArrayValue(ay, 1, 0, 0, 1.0);
  setArrayValue(ay, 1, 1, 0, 1.5);
  setArrayValue(ay, 1, 2, 0, 2.0);


  NewArray(&axp, 2);
  NewArrayMatrix(&axp, 0, 1, 4);
  NewArrayMatrix(&axp, 1, 1, 4);

  setArrayValue(axp, 0, 0, 0, 0.5);  setArrayValue(axp, 0, 0, 1, 0.6);  setArrayValue(axp, 0, 0, 2, 0.45);  setArrayValue(axp, 0, 0, 3, 0.55);
  setArrayValue(axp, 1, 0, 0, 0.6);  setArrayValue(axp, 1, 0, 1, 0.4);  setArrayValue(axp, 1, 0, 2, 0.55);  setArrayValue(axp, 1, 0, 3, 0.45);


  NewUPLSModel(&m);

  UPLS(ax, ay, 3, 1, 1, m, NULL);


  initArray(&py);
  initMatrix(&pscores);

  UPLSScorePredictor(axp, m, 3, &pscores);

  UPLSYPredictor(pscores, m, 3, &py);

  PrintUPLSModel(m);

  puts("Prediction Scores");
  PrintMatrix(pscores);

  puts("Prediction Y");
  PrintArray(py);

  DelArray(&py);
  DelMatrix(&pscores);
  DelArray(&axp);
  DelUPLSModel(&m);
  DelArray(&ay);
  DelArray(&ax);
}

/*Test from:
 *
 * PCA Nipals algorithm for Multi-Way Principal Component Analysis
 * Journal of Chemometrics, Vol 1, 41-56 (1987)
 * Svante Woold, Paul Geladi and Kim Esbensen
 */
void test3()
{
  array *ax;
  array *ay;

  array *axp;

  UPLSMODEL *m;


  array *py;
  matrix *pscores;

  puts(">>>>>>> Test 3: Compute Multi Way PLS and Prediction");

  NewArray(&ax, 2);

  NewArrayMatrix(&ax, 0, 3, 2);
  NewArrayMatrix(&ax, 1, 3, 2);

  setArrayValue(ax, 0, 0, 0, 0.424264);  setArrayValue(ax, 0, 0, 1, 0.565685);
  setArrayValue(ax, 0, 1, 0, 0.565685);  setArrayValue(ax, 0, 1, 1, 0.424264);
  setArrayValue(ax, 0, 2, 0, 0.707101);  setArrayValue(ax, 0, 2, 1, 0.707101);

  setArrayValue(ax, 1, 0, 0, 0.565685);  setArrayValue(ax, 1, 0, 1, 0.424264);
  setArrayValue(ax, 1, 1, 0, 0.424264);  setArrayValue(ax, 1, 1, 1, 0.565685);
  setArrayValue(ax, 1, 2, 0, 0.707101);  setArrayValue(ax, 1, 2, 1, 0.707101);


  NewArray(&ay, 2);

  NewArrayMatrix(&ay, 0, 3, 1);
  NewArrayMatrix(&ay, 1, 3, 1);

  setArrayValue(ay, 0, 0, 0, 1.0);
  setArrayValue(ay, 0, 1, 0, 2.0);
  setArrayValue(ay, 0, 2, 0, 3.0);


  setArrayValue(ay, 1, 0, 0, 1.0);
  setArrayValue(ay, 1, 1, 0, 1.5);
  setArrayValue(ay, 1, 2, 0, 2.0);


  NewArray(&axp, 2);
  NewArrayMatrix(&axp, 0, 1, 2);
  NewArrayMatrix(&axp, 1, 1, 2);

  setArrayValue(axp, 0, 0, 0, 0.5);  setArrayValue(axp, 0, 0, 1, 0.6);
  setArrayValue(axp, 1, 0, 0, 0.6);  setArrayValue(axp, 1, 0, 1, 0.4);


  NewUPLSModel(&m);


  UPLS(ax, ay, 3, 1, 0, m, NULL);


  initArray(&py);
  initMatrix(&pscores);

  UPLSScorePredictor(axp, m, 3, &pscores);

  UPLSYPredictor(pscores, m, 3, &py);

  PrintUPLSModel(m);

  puts("Prediction Scores");
  PrintMatrix(pscores);

  puts("Prediction Y");
  PrintArray(py);

  DelArray(&py);
  DelMatrix(&pscores);
  DelArray(&axp);

  DelUPLSModel(&m);
  DelArray(&ay);
  DelArray(&ax);
}



/*
 * Test with Random numbers
 */
void test2()
{
  size_t i, j, k;
  array *ax;
  array *ay;

  UPLSMODEL *m;

  puts(">>>>>>> Test 2: Compute Multi Way PLS with random variables");

  NewArray(&ax, 3);

  for(i = 0; i < ax->order; i++)
    NewArrayMatrix(&ax, i, 30, 10);

  srand(time(0));
  for(i = 0; i < ax->order; i++){
    for(j = 0; j < ax->m[i]->row; j++){
      for(k = 0; k < ax->m[i]->col; k++)
          setArrayValue(ax, i, j, k, rand()/((double)(RAND_MAX)+1));
    }
  }



  NewArray(&ay, 1);
  NewArrayMatrix(&ay, 0, 30, 1);

  srand(time(0));
  for(i = 0; i < ay->order; i++){
    for(j = 0; j < ay->m[i]->row; j++){
      for(k = 0; k < ay->m[i]->col; k++)
          setArrayValue(ay, i, j, k, rand()/((double)(RAND_MAX)+1));
    }
  }

  NewUPLSModel(&m);

  UPLS(ax, ay, 3, 1, 0, m, NULL);

  PrintUPLSModel(m);

  DelUPLSModel(&m);
  DelArray(&ay);
  DelArray(&ax);
}




/*Test from:
 *
 * PCA Nipals algorithm for Multi-Way Principal Component Analysis
 * Journal of Chemometrics, Vol 1, 41-56 (1987)
 * Svante Woold, Paul Geladi and Kim Esbensen
 */
void test1()
{
  array *ax;
  array *ay;

  UPLSMODEL *m;

  puts(">>>>>>> Test 1: Compute Multi Way PLS");

  NewArray(&ax, 2);

  NewArrayMatrix(&ax, 0, 3, 2);
  NewArrayMatrix(&ax, 1, 3, 2);

  setArrayValue(ax, 0, 0, 0, 0.424264);  setArrayValue(ax, 0, 0, 1, 0.565685);
  setArrayValue(ax, 0, 1, 0, 0.565685);  setArrayValue(ax, 0, 1, 1, 0.424264);
  setArrayValue(ax, 0, 2, 0, 0.707101);  setArrayValue(ax, 0, 2, 1, 0.707101);

  setArrayValue(ax, 1, 0, 0, 0.565685);  setArrayValue(ax, 1, 0, 1, 0.424264);
  setArrayValue(ax, 1, 1, 0, 0.424264);  setArrayValue(ax, 1, 1, 1, 0.565685);
  setArrayValue(ax, 1, 2, 0, 0.707101);  setArrayValue(ax, 1, 2, 1, 0.707101);


  NewArray(&ay, 2);

  NewArrayMatrix(&ay, 0, 3, 1);
  NewArrayMatrix(&ay, 1, 3, 1);

  setArrayValue(ay, 0, 0, 0, 1.0);
  setArrayValue(ay, 0, 1, 0, 2.0);
  setArrayValue(ay, 0, 2, 0, 3.0);


  setArrayValue(ay, 1, 0, 0, 1.0);
  setArrayValue(ay, 1, 1, 0, 1.5);
  setArrayValue(ay, 1, 2, 0, 2.0);


  NewUPLSModel(&m);

  UPLS(ax, ay, 3, 1, 1, m, NULL);

  PrintUPLSModel(m);

  DelUPLSModel(&m);
  DelArray(&ay);
  DelArray(&ax);
}



int main(void)
{
  /* test 1-5 */
  test1();
  test2();
  test3();
  test4();
  test5();
  /*test 6-9 */
  test6();
  test7();
  test8();
  test9();
  return 0;
}

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
#include "tensor.h"
#include "upls.h"
#include "upca.h"

void test9()
{
  tensor *ax;
  tensor *ay;

  dvector *q2x;
  tensor *q2y;
  tensor *sdep;
  tensor *yscrabl_q2, *yscrabl_sdep;

  UPLSMODEL *m;

  puts(">>>>>>> Test 9: Compute Multi Way LOO Cross Validation");

  NewTensor(&ax, 2);

  NewTensorMatrix(ax, 0, 10, 3);
  NewTensorMatrix(ax, 1, 10, 3);

  setTensorValue(ax, 0, 0, 0, 37);  setTensorValue(ax, 0, 0, 1, 12);  setTensorValue(ax, 0, 0, 2, -1);
  setTensorValue(ax, 0, 1, 0, 62);  setTensorValue(ax, 0, 1, 1, 40);  setTensorValue(ax, 0, 1, 2, -0.331);
  setTensorValue(ax, 0, 2, 0, 13);  setTensorValue(ax, 0, 2, 1, 2);   setTensorValue(ax, 0, 2, 2, -0.731);
  setTensorValue(ax, 0, 3, 0, 62);  setTensorValue(ax, 0, 3, 1, 62);  setTensorValue(ax, 0, 3, 2, -0.893);
  setTensorValue(ax, 0, 4, 0, 28);  setTensorValue(ax, 0, 4, 1, 46);  setTensorValue(ax, 0, 4, 2, 0.283);
  setTensorValue(ax, 0, 5, 0, 15);  setTensorValue(ax, 0, 5, 1, 53);  setTensorValue(ax, 0, 5, 2, 0.940);
  setTensorValue(ax, 0, 6, 0, 49);  setTensorValue(ax, 0, 6, 1, 99);  setTensorValue(ax, 0, 6, 2, 0.231);
  setTensorValue(ax, 0, 7, 0, 88);  setTensorValue(ax, 0, 7, 1, 48);  setTensorValue(ax, 0, 7, 2, -1);
  setTensorValue(ax, 0, 8, 0, 15);  setTensorValue(ax, 0, 8, 1, 29);  setTensorValue(ax, 0, 8, 2, 0.88);
  setTensorValue(ax, 0, 9, 0, 63);  setTensorValue(ax, 0, 9, 1, 78);  setTensorValue(ax, 0, 9, 2, 0.321);

  setTensorValue(ax, 1, 0, 0, 4);   setTensorValue(ax, 1, 0, 1, 3);   setTensorValue(ax, 1, 0, 2, 0.30);
  setTensorValue(ax, 1, 1, 0, 2);   setTensorValue(ax, 1, 1, 1, 2);   setTensorValue(ax, 1, 1, 2, 0.12);
  setTensorValue(ax, 1, 2, 0, 5);   setTensorValue(ax, 1, 2, 1, 2);   setTensorValue(ax, 1, 2, 2, 0.53);
  setTensorValue(ax, 1, 3, 0, 1);   setTensorValue(ax, 1, 3, 1, 1);   setTensorValue(ax, 1, 3, 2, -1);
  setTensorValue(ax, 1, 4, 0, 7);   setTensorValue(ax, 1, 4, 1, 3);   setTensorValue(ax, 1, 4, 2, -13);
  setTensorValue(ax, 1, 5, 0, 5);   setTensorValue(ax, 1, 5, 1, 4);   setTensorValue(ax, 1, 5, 2, -56.9);
  setTensorValue(ax, 1, 6, 0, 12);  setTensorValue(ax, 1, 6, 1, 3);   setTensorValue(ax, 1, 6, 2, 0.36);
  setTensorValue(ax, 1, 7, 0, 3);   setTensorValue(ax, 1, 7, 1, 3);   setTensorValue(ax, 1, 7, 2, 0.4);
  setTensorValue(ax, 1, 8, 0, 26);  setTensorValue(ax, 1, 8, 1, 7);   setTensorValue(ax, 1, 8, 2, 0.1);
  setTensorValue(ax, 1, 9, 0, 5);   setTensorValue(ax, 1, 9, 1, 9);   setTensorValue(ax, 1, 9, 2, 0.37);


  NewTensor(&ay, 2);

  NewTensorMatrix(ay, 0, 10, 1);
  NewTensorMatrix(ay, 1, 10, 1);

  setTensorValue(ay, 0, 0, 0, 50);
  setTensorValue(ay, 0, 1, 0, 86);
  setTensorValue(ay, 0, 2, 0, 20);
  setTensorValue(ay, 0, 3, 0, 95);
  setTensorValue(ay, 0, 4, 0, 61);
  setTensorValue(ay, 0, 5, 0, 50.5);
  setTensorValue(ay, 0, 6, 0, 113.5);
  setTensorValue(ay, 0, 7, 0, 119);
  setTensorValue(ay, 0, 8, 0, 62);
  setTensorValue(ay, 0, 9, 0, 116);

  setTensorValue(ay, 1, 0, 0, 50);
  setTensorValue(ay, 1, 1, 0, 86);
  setTensorValue(ay, 1, 2, 0, 20);
  setTensorValue(ay, 1, 3, 0, 95);
  setTensorValue(ay, 1, 4, 0, 61);
  setTensorValue(ay, 1, 5, 0, 50.5);
  setTensorValue(ay, 1, 6, 0, 113.5);
  setTensorValue(ay, 1, 7, 0, 119);
  setTensorValue(ay, 1, 8, 0, 62);
  setTensorValue(ay, 1, 9, 0, 116);

  initDVector(&q2x);
  initTensor(&q2y);
  initTensor(&sdep);

  NewUPLSModel(&m);

  ssignal run = SIGSCIENTIFICRUN;

  UPLS(ax, ay, 99, 1, 0, m, &run);

  /* we divide all the dataset in 5 groups */
  UPLSRandomGroupsCV(ax, ay, 1, 0, 3, 5, 20, &q2x, &q2y, &sdep, &m->predicted_y, &m->pred_residuals, &run);

  /*UPLSLOOCV(ax, ay, 1, 1, 3,&q2x, &q2y, &sdep, &m->predicted_y, &m->pred_residuals, &run);*/


  initTensor(&yscrabl_q2);
  initTensor(&yscrabl_sdep);
  UPLSYScrambling(ax, ay, 1, 1, 3, 3, 0, 0, 0, &yscrabl_q2, &yscrabl_sdep, NULL);

  puts("YSCRAMBL Q2");
  PrintTensor(yscrabl_q2);
  puts("YSCRAMBL SDEP");
  PrintTensor(yscrabl_sdep);

  puts("Q2Y");
  PrintTensor(q2y);

  puts("SDEP");
  PrintTensor(m->sdep);

  puts("Predicted Y");
  PrintTensor(m->predicted_y);

  puts("RealY");
  PrintTensor(ay);

  puts("Predicted Residuals");
  PrintTensor(m->pred_residuals);

  if(q2y->order > 0){
    printf("CUTOFF: %u\n", (unsigned int)UPLSGetPCModelCutOff(q2y));
  }

  DelTensor(&yscrabl_q2);
  DelTensor(&yscrabl_sdep);
  DelUPLSModel(&m);
  DelTensor(&sdep);
  DelDVector(&q2x);
  DelTensor(&q2y);
  DelTensor(&ay);
  DelTensor(&ax);
}

void test8()
{
  tensor *ax;
  tensor *ay;

  dvector *q2x;
  tensor *q2y;
  tensor *sdep;
  tensor *predicted_y;

  puts(">>>>>>> Test 8: Compute Multi Way Random Group Cross Validation");

  NewTensor(&ax, 2);

  NewTensorMatrix(ax, 0, 10, 3);
  NewTensorMatrix(ax, 1, 10, 3);

  setTensorValue(ax, 0, 0, 0, 37);  setTensorValue(ax, 0, 0, 1, 12);  setTensorValue(ax, 0, 0, 2, -1);
  setTensorValue(ax, 0, 1, 0, 62);  setTensorValue(ax, 0, 1, 1, 40);  setTensorValue(ax, 0, 1, 2, -0.331);
  setTensorValue(ax, 0, 2, 0, 13);  setTensorValue(ax, 0, 2, 1, 2);   setTensorValue(ax, 0, 2, 2, -0.731);
  setTensorValue(ax, 0, 3, 0, 62);  setTensorValue(ax, 0, 3, 1, 62);  setTensorValue(ax, 0, 3, 2, -0.893);
  setTensorValue(ax, 0, 4, 0, 28);  setTensorValue(ax, 0, 4, 1, 46);  setTensorValue(ax, 0, 4, 2, 0.283);
  setTensorValue(ax, 0, 5, 0, 15);  setTensorValue(ax, 0, 5, 1, 53);  setTensorValue(ax, 0, 5, 2, 0.940);
  setTensorValue(ax, 0, 6, 0, 49);  setTensorValue(ax, 0, 6, 1, 99);  setTensorValue(ax, 0, 6, 2, 0.231);
  setTensorValue(ax, 0, 7, 0, 88);  setTensorValue(ax, 0, 7, 1, 48);  setTensorValue(ax, 0, 7, 2, -1);
  setTensorValue(ax, 0, 8, 0, 15);  setTensorValue(ax, 0, 8, 1, 29);  setTensorValue(ax, 0, 8, 2, 0.88);
  setTensorValue(ax, 0, 9, 0, 63);  setTensorValue(ax, 0, 9, 1, 78);  setTensorValue(ax, 0, 9, 2, 0.321);

  setTensorValue(ax, 1, 0, 0, 4);   setTensorValue(ax, 1, 0, 1, 3);   setTensorValue(ax, 1, 0, 2, 0.30);
  setTensorValue(ax, 1, 1, 0, 2);   setTensorValue(ax, 1, 1, 1, 2);   setTensorValue(ax, 1, 1, 2, 0.12);
  setTensorValue(ax, 1, 2, 0, 5);   setTensorValue(ax, 1, 2, 1, 2);   setTensorValue(ax, 1, 2, 2, 0.53);
  setTensorValue(ax, 1, 3, 0, 1);   setTensorValue(ax, 1, 3, 1, 1);   setTensorValue(ax, 1, 3, 2, -1);
  setTensorValue(ax, 1, 4, 0, 7);   setTensorValue(ax, 1, 4, 1, 3);   setTensorValue(ax, 1, 4, 2, -13);
  setTensorValue(ax, 1, 5, 0, 5);   setTensorValue(ax, 1, 5, 1, 4);   setTensorValue(ax, 1, 5, 2, -56.9);
  setTensorValue(ax, 1, 6, 0, 12);  setTensorValue(ax, 1, 6, 1, 3);   setTensorValue(ax, 1, 6, 2, 0.36);
  setTensorValue(ax, 1, 7, 0, 3);   setTensorValue(ax, 1, 7, 1, 3);   setTensorValue(ax, 1, 7, 2, 0.4);
  setTensorValue(ax, 1, 8, 0, 26);  setTensorValue(ax, 1, 8, 1, 7);   setTensorValue(ax, 1, 8, 2, 0.1);
  setTensorValue(ax, 1, 9, 0, 5);   setTensorValue(ax, 1, 9, 1, 9);   setTensorValue(ax, 1, 9, 2, 0.37);


  NewTensor(&ay, 2);

  NewTensorMatrix(ay, 0, 10, 1);
  NewTensorMatrix(ay, 1, 10, 1);

  setTensorValue(ay, 0, 0, 0, 50);
  setTensorValue(ay, 0, 1, 0, 86);
  setTensorValue(ay, 0, 2, 0, 20);
  setTensorValue(ay, 0, 3, 0, 95);
  setTensorValue(ay, 0, 4, 0, 61);
  setTensorValue(ay, 0, 5, 0, 50.5);
  setTensorValue(ay, 0, 6, 0, 113.5);
  setTensorValue(ay, 0, 7, 0, 119);
  setTensorValue(ay, 0, 8, 0, 62);
  setTensorValue(ay, 0, 9, 0, 116);

  setTensorValue(ay, 1, 0, 0, 50);
  setTensorValue(ay, 1, 1, 0, 86);
  setTensorValue(ay, 1, 2, 0, 20);
  setTensorValue(ay, 1, 3, 0, 95);
  setTensorValue(ay, 1, 4, 0, 61);
  setTensorValue(ay, 1, 5, 0, 50.5);
  setTensorValue(ay, 1, 6, 0, 113.5);
  setTensorValue(ay, 1, 7, 0, 119);
  setTensorValue(ay, 1, 8, 0, 62);
  setTensorValue(ay, 1, 9, 0, 116);

  initDVector(&q2x);
  initTensor(&q2y);
  initTensor(&sdep);
  initTensor(&predicted_y);

  ssignal run = SIGSCIENTIFICRUN;
  /* we divide all the dataset in 5 groups */
  UPLSRandomGroupsCV(ax, ay, 1, 0, 3, 3, 20, &q2x, &q2y, &sdep, &predicted_y, NULL, &run);

  /*
  puts("Q2X");
  PrintDVector(q2x);
  */

  puts("Q2Y");
  PrintTensor(q2y);

  puts("SDEP");
  PrintTensor(sdep);

  puts("Predicted Y");
  PrintTensor(predicted_y);

  puts("RealY");
  PrintTensor(ay);

  if(q2y->order > 0){
    printf("CUTOFF: %u\n", (unsigned int)UPLSGetPCModelCutOff(q2y));
  }

  DelTensor(&predicted_y);
  DelTensor(&sdep);
  DelDVector(&q2x);
  DelTensor(&q2y);
  DelTensor(&ay);
  DelTensor(&ax);
}

/*
 * Test from Fabio Broccatelli Calculation R(y)^2
 */
void test7()
{
  tensor *ax;
  tensor *ay;

  dvector *q2x;
  tensor *q2y;
  tensor *sdep;

  puts(">>>>>>> Test 7: Compute Multi Way Random Group Cross Validation");

  NewTensor(&ax, 2);

  NewTensorMatrix(ax, 0, 10, 2);
  NewTensorMatrix(ax, 1, 10, 2);

  setTensorValue(ax, 0, 0, 0, 37);  setTensorValue(ax, 0, 0, 1, 12);
  setTensorValue(ax, 0, 1, 0, 62);  setTensorValue(ax, 0, 1, 1, 40);
  setTensorValue(ax, 0, 2, 0, 13);  setTensorValue(ax, 0, 2, 1, 2);
  setTensorValue(ax, 0, 3, 0, 62);  setTensorValue(ax, 0, 3, 1, 62);
  setTensorValue(ax, 0, 4, 0, 28);  setTensorValue(ax, 0, 4, 1, 46);
  setTensorValue(ax, 0, 5, 0, 15);  setTensorValue(ax, 0, 5, 1, 53);
  setTensorValue(ax, 0, 6, 0, 49);  setTensorValue(ax, 0, 6, 1, 99);
  setTensorValue(ax, 0, 7, 0, 88);  setTensorValue(ax, 0, 7, 1, 48);
  setTensorValue(ax, 0, 8, 0, 15);  setTensorValue(ax, 0, 8, 1, 29);
  setTensorValue(ax, 0, 9, 0, 63);  setTensorValue(ax, 0, 9, 1, 78);

  setTensorValue(ax, 1, 0, 0, 4);  setTensorValue(ax, 1, 0, 1, 3);
  setTensorValue(ax, 1, 1, 0, 2);  setTensorValue(ax, 1, 1, 1, 2);
  setTensorValue(ax, 1, 2, 0, 5);  setTensorValue(ax, 1, 2, 1, 2);
  setTensorValue(ax, 1, 3, 0, 1);  setTensorValue(ax, 1, 3, 1, 1);
  setTensorValue(ax, 1, 4, 0, 7);  setTensorValue(ax, 1, 4, 1, 3);
  setTensorValue(ax, 1, 5, 0, 5);  setTensorValue(ax, 1, 5, 1, 4);
  setTensorValue(ax, 1, 6, 0, 12);  setTensorValue(ax, 1, 6, 1, 3);
  setTensorValue(ax, 1, 7, 0, 3);  setTensorValue(ax, 1, 7, 1, 3);
  setTensorValue(ax, 1, 8, 0, 26);  setTensorValue(ax, 1, 8, 1, 7);
  setTensorValue(ax, 1, 9, 0, 5);  setTensorValue(ax, 1, 9, 1, 9);

  NewTensor(&ay, 1);

  NewTensorMatrix(ay, 0, 10, 1);

  setTensorValue(ay, 0, 0, 0, 50);
  setTensorValue(ay, 0, 1, 0, 86);
  setTensorValue(ay, 0, 2, 0, 20);
  setTensorValue(ay, 0, 3, 0, 95);
  setTensorValue(ay, 0, 4, 0, 61);
  setTensorValue(ay, 0, 5, 0, 50.5);
  setTensorValue(ay, 0, 6, 0, 113.5);
  setTensorValue(ay, 0, 7, 0, 119);
  setTensorValue(ay, 0, 8, 0, 62);
  setTensorValue(ay, 0, 9, 0, 116);



  puts("X:");
  PrintTensor(ax);
  puts("Y:");
  PrintTensor(ay);

  initDVector(&q2x);
  initTensor(&q2y);
  initTensor(&sdep);

  /* we divide all the dataset in 5 groups */
  UPLSRandomGroupsCV(ax, ay, 1, 0, 2, 5, 20, &q2x, &q2y, &sdep, NULL, NULL, NULL);

  puts("Q2X");
  PrintDVector(q2x);

  puts("Q2Y");
  PrintTensor(q2y);

  puts("SDEP");
  PrintTensor(sdep);

  DelTensor(&sdep);
  DelDVector(&q2x);
  DelTensor(&q2y);
  DelTensor(&ay);
  DelTensor(&ax);
}


/*
 * Test from Fabio Broccatelli Calculation R(y)^2
 */
void test6()
{
  tensor *ax;
  tensor *ay;

  tensor *axp;

  UPLSMODEL *m;

  tensor *py;
  matrix *pscores;

  dvector *r2x;
  tensor *r2y;
  tensor *sdec;


  puts(">>>>>>> Test 6: Compute Multi Way PLS and Prediction");

  NewTensor(&ax, 2);

  NewTensorMatrix(ax, 0, 3, 2);
  NewTensorMatrix(ax, 1, 3, 2);

  setTensorValue(ax, 0, 0, 0, 37);  setTensorValue(ax, 0, 0, 1, 12);
  setTensorValue(ax, 0, 1, 0, 62);  setTensorValue(ax, 0, 1, 1, 40);
  setTensorValue(ax, 0, 2, 0, 13);  setTensorValue(ax, 0, 2, 1, 2);

  setTensorValue(ax, 1, 0, 0, 4);  setTensorValue(ax, 1, 0, 1, 3);
  setTensorValue(ax, 1, 1, 0, 2);  setTensorValue(ax, 1, 1, 1, 2);
  setTensorValue(ax, 1, 2, 0, 5);  setTensorValue(ax, 1, 2, 1, 2);


  NewTensor(&ay, 2);

  NewTensorMatrix(ay, 0, 3, 1);
  NewTensorMatrix(ay, 1, 3, 1);

  setTensorValue(ay, 0, 0, 0, 50);
  setTensorValue(ay, 0, 1, 0, 86);
  setTensorValue(ay, 0, 2, 0, 20);


  setTensorValue(ay, 1, 0, 0, 50);
  setTensorValue(ay, 1, 1, 0, 86);
  setTensorValue(ay, 1, 2, 0, 20);


  NewTensor(&axp, 2);
  NewTensorMatrix(axp, 0, 1, 2);
  NewTensorMatrix(axp, 1, 1, 2);

  setTensorValue(axp, 0, 0, 0, 62);  setTensorValue(axp, 0, 0, 1, 62);
  setTensorValue(axp, 1, 0, 0, 1);  setTensorValue(axp, 1, 0, 1, 1);


  NewUPLSModel(&m);

  /* Compute the model */
  ssignal run = SIGSCIENTIFICRUN;
  UPLS(ax, ay, 2, 1, 0, m, &run);


  initDVector(&r2x);
  initTensor(&r2y);
  initTensor(&sdec);

  /* calculating the r^2 for x and y model*/
  UPLSRSquared(ax, ay, m, 2, r2x, r2y, sdec);


  puts("Data\nX:");
  PrintTensor(ax);
  puts("Y:");
  PrintTensor(ay);

  initTensor(&py);
  initMatrix(&pscores);


  /*compute the xscore prediction*/
  UPLSScorePredictor(axp,m, 2, pscores);

  /*compute the Y prediction*/
  UPLSYPredictor(pscores, m, 2, py);


  PrintUPLSModel(m);

  puts("R^2 for Y");
  PrintTensor(r2y);

  puts("R^2 for X");
  PrintDVector(r2x);

  puts("SDEC:");
  PrintTensor(sdec);

  puts("PREDICTIONS");
  PrintTensor(py);

  DelTensor(&py);
  DelMatrix(&pscores);
  DelTensor(&axp);

  DelTensor(&sdec);
  DelDVector(&r2x);
  DelTensor(&r2y);

  DelUPLSModel(&m);
  DelTensor(&ay);
  DelTensor(&ax);
}


/*
 * Test from Fabio Broccatelli
 */
void test5()
{
  tensor *ax;
  tensor *ay;

  tensor *axp;

  UPLSMODEL *m;

  tensor *py;
  matrix *pscores;

  puts(">>>>>>> Test 5: Compute Multi Way PLS and Prediction");

  NewTensor(&ax, 2);

  NewTensorMatrix(ax, 0, 3, 2);
  NewTensorMatrix(ax, 1, 3, 2);

  setTensorValue(ax, 0, 0, 0, 37);  setTensorValue(ax, 0, 0, 1, 12);
  setTensorValue(ax, 0, 1, 0, 62);  setTensorValue(ax, 0, 1, 1, 40);
  setTensorValue(ax, 0, 2, 0, 13);  setTensorValue(ax, 0, 2, 1, 2);

  setTensorValue(ax, 1, 0, 0, 4);  setTensorValue(ax, 1, 0, 1, 3);
  setTensorValue(ax, 1, 1, 0, 2);  setTensorValue(ax, 1, 1, 1, 2);
  setTensorValue(ax, 1, 2, 0, 5);  setTensorValue(ax, 1, 2, 1, 2);


  NewTensor(&ay, 2);

  NewTensorMatrix(ay, 0, 3, 1);
  NewTensorMatrix(ay, 1, 3, 1);

  setTensorValue(ay, 0, 0, 0, 50);
  setTensorValue(ay, 0, 1, 0, 86);
  setTensorValue(ay, 0, 2, 0, 20);


  setTensorValue(ay, 1, 0, 0, 50);
  setTensorValue(ay, 1, 1, 0, 86);
  setTensorValue(ay, 1, 2, 0, 20);


  NewTensor(&axp, 2);
  NewTensorMatrix(axp, 0, 1, 2);
  NewTensorMatrix(axp, 1, 1, 2);

  setTensorValue(axp, 0, 0, 0, 62);  setTensorValue(axp, 0, 0, 1, 62);
  setTensorValue(axp, 1, 0, 0, 1);  setTensorValue(axp, 1, 0, 1, 1);


  NewUPLSModel(&m);

  UPLS(ax, ay, 2, 1, 0, m, NULL);

  puts("Data\nX:");
  PrintTensor(ax);
  puts("Y:");
  PrintTensor(ay);


  initTensor(&py);
  initMatrix(&pscores);


  UPLSScorePredictor(axp, m, 2, pscores);

  UPLSYPredictor(pscores, m, 2, py);

  PrintUPLSModel(m);

  puts("Prediction Scores");
  PrintMatrix(pscores);

  puts("Prediction Y");
  PrintTensor(py);

  DelTensor(&py);
  DelMatrix(&pscores);
  DelTensor(&axp);

  DelUPLSModel(&m);
  DelTensor(&ay);
  DelTensor(&ax);
}

void test4()
{
  tensor *ax;
  tensor *ay;

  tensor *axp;

  UPLSMODEL *m;

  tensor *py;
  matrix *pscores;

  puts(">>>>>>> Test 4: Compute Multi Way PLS and Prediction");

  NewTensor(&ax, 2);

  NewTensorMatrix(ax, 0, 3, 4);
  NewTensorMatrix(ax, 1, 3, 4);

  setTensorValue(ax, 0, 0, 0, 0.424264);  setTensorValue(ax, 0, 0, 1, 0.565685);  setTensorValue(ax, 0, 0, 2, 0.4);  setTensorValue(ax, 0, 0, 3, 0.6);
  setTensorValue(ax, 0, 1, 0, 0.565685);  setTensorValue(ax, 0, 1, 1, 0.424264);  setTensorValue(ax, 0, 1, 2, 0.6);  setTensorValue(ax, 0, 1, 3, 0.4);
  setTensorValue(ax, 0, 2, 0, 0.707101);  setTensorValue(ax, 0, 2, 1, 0.707101);  setTensorValue(ax, 0, 2, 2, 0.1);  setTensorValue(ax, 0, 2, 3, 0.1);

  setTensorValue(ax, 1, 0, 0, 0.565685);  setTensorValue(ax, 1, 0, 1, 0.424264); setTensorValue(ax, 1, 0, 2, 0.6);  setTensorValue(ax, 1, 0, 3, 0.4);
  setTensorValue(ax, 1, 1, 0, 0.424264);  setTensorValue(ax, 1, 1, 1, 0.565685);  setTensorValue(ax, 1, 1, 2, 0.4);  setTensorValue(ax, 1, 1, 3, 0.6);
  setTensorValue(ax, 1, 2, 0, 0.707101);  setTensorValue(ax, 1, 2, 1, 0.707101);  setTensorValue(ax, 1, 2, 2, 0.1);  setTensorValue(ax, 1, 2, 3, 0.1);


  NewTensor(&ay, 2);

  NewTensorMatrix(ay, 0, 3, 1);
  NewTensorMatrix(ay, 1, 3, 1);

  setTensorValue(ay, 0, 0, 0, 1.0);
  setTensorValue(ay, 0, 1, 0, 2.0);
  setTensorValue(ay, 0, 2, 0, 3.0);


  setTensorValue(ay, 1, 0, 0, 1.0);
  setTensorValue(ay, 1, 1, 0, 1.5);
  setTensorValue(ay, 1, 2, 0, 2.0);


  NewTensor(&axp, 2);
  NewTensorMatrix(axp, 0, 1, 4);
  NewTensorMatrix(axp, 1, 1, 4);

  setTensorValue(axp, 0, 0, 0, 0.5);  setTensorValue(axp, 0, 0, 1, 0.6);  setTensorValue(axp, 0, 0, 2, 0.45);  setTensorValue(axp, 0, 0, 3, 0.55);
  setTensorValue(axp, 1, 0, 0, 0.6);  setTensorValue(axp, 1, 0, 1, 0.4);  setTensorValue(axp, 1, 0, 2, 0.55);  setTensorValue(axp, 1, 0, 3, 0.45);


  NewUPLSModel(&m);

  UPLS(ax, ay, 3, 1, 1, m, NULL);


  initTensor(&py);
  initMatrix(&pscores);

  UPLSScorePredictor(axp, m, 3, pscores);

  UPLSYPredictor(pscores, m, 3, py);

  PrintUPLSModel(m);

  puts("Prediction Scores");
  PrintMatrix(pscores);

  puts("Prediction Y");
  PrintTensor(py);

  DelTensor(&py);
  DelMatrix(&pscores);
  DelTensor(&axp);
  DelUPLSModel(&m);
  DelTensor(&ay);
  DelTensor(&ax);
}

/*Test from:
 *
 * PCA Nipals algorithm for Multi-Way Principal Component Analysis
 * Journal of Chemometrics, Vol 1, 41-56 (1987)
 * Svante Woold, Paul Geladi and Kim Esbensen
 */
void test3()
{
  tensor *ax;
  tensor *ay;

  tensor *axp;

  UPLSMODEL *m;


  tensor *py;
  matrix *pscores;

  puts(">>>>>>> Test 3: Compute Multi Way PLS and Prediction");

  NewTensor(&ax, 2);

  NewTensorMatrix(ax, 0, 3, 2);
  NewTensorMatrix(ax, 1, 3, 2);

  setTensorValue(ax, 0, 0, 0, 0.424264);  setTensorValue(ax, 0, 0, 1, 0.565685);
  setTensorValue(ax, 0, 1, 0, 0.565685);  setTensorValue(ax, 0, 1, 1, 0.424264);
  setTensorValue(ax, 0, 2, 0, 0.707101);  setTensorValue(ax, 0, 2, 1, 0.707101);

  setTensorValue(ax, 1, 0, 0, 0.565685);  setTensorValue(ax, 1, 0, 1, 0.424264);
  setTensorValue(ax, 1, 1, 0, 0.424264);  setTensorValue(ax, 1, 1, 1, 0.565685);
  setTensorValue(ax, 1, 2, 0, 0.707101);  setTensorValue(ax, 1, 2, 1, 0.707101);


  NewTensor(&ay, 2);

  NewTensorMatrix(ay, 0, 3, 1);
  NewTensorMatrix(ay, 1, 3, 1);

  setTensorValue(ay, 0, 0, 0, 1.0);
  setTensorValue(ay, 0, 1, 0, 2.0);
  setTensorValue(ay, 0, 2, 0, 3.0);


  setTensorValue(ay, 1, 0, 0, 1.0);
  setTensorValue(ay, 1, 1, 0, 1.5);
  setTensorValue(ay, 1, 2, 0, 2.0);


  NewTensor(&axp, 2);
  NewTensorMatrix(axp, 0, 1, 2);
  NewTensorMatrix(axp, 1, 1, 2);

  setTensorValue(axp, 0, 0, 0, 0.5);  setTensorValue(axp, 0, 0, 1, 0.6);
  setTensorValue(axp, 1, 0, 0, 0.6);  setTensorValue(axp, 1, 0, 1, 0.4);


  NewUPLSModel(&m);


  UPLS(ax, ay, 3, 1, 0, m, NULL);


  initTensor(&py);
  initMatrix(&pscores);

  UPLSScorePredictor(axp, m, 3, pscores);

  UPLSYPredictor(pscores, m, 3, py);

  PrintUPLSModel(m);

  puts("Prediction Scores");
  PrintMatrix(pscores);

  puts("Prediction Y");
  PrintTensor(py);

  DelTensor(&py);
  DelMatrix(&pscores);
  DelTensor(&axp);

  DelUPLSModel(&m);
  DelTensor(&ay);
  DelTensor(&ax);
}



/*
 * Test with Random numbers
 */
void test2()
{
  size_t i, j, k;
  tensor *ax;
  tensor *ay;

  UPLSMODEL *m;

  puts(">>>>>>> Test 2: Compute Multi Way PLS with random variables");

  NewTensor(&ax, 3);

  for(i = 0; i < ax->order; i++)
    NewTensorMatrix(ax, i, 30, 10);

  srand(time(0));
  for(i = 0; i < ax->order; i++){
    for(j = 0; j < ax->m[i]->row; j++){
      for(k = 0; k < ax->m[i]->col; k++)
          setTensorValue(ax, i, j, k, rand()/((double)(RAND_MAX)+1));
    }
  }



  NewTensor(&ay, 1);
  NewTensorMatrix(ay, 0, 30, 1);

  srand(time(0));
  for(i = 0; i < ay->order; i++){
    for(j = 0; j < ay->m[i]->row; j++){
      for(k = 0; k < ay->m[i]->col; k++)
          setTensorValue(ay, i, j, k, rand()/((double)(RAND_MAX)+1));
    }
  }

  NewUPLSModel(&m);

  UPLS(ax, ay, 3, 1, 0, m, NULL);

  PrintUPLSModel(m);

  DelUPLSModel(&m);
  DelTensor(&ay);
  DelTensor(&ax);
}




/*Test from:
 *
 * PCA Nipals algorithm for Multi-Way Principal Component Analysis
 * Journal of Chemometrics, Vol 1, 41-56 (1987)
 * Svante Woold, Paul Geladi and Kim Esbensen
 */
void test1()
{
  tensor *ax;
  tensor *ay;

  UPLSMODEL *m;

  puts(">>>>>>> Test 1: Compute Multi Way PLS");

  NewTensor(&ax, 2);

  NewTensorMatrix(ax, 0, 3, 2);
  NewTensorMatrix(ax, 1, 3, 2);

  setTensorValue(ax, 0, 0, 0, 0.424264);  setTensorValue(ax, 0, 0, 1, 0.565685);
  setTensorValue(ax, 0, 1, 0, 0.565685);  setTensorValue(ax, 0, 1, 1, 0.424264);
  setTensorValue(ax, 0, 2, 0, 0.707101);  setTensorValue(ax, 0, 2, 1, 0.707101);

  setTensorValue(ax, 1, 0, 0, 0.565685);  setTensorValue(ax, 1, 0, 1, 0.424264);
  setTensorValue(ax, 1, 1, 0, 0.424264);  setTensorValue(ax, 1, 1, 1, 0.565685);
  setTensorValue(ax, 1, 2, 0, 0.707101);  setTensorValue(ax, 1, 2, 1, 0.707101);


  NewTensor(&ay, 2);

  NewTensorMatrix(ay, 0, 3, 1);
  NewTensorMatrix(ay, 1, 3, 1);

  setTensorValue(ay, 0, 0, 0, 1.0);
  setTensorValue(ay, 0, 1, 0, 2.0);
  setTensorValue(ay, 0, 2, 0, 3.0);


  setTensorValue(ay, 1, 0, 0, 1.0);
  setTensorValue(ay, 1, 1, 0, 1.5);
  setTensorValue(ay, 1, 2, 0, 2.0);


  NewUPLSModel(&m);

  UPLS(ax, ay, 3, 1, 1, m, NULL);

  PrintUPLSModel(m);

  DelUPLSModel(&m);
  DelTensor(&ay);
  DelTensor(&ax);
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
  /*test8(); CRASH! */
    test9();
  return 0;
}

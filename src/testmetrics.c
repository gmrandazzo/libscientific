/* testmetrics.c
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
#include "metricspace.h"
#include "numeric.h"
#include "matrix.h"

void test4()
{
  puts("Test4: cubic spline interpolation.");
  matrix *xy;
  NewMatrix(&xy, 4, 2);
  xy->data[0][0] = 0.1; xy->data[0][1] = -0.62049958;
  xy->data[1][0] = 0.2; xy->data[1][1] = -0.28398668;
  xy->data[2][0] = 0.3; xy->data[2][1] = -0.00660095;
  xy->data[3][0] = 0.4; xy->data[3][1] = -0.24842440;

  matrix *S;
  initMatrix(&S);
  cubic_spline_interpolation(xy, &S);


  PrintMatrix(S);
  if(FLOAT_EQ(S->data[0][1], -0.620, 1e-3) == 1 &&
      FLOAT_EQ(S->data[0][2], 3.177, 1e-3) == 1 &&
      FLOAT_EQ(S->data[0][3], 0.000, 1e-3) == 1 &&
      FLOAT_EQ(S->data[0][4], 18.847, 1e-3) == 1 &&

      FLOAT_EQ(S->data[1][1], -0.284, 1e-3) == 1 &&
      FLOAT_EQ(S->data[1][2], 3.742, 1e-3) == 1 &&
      FLOAT_EQ(S->data[1][3], 5.654, 1e-3) == 1 &&
      FLOAT_EQ(S->data[1][4], -153.361, 1e-3) == 1 &&

      FLOAT_EQ(S->data[2][1], -0.007, 1e-3) == 1 &&
      FLOAT_EQ(S->data[2][2], 0.272, 1e-3) == 1 &&
      FLOAT_EQ(S->data[2][3], -40.354, 1e-3) == 1 &&
      FLOAT_EQ(S->data[2][4], 134.514, 1e-3) == 1){
    printf("OK!\n");
  }
  else{
    printf("ERROR\n");
  }

  DelMatrix(&S);
  DelMatrix(&xy);
}

void test3()
{
  puts("Test3: ROC and Precision-Recall test.");
  /* ROC curve test */
  dvector *y_score, *y_true;
  matrix *roc, *pr;
  double auc, ap;
  NewDVector(&y_score, 13);
  NewDVector(&y_true, 13);

  y_score->data[0] = 0.1;
  y_score->data[1] = 0.2;
  y_score->data[2] = 0.3;
  y_score->data[3] = 0.4;
  y_score->data[4] = 0.5;
  y_score->data[5] = 0.6;
  y_score->data[6] = 0.7;
  y_score->data[7] = 0.8;
  y_score->data[8] = 0.9;
  y_score->data[9] = 1.0;
  y_score->data[10] = 2.0;
  y_score->data[11] = 3.0;
  y_score->data[12] = 3.4;


  y_true->data[0] = 0;
  y_true->data[1] = 0;
  y_true->data[2] = 0;
  y_true->data[3] = 0;
  y_true->data[4] = 0;
  y_true->data[5] = 1;
  y_true->data[6] = 1;
  y_true->data[7] = 1;
  y_true->data[8] = 1;
  y_true->data[9] = 0;
  y_true->data[10] = 1;
  y_true->data[11] = 1;
  y_true->data[12] = 1;

  initMatrix(&roc);
  ROC(y_true, y_score,  &roc, &auc);
  PrintMatrix(roc);
  if(FLOAT_EQ(auc, 0.904762, 1e-6)){
    printf("AUC OK!\n");
  }else{
    printf("AUC ERROR!\n");
  }
  printf("AUC: %f\n", auc);

  initMatrix(&pr);
  PrecisionRecall(y_true, y_score,  &pr, &ap);
  PrintMatrix(pr);

  if(FLOAT_EQ(ap, 0.900425, 1e-6)){
    printf("AVERAGE PRECISION-RECALL OK!\n");
  }else{
    printf("AVERAGE PRECISION-RECALL ERROR!\n");
  }
  printf("AVERAGE PRECISION-RECALL: %f\n", ap);

  DelMatrix(&pr);
  DelMatrix(&roc);
  DelDVector(&y_true);
  DelDVector(&y_score);
}

void test2()
{
  puts("Test2: Natural Cubic spline interpolation and area under the curve tests.");
  /* Natural Cubic spline interpolation and area under the curve tests */
  matrix *xy, *interp_xy;

  NewMatrix(&xy, 6, 2);

  xy->data[0][0] = 1;
  xy->data[1][0] = 2;
  xy->data[2][0] = 3;
  xy->data[3][0] = 4;
  xy->data[4][0] = 5;
  xy->data[5][0] = 6;


  xy->data[0][1] = 1;
  xy->data[1][1] = 0.5;
  xy->data[2][1] = 0.33;
  xy->data[3][1] = 0.225;
  xy->data[4][1] = 0.2;
  xy->data[5][1] = 0.19;

  initMatrix(&interp_xy);
  interpolate(xy, 10, &interp_xy);
  double area = curve_area(xy, 10);
  if(FLOAT_EQ(area, 1.822294, 1e-6)){
    printf("Area calculation OK!\n");
  }else{
    printf("Area calculation ERROR!\n");
  }
  printf("area under the curve: %f\n", area);

  puts("X and Y");
  PrintMatrix(xy);
  puts("Interpolated X and Y");
  PrintMatrix(interp_xy);
  DelMatrix(&interp_xy);
  DelMatrix(&xy);
}

void test1()
{
  printf("Test1: Euclidean distance between two matrix.");
  /* Euclidean distance between two matrix test */
  matrix *m;
  matrix *dist;

  NewMatrix(&m, 4, 3);
  m->data[0][0] = 1; m->data[0][1] = 2; m->data[0][2] = 3;
  m->data[1][0] = 4; m->data[1][1] = 5; m->data[1][2] = 6;
  m->data[2][0] = 7; m->data[2][1] = 8; m->data[2][2] = 11;
  m->data[3][0] = 9; m->data[3][1] = 10; m->data[3][2] = 12;

  initMatrix(&dist);
  EuclideanDistance(m, m, &dist);

  puts("Matrix");
  PrintMatrix(m);
  puts("Distance Matrix");
  PrintMatrix(dist);

  DelMatrix(&dist);
  DelMatrix(&m);
}

int main(void)
{
  test1();
  test2();
  test3();
  test4();
}

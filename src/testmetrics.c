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

void test3()
{
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
  cubic_spline_interpolation(xy, 10, &interp_xy);
  double area = curve_area(xy, 10);
  if(FLOAT_EQ(area, 1.816975, 1e-6)){
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


void test2()
{
  matrix *xy, *interp_xy;

   NewMatrix(&xy, 4, 2);
   xy->data[0][0] = -1;
   xy->data[1][0] = -0.5;
   xy->data[2][0] = 0.0;
   xy->data[3][0] = 0.5;

   xy->data[0][1] = 0.86199480;
   xy->data[1][1] = 0.95802009;
   xy->data[2][1] = 1.0986123;
   xy->data[3][1] = 1.2943767;

   /* results
   0  0.861995  0.175638  0.0 0.0656509
   1  0.95802 0.224876  0.0984763 0.028281
   2  1.09861 0.344563  0.140898  -0.093932
   */
  initMatrix(&interp_xy);
  cubic_spline_interpolation(xy, 10, &interp_xy);
  printf("area: %f\n", curve_area(xy, 10));
  /* result

  0 0 75 −0.659292 0.219764
  3 225 76.9779 1.31858 −0.153761
  5 383 80.4071 0.396018 −0.177237
  8 623 77.9978 −1.19912 0.0799115

  */
  puts("X and Y");
  PrintMatrix(xy);
  puts("Interpolated X and Y");
  PrintMatrix(interp_xy);
  DelMatrix(&interp_xy);
  DelMatrix(&xy);
}


void test1()
{
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
}

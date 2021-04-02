/* testinterpolate.c
*
* Copyright (C) <2018>  Giuseppe Marco Randazzo
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
#include <math.h>
#include "interpolate.h"
#include "matrix.h"
#include "numeric.h"


void test2()
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


void test1()
{
  puts("Test1: Natural Cubic spline interpolation and area under the curve tests.");
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

int main(void)
{
  test1();
  test2();
}

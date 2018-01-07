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
#include "matrix.h"

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
}
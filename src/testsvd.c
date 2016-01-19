/* testsvd.c
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
#include "matrix.h"
#include "vector.h"
#include "svd.h"

void test1()
{
  matrix *mx, *u, *v;
  dvector *w;
  NewMatrix(&mx, 4, 5);
  MatrixSet(mx, 0.f);
  mx->data[0][0] = 1;  mx->data[0][4] = 2;
  mx->data[1][2] = 3;
  mx->data[3][1] = 4;
  
  PrintMatrix(mx);
  
  initMatrix(&u);
  initMatrix(&v);
  initDVector(&w);
  svd(mx, &u, &w, &v);
  
  puts("U");
  PrintMatrix(u);
  puts("W");
  PrintDVector(w);
  puts("V");
  PrintMatrix(v);
  
  DelMatrix(&u);
  DelMatrix(&v);
  DelDVector(&w);
  DelMatrix(&mx);
}

int main(void)
{
  test1();
  return 0;
}
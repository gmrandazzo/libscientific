/* testica.c
*
* Copyright (C) <2022>  Giuseppe Marco Randazzo
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

#include <stdlib.h>
#include <stdlib.h>
#include "numeric.h"
#include "ica.h"
#include "datasets.h"
#include "numeric.h"
#include "scientificinfo.h"

void test2()
{
  puts("Test ICA 2: ICA on iris dataset");
  matrix *m, *_;
  ICAMODEL *model;

  initMatrix(&m);
  initMatrix(&_);
  iris(m, _);

  NewICAModel(&model);

  ICA(m, 1, 2, model);

  PrintICA(model);
  DelICAModel(&model);
  DelMatrix(&m);
  DelMatrix(&_);
}


void test1()
{
  puts("Test ICA 1: ICA on iris dataset");
  matrix *m, *_;
  ICAMODEL *model;
  size_t i;
  size_t j;

  NewMatrix(&m, 10, 3);
  m->data[0][0] = -1.; m->data[0][1] = -1.; m->data[0][2] = -2.;
  m->data[1][0] = 2.75643348; m->data[1][1] = 3.26710563; m->data[1][2] = 4.02353911;
  m->data[2][0] = -0.84668509; m->data[2][1] = -1.64556477; m->data[2][2] = -0.49224986;
  m->data[3][0] = 0.52000394; m->data[3][1] = 1.92666864; m->data[3][2] = 0.44667258;
  m->data[4][0] = -0.15235893; m->data[4][1] = -1.52062391;  m->data[4][2] = 0.32701716;
  m->data[5][0] = 1.39949457; m->data[5][1] = 2.14419173; m->data[5][2] = 1.5436863;
  m->data[6][0] = -2.27972909; m->data[6][1] = -2.80653121; m->data[6][2] = -3.0862603;
  m->data[7][0] = -1.67717986; m->data[7][1] = -2.61636771; m->data[7][2] = -2.29354756;
  m->data[8][0] = 1.2186072;  m->data[8][1] = 1.72041471; m->data[8][2] = 0.93902191;
  m->data[9][0] = -2.28790332; m->data[9][1] = -3.14395166; m->data[9][2] = -3.43185497;


  NewICAModel(&model);

  ICA(m, 0, 2, model);

  PrintICA(model);
  DelICAModel(&model);
  DelMatrix(&m);
  DelMatrix(&_);
}


int main(void)
{
  test1();
  return 0;
}

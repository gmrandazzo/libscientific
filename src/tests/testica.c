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
  for(i = 0; i < 10; i++){
    for(j = 0; j < 3; j++){
      m->data[i][j] = randDouble(-2, 2);
    }
  }

  NewICAModel(&model);

  ICA(m, 1, 2, model);

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

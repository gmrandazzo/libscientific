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
#include "scientificinfo.h"

void test1()
{
  puts("Test ICA 1: ICA on iris dataset");
  matrix *m, *_;
  PCAMODEL *model;

  initMatrix(&m);
  initMatrix(&_);
  iris(m, _);

  NewICAModel(&model);

  ICA(m, 1, 2, model, NULL);
  PrintICA(model);
  DelICAModel(&model);
  DelMatrix(&m);
  DelMatrix(&_);
}


int main(void)
{
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
  return 0;
}

/* Unit tests for the list module.
 * Copyright (C) 2017-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "vector.h"
#include "list.h"

/* List tests using:
 * void initDVectorList(dvectorlist **lst);
 * void NewDVectorList(dvectorlist **lst, size_t size_);
 * void DelDVectorList(dvectorlist **lst);
 * void DVectorListAppend(dvectorlist *lst, dvector *d);
 */
void test1()
{
  int i, j;
  dvector *v;
  dvectorlist *lst;

  printf("Test 1\n");
  NewDVector(&v, 100);

  for(i = 0; i < 100; i++){
    v->data[i] = i;
  }

  printf("Appending dvector to dvectorlist\n");
  initDVectorList(&lst);
  DVectorListAppend(lst, v);

  printf("Final output\n");
  for(i = 0; i < lst->size; i++){
    for(j = 0; j < lst->d[i]->size; j++)
      printf("%f\n", lst->d[i]->data[j]);
    printf("--------\n");
  }

  DelDVector(&v);
  DelDVectorList(&lst);
}

int main(void)
{
  test1();
  return 0;
}

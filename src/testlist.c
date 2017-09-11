/* testlist.c
*
* Copyright (C) <2017>  Giuseppe Marco Randazzo
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

#include "vector.h"
#include "list.h"


/* Allocate the vector by using the NewDVector function */
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
  DVectorListAppend(&lst, v);

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

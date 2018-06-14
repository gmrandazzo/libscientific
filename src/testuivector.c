/* testuivector.c
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

#include <stdlib.h>
#include "vector.h"


/* Extend two vector into one vector*/
void test3()
{
   int i;
  uivector *v1;
  uivector *v2;
  uivector *v1v2;

  printf("Test 3\n");

  initUIVector(&v1);

  printf("Creating v1\n");
  for(i = 1; i < 100; i++){
    UIVectorAppend(&v1, i);
  }


  printf("Creating v2\n");
  NewUIVector(&v2, 100);

  for(i = 0; i < 100; i++){
    v2->data[i] = i;
  }

  printf("Appending v1 to v2\n");
  v1v2 = UIVectorExtend(v1,v2);

  printf("Final output\n");
  for(i = 0; i < v1v2->size; i++){
   printf("%u\n", (unsigned int)v1v2->data[i]);
  }

  UIVectorRemoveAt(&v1, 48);

  DelUIVector(&v1v2);
  DelUIVector(&v2);
  DelUIVector(&v1);
}

/* Initialize the uivector by using initUIVector function and then append values with UIVectorAppend function*/
void test2()
{
  int i;
  uivector *v;

  printf("Test 2\n");

  initUIVector(&v);

  printf("Appending 100 value\n");
  for(i = 0; i < 100; i++){
    UIVectorAppend(&v, i);
  }

  printf("Final output\n");
  for(i = 0; i < v->size; i++){
   printf("%u\n", (unsigned int)v->data[i]);
  }

  DelUIVector(&v);

}



/* Allocate the vector by using the NewDVector function */
void test1()
{
  int i;
  uivector *v;

  printf("Test 1\n");

  NewUIVector(&v, 100);

  for(i = 0; i < 100; i++){
    v->data[i] = i;
  }

  printf("Appending 123 to vector\n");
  UIVectorAppend(&v, 123);

  printf("Final output\n");
  for(i = 0; i < v->size; i++){
   printf("%u\n", (unsigned int)v->data[i]);
  }

  DelUIVector(&v);

}

int main(void)
{
  test1();
  test2();
  test3();
  return 0;
}

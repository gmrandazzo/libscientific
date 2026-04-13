/* Unit tests for the ivector module.
 * Copyright (C) 2016-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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

#include <stdlib.h>
#include "vector.h"


/* Extend two vector into one vector*/
void test3()
{
   int i;
  ivector *v1;
  ivector *v2;
  ivector *v1v2;
  
  printf("Test 3\n");
  
  initIVector(&v1);
  
  printf("Creating v1\n");
  for(i = 1; i < 100; i++){
    IVectorAppend(v1, i);
  }
  

    
  printf("Creating v2\n");
  NewIVector(&v2, 100);

  for(i = 0; i < 100; i++){
    v2->data[i] = i;
  }
 
  printf("Appending v1 to v2\n");
  v1v2 = IVectorExtend(v1,v2);
  
  printf("Final output\n");
  for(i = 0; i < v1v2->size; i++){
   printf("%d\n", v1v2->data[i]);
  }
  
  DelIVector(&v1v2);
  DelIVector(&v2);
  DelIVector(&v1);
}

/* Initialize the ivector by using initIVector function and then append values with IVectorAppend function*/
void test2()
{
  int i;
  ivector *v;
  
  printf("Test 2\n");
  
  initIVector(&v);
  
  printf("Appending 100 value\n");
  for(i = 0; i < 100; i++){
    IVectorAppend(v, i);
  }
  
  printf("Final output\n");
  for(i = 0; i < v->size; i++){
   printf("%d\n", v->data[i]);
  }
  
  DelIVector(&v);

}



/* Allocate the vector by using the NewDVector function */
void test1()
{
  int i;
  ivector *v;
  
  printf("Test 1\n");
  
  NewIVector(&v, 100);
  
  for(i = 0; i < 100; i++){
    v->data[i] = i;
  }
  
  printf("Appending 123 to vector\n");
  IVectorAppend(v, 123);
  
  printf("Final output\n");
  for(i = 0; i < v->size; i++){
   printf("%d\n", v->data[i]);
  }
  
  DelIVector(&v);

}

int main(void)
{  
  test1();
  test2();
  test3();
  return 0;
}

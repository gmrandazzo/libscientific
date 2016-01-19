/* testarray.c
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
#include <stdlib.h>
#include <time.h>
#include "array.h"
#include "vector.h"

void test6()
{
  array *t1;
  
  matrix *p;
  
  dvector *v;
  
  puts("Test 6 - the result of this vector is 455.000");
  
  NewArray(&t1, 2);
  NewArrayMatrix(&t1, 0, 1, 3);  /* i = 1 j = 3 */
  NewArrayMatrix(&t1, 1, 1, 3);
  
  setArrayValue(t1, 0, 0, 0, 1); setArrayValue(t1, 0, 0, 1, 2);   setArrayValue(t1, 0, 0, 2, 3);
  setArrayValue(t1, 1, 0, 0, 4); setArrayValue(t1, 1, 0, 1, 5);   setArrayValue(t1, 1, 0, 2, 6);
  
  NewMatrix(&p, 3, 2);
  
  setMatrixValue(p, 0, 0, 5.0); setMatrixValue(p, 0, 1, 20);
  setMatrixValue(p, 1, 0, 10);  setMatrixValue(p, 1, 1, 25);
  setMatrixValue(p, 2, 0, 15);  setMatrixValue(p, 2, 1, 30);
  
  
  NewDVector(&v, 1);
  
  ArrayMatrixDotProduct(t1, p, v);
  
  PrintDVector(v);
  
  DelDVector(&v);
  DelMatrix(&p);
  DelArray(&t1);
}

void test5()
{
  array *t1, *t2;
  
  matrix *p;
  
  dvector *v;
  
  puts("Test 5 - Computing Array Dvector dot Product");
  
  NewArray(&t1, 2); /*k = 2 */
  NewArrayMatrix(&t1, 0, 1, 3);  /* i = 1 j = 3 */
  NewArrayMatrix(&t1, 1, 1, 3);
  
  NewArray(&t2, 2);
  NewArrayMatrix(&t2, 0, 3, 1);
  NewArrayMatrix(&t2, 1, 3, 1);
  
  
  setArrayValue(t1, 0, 0, 0, 1); setArrayValue(t1, 0, 0, 1, 2);   setArrayValue(t1, 0, 0, 2, 3);
  setArrayValue(t1, 1, 0, 0, 4); setArrayValue(t1, 1, 0, 1, 5);   setArrayValue(t1, 1, 0, 2, 6);
  
  
  NewDVector(&v, 1);
  
  setDVectorValue(v, 0, 5);
  
  
  PrintArray(t1);
  
  ArrayTranspose(t1, t2);
  
  PrintArray(t2);
 
 
  puts("TransposedArrayDVectorProduct p[k][j] where k is the order and j is the column ");
  NewMatrix(&p, 2, 3); /* k = 2  j = 3    p[k][j]*/
  TransposedArrayDVectorProduct(t2, v, p);
  PrintMatrix(p);
  DelMatrix(&p);
  
  NewMatrix(&p, 3, 2); /* p[j][k] */
  puts("DvectorArrayDotProduct p[j][k] where k is the order and j is the column ");
  DvectorArrayDotProduct(t1, v, p);
  PrintMatrix(p);
  DelMatrix(&p);
  
  DelDVector(&v);
  DelArray(&t1);
  DelArray(&t2);
  
}

/* This test show how to use NewArray() with MeanCenteredArray() and ArrayAutoScaling() */
void test4()
{
  size_t i, j, k;
  array *t, *tc;

  puts(">>>>>>>> Test 4: Create a Array of order 4 with 5 objects, meancentering and autoscaling this by using NewArray()");
  
  srand(time(0));
  NewArray(&t, 4);
  NewArray(&tc, 4);
      
  for(i = 0; i < t->order; i ++){
    NewArrayMatrix(&t, i, 5, i+1);
    
    NewArrayMatrix(&tc, i, 5, i+1);
    
    for(j = 0; j < t->m[i]->col; j++){
      for(k = 0; k < t->m[i]->row; k++){
        setArrayValue(t, i, k, j, rand() % 300);
      }
    }
  }
  
  MeanCenteredArray(t, tc);
  puts("Original Array");
  PrintArray(t);
  puts("Mean Centered Array");
  PrintArray(tc);
  
  DelArray(&t);
  DelArray(&tc);
  
}

/* This test show how to use NewArray() with MeanCenteredArray()*/
void test3()
{
  size_t i, j, k;
  array *t, *tc;

  puts(">>>>>>>> Test 3: Create a Array of order 4 with 5 objects and meancentering this by using NewArray()");
  
  srand(time(0));
  NewArray(&t, 4);
  NewArray(&tc, 4);
      
  for(i = 0; i < t->order; i ++){
    NewArrayMatrix(&t, i, 5, i+1);
    NewArrayMatrix(&tc, i, 5, i+1);
    for(j = 0; j < t->m[i]->col; j++){
      for(k = 0; k < t->m[i]->row; k++){
        setArrayValue(t, i, k, j, rand() % 300);
      }
    }
  }
  
  MeanCenteredArray(t, tc);
  puts("Original Array");
  PrintArray(t);
  puts("Mean Centered Array");
  PrintArray(tc);
  
  DelArray(&t);
  DelArray(&tc);
}

/* This test show how to use NewArray with setArrayValue()*/
void test2()
{
  size_t i, j, k;
  array *t;

  puts(">>>>>>>> Test 2: Create a Array of order 4 with 20 objects by using NewArray()");
  
  srand(time(0));
  NewArray(&t, 4);
      
  for(i = 0; i < t->order; i ++){
    NewArrayMatrix(&t, i, 20, i+1);
    for(j = 0; j < t->m[i]->col; j++){
      for(k = 0; k < t->m[i]->row; k++){
        setArrayValue(t, i, k, j, rand() % 300);
      }
    }
  }
  
  PrintArray(t);
  DelArray(&t);
}


void test1_2()
{
  size_t i, j, k;
  array *t;
  matrix *m;

  puts(">>>>>>>> Test 1: Create a Array of order 4 with 10 objects by using initArray()");
  
  srand(time(0));
  initArray(&t);
  NewMatrix(&m, 10, 5);
  
  for(i = 0; i < 4; i ++){
    for(j = 0; j < 10; j++){
      for(k = 0; k < 5; k++){
        setMatrixValue(m, j , k, rand() % 300);
      }
    }
    /*PrintMatrix(m);*/
    ArrayAppendMatrix(&t, m); /*only with initArray*/
  }
  DelMatrix(&m);
  
  PrintArray(t);
  DelArray(&t);
  
}

void test1_1()
{
  size_t i, j, k;
  array *t;
  matrix *m;

  puts(">>>>>>>> Test 1: Create a Array of order 4 with 10 objects by using NewArray()");
  
  srand(time(0));
  NewArray(&t, 4);
  
  for(i = 0; i < 4; i ++){
    
    NewMatrix(&m, 10, i+1);
    for(j = 0; j < 10; j++){
      for(k = 0; k < i+1; k++){
        setMatrixValue(m, j , k, rand() % 300);
      }
    }
    
    /*PrintMatrix(m);*/
    ArrayAppendMatrixAt(&t, i, m); /*only with initArray*/
    DelMatrix(&m);
  }
  
  PrintArray(t);
  DelArray(&t);
  
}

/* This test show how to use initArray with ArrayAppendMatrix */
void test1()
{
  size_t i, j, k;
  array *t;
  matrix *m;

  puts(">>>>>>>> Test 1: Create a Array of order 4 with 10 objects by using initArray()");
  
  srand(time(0));
  initArray(&t);
  
  for(i = 0; i < 4; i ++){
    
    NewMatrix(&m, 10, i+1);
    for(j = 0; j < 10; j++){
      for(k = 0; k < i+1; k++){
        setMatrixValue(m, j , k, rand() % 300);
      }
    }
    
    /*PrintMatrix(m);*/
    ArrayAppendMatrix(&t, m); /*only with initArray*/
    DelMatrix(&m);
  }
  
  PrintArray(t);
  DelArray(&t);
  
}

int main(void)
{
  test1();
  test1_2();
  test2();
  test3();
  test4();
  test5();
  test6();
  return 0;
}

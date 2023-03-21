/* testtensor.c
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
#include "tensor.h"
#include "vector.h"

void test6()
{
  tensor *t1;

  matrix *p;

  dvector *v;

  puts("Test 6 - the result of this vector is 455.000");

  NewTensor(&t1, 2);
  NewTensorMatrix(t1, 0, 1, 3);  /* i = 1 j = 3 */
  NewTensorMatrix(t1, 1, 1, 3);

  setTensorValue(t1, 0, 0, 0, 1); setTensorValue(t1, 0, 0, 1, 2);   setTensorValue(t1, 0, 0, 2, 3);
  setTensorValue(t1, 1, 0, 0, 4); setTensorValue(t1, 1, 0, 1, 5);   setTensorValue(t1, 1, 0, 2, 6);

  NewMatrix(&p, 3, 2);

  setMatrixValue(p, 0, 0, 5.0); setMatrixValue(p, 0, 1, 20);
  setMatrixValue(p, 1, 0, 10);  setMatrixValue(p, 1, 1, 25);
  setMatrixValue(p, 2, 0, 15);  setMatrixValue(p, 2, 1, 30);


  NewDVector(&v, 1);

  TensorMatrixDotProduct(t1, p, v);

  PrintDVector(v);

  DelDVector(&v);
  DelMatrix(&p);
  DelTensor(&t1);
}

void test5()
{
  tensor *t1, *t2;

  matrix *p;

  dvector *v;

  puts("Test 5 - Computing Tensor Dvector dot Product");

  NewTensor(&t1, 2); /*k = 2 */
  NewTensorMatrix(t1, 0, 1, 3);  /* i = 1 j = 3 */
  NewTensorMatrix(t1, 1, 1, 3);

  NewTensor(&t2, 2);
  NewTensorMatrix(t2, 0, 3, 1);
  NewTensorMatrix(t2, 1, 3, 1);


  setTensorValue(t1, 0, 0, 0, 1); setTensorValue(t1, 0, 0, 1, 2);   setTensorValue(t1, 0, 0, 2, 3);
  setTensorValue(t1, 1, 0, 0, 4); setTensorValue(t1, 1, 0, 1, 5);   setTensorValue(t1, 1, 0, 2, 6);


  NewDVector(&v, 1);

  setDVectorValue(v, 0, 5);


  PrintTensor(t1);

  TensorTranspose(t1, t2);

  PrintTensor(t2);


  puts("TransposedTensorDVectorProduct p[k][j] where k is the order and j is the column ");
  NewMatrix(&p, 2, 3); /* k = 2  j = 3    p[k][j]*/
  TransposedTensorDVectorProduct(t2, v, p);
  PrintMatrix(p);
  DelMatrix(&p);

  NewMatrix(&p, 3, 2); /* p[j][k] */
  puts("DvectorTensorDotProduct p[j][k] where k is the order and j is the column ");
  DvectorTensorDotProduct(t1, v, p);
  PrintMatrix(p);
  DelMatrix(&p);

  DelDVector(&v);
  DelTensor(&t1);
  DelTensor(&t2);

}

/* This test show how to use NewTensor() with MeanCenteredTensor() and TensorAutoScaling() */
void test4()
{
  size_t i, j, k;
  tensor *t, *tc;

  puts(">>>>>>>> Test 4: Create a Tensor of order 4 with 5 objects, meancentering and autoscaling this by using NewTensor()");

  srand(time(0));
  NewTensor(&t, 4);
  NewTensor(&tc, 4);

  for(i = 0; i < t->order; i ++){
    NewTensorMatrix(t, i, 5, i+1);

    NewTensorMatrix(tc, i, 5, i+1);

    for(j = 0; j < t->m[i]->col; j++){
      for(k = 0; k < t->m[i]->row; k++){
        setTensorValue(t, i, k, j, rand() % 300);
      }
    }
  }

  MeanCenteredTensor(t, tc);
  puts("Original Tensor");
  PrintTensor(t);
  puts("Mean Centered Tensor");
  PrintTensor(tc);

  DelTensor(&t);
  DelTensor(&tc);

}

/* This test show how to use NewTensor() with MeanCenteredTensor()*/
void test3()
{
  size_t i, j, k;
  tensor *t, *tc;

  puts(">>>>>>>> Test 3: Create a Tensor of order 4 with 5 objects and meancentering this by using NewTensor()");

  srand(time(0));
  NewTensor(&t, 4);
  NewTensor(&tc, 4);

  for(i = 0; i < t->order; i ++){
    NewTensorMatrix(t, i, 5, i+1);
    NewTensorMatrix(tc, i, 5, i+1);
    for(j = 0; j < t->m[i]->col; j++){
      for(k = 0; k < t->m[i]->row; k++){
        setTensorValue(t, i, k, j, rand() % 300);
      }
    }
  }

  MeanCenteredTensor(t, tc);
  puts("Original Tensor");
  PrintTensor(t);
  puts("Mean Centered Tensor");
  PrintTensor(tc);

  DelTensor(&t);
  DelTensor(&tc);
}

/* This test show how to use NewTensor with setTensorValue()*/
void test2()
{
  size_t i, j, k;
  tensor *t;

  puts(">>>>>>>> Test 2: Create a Tensor of order 4 with 20 objects by using NewTensor()");

  srand(time(0));
  NewTensor(&t, 4);

  for(i = 0; i < t->order; i ++){
    NewTensorMatrix(t, i, 20, i+1);
    for(j = 0; j < t->m[i]->col; j++){
      for(k = 0; k < t->m[i]->row; k++){
        setTensorValue(t, i, k, j, rand() % 300);
      }
    }
  }

  PrintTensor(t);
  DelTensor(&t);
}


void test1_2()
{
  size_t i, j, k;
  tensor *t;
  matrix *m;

  puts(">>>>>>>> Test 1: Create a Tensor of order 4 with 10 objects by using initTensor()");

  srand(time(0));
  initTensor(&t);
  NewMatrix(&m, 10, 5);

  for(i = 0; i < 4; i ++){
    for(j = 0; j < 10; j++){
      for(k = 0; k < 5; k++){
        setMatrixValue(m, j , k, rand() % 300);
      }
    }
    /*PrintMatrix(m);*/
    TensorAppendMatrix(t, m); /*only with initTensor*/
  }
  DelMatrix(&m);

  PrintTensor(t);
  DelTensor(&t);

}

void test1_1()
{
  size_t i, j, k;
  tensor *t;
  matrix *m;

  puts(">>>>>>>> Test 1: Create a Tensor of order 4 with 10 objects by using NewTensor()");

  srand(time(0));
  NewTensor(&t, 4);

  for(i = 0; i < 4; i ++){

    NewMatrix(&m, 10, i+1);
    for(j = 0; j < 10; j++){
      for(k = 0; k < i+1; k++){
        setMatrixValue(m, j , k, rand() % 300);
      }
    }

    /*PrintMatrix(m);*/
    TensorAppendMatrixAt(t, i, m); /*only with initTensor*/
    DelMatrix(&m);
  }

  PrintTensor(t);
  DelTensor(&t);

}

/* This test show how to use initTensor with TensorAppendMatrix */
void test1()
{
  size_t i, j, k;
  tensor *t;
  matrix *m;

  puts(">>>>>>>> Test 1: Create a Tensor of order 4 with 10 objects by using initTensor()");

  srand(time(0));
  initTensor(&t);

  for(i = 0; i < 4; i ++){

    NewMatrix(&m, 10, i+1);
    for(j = 0; j < 10; j++){
      for(k = 0; k < i+1; k++){
        setMatrixValue(m, j , k, rand() % 300);
      }
    }

    /*PrintMatrix(m);*/
    TensorAppendMatrix(t, m); /*only with initTensor*/
    DelMatrix(&m);
  }

  PrintTensor(t);
  DelTensor(&t);

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

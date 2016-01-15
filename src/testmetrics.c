#include <stdio.h>
#include "metricspace.h"
#include "matrix.h"

void test1()
{
  matrix *m;
  matrix *dist;
  
  NewMatrix(&m, 4, 3);
  m->data[0][0] = 1; m->data[0][1] = 2; m->data[0][2] = 3;
  m->data[1][0] = 4; m->data[1][1] = 5; m->data[1][2] = 6;
  m->data[2][0] = 7; m->data[2][1] = 8; m->data[2][2] = 11;
  m->data[3][0] = 9; m->data[3][1] = 10; m->data[3][2] = 12;
  
  initMatrix(&dist);
  EuclideanDistance(m, m, &dist);

  puts("Matrix");
  PrintMatrix(m);
  puts("Distance Matrix");
  PrintMatrix(dist);
  
  DelMatrix(&dist);
  DelMatrix(&m);
}

int main(void)
{
  test1();
}
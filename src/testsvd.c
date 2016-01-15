#include <stdio.h>
#include "matrix.h"
#include "vector.h"
#include "svd.h"

void test1()
{
  matrix *mx, *u, *v;
  dvector *w;
  NewMatrix(&mx, 4, 5);
  MatrixSet(mx, 0.f);
  mx->data[0][0] = 1;  mx->data[0][4] = 2;
  mx->data[1][2] = 3;
  mx->data[3][1] = 4;
  
  PrintMatrix(mx);
  
  initMatrix(&u);
  initMatrix(&v);
  initDVector(&w);
  svd(mx, &u, &w, &v);
  
  puts("U");
  PrintMatrix(u);
  puts("W");
  PrintDVector(w);
  puts("V");
  PrintMatrix(v);
  
  DelMatrix(&u);
  DelMatrix(&v);
  DelDVector(&w);
  DelMatrix(&mx);
}

int main(void)
{
  test1();
  return 0;
}
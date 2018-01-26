#include <stdio.h>
#include "optimization.h"
#include "vector.h"

double square(double x){
  return x*x;
}

double rosen(dvector *x)
{
  size_t i;
  double res = 0.f;
  i = x->size-1;
  while(i > 0){
    res += (100 * square((x->data[i]-square(x->data[i-1])))) +square(1-x->data[i-1]);
    i--;
  }
  return res;
}


void test1(){
  dvector *x0, *best;
  NewDVector(&x0, 5);
  /*1.3, 0.7, 0.8, 1.9, 1.2*/
  x0->data[0] = 1.3;
  x0->data[1] = 0.7;
  x0->data[2] = 0.8;
  x0->data[3] = 1.9;
  x0->data[4] = 1.2;

  printf("Init ROSEN : %f\n", rosen(x0));

  initDVector(&best);
  double min = NelderMeadSimplex(&rosen, x0, NULL, 1e-10, 1000, &best);
  printf("min results obtained: %f\n",  min);
  puts("Best results obtained vector");
  PrintDVector(best);
  /*litterature*/
  best->data[0] = 0.995260;
  best->data[1] = 0.990385;
  printf("min results litterature: %f\n",rosen(best));
  puts("Best results litterature vector");
  PrintDVector(best);
  DelDVector(&x0);
  DelDVector(&best);
}

int main(void){
  test1();
  return 0;
}

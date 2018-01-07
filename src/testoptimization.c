#include <stdio.h>
#include "optimization.h"
#include "vector.h"

double rosen(dvector *x)
{
  return (100*(x->data[1]-x->data[0]*x->data[0])*(x->data[1]-x->data[0]*x->data[0])+(1.0-x->data[0])*(1.0-x->data[0]));
}


void test1(){
  dvector *x0, *best;
  NewDVector(&x0, 2);
  x0->data[0] = -1.2; x0->data[1] = 1.0;
  initDVector(&best);
  double min = NelderMeadSimplex(&rosen, x0, NULL, 1e-3, 1000, &best, minimization);
  printf("Best results: %.3f\n",  min);
  printf("rosen: %f\n", rosen(best));
  best->data[0] = 0.995260;
  best->data[1] = 0.990385;
  printf("rosen litterature: %f\n", rosen(best));
  PrintDVector(best);
  DelDVector(&x0);
  DelDVector(&best);
}

int main(void){
  test1();
  return 0;
}

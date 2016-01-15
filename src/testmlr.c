#include <stdio.h>
#include "mlr.h"
#include "numeric.h"

void test2()
{
  matrix *mx, *my;
  MLRMODEL *m;
  ssignal s = SIGSCIENTIFICRUN;
  puts("Test1: Simple Calculation MLR Model with LOOCV and Model Consistency");


  NewMatrix(&mx, 17, 1);
  NewMatrix(&my, 17, 1);

  srand(17);
  for(size_t i = 0; i < 17; i++){
    setMatrixValue(mx, i, 0, randDouble(0,20));
    setMatrixValue(my, i, 0, randDouble(0,1));
  }

  puts("X");
  PrintMatrix(mx);
  puts("Y");
  PrintMatrix(my);

  NewMLRModel(&m);
  MLR(mx, my, m, &s);

  MLRYScrambling(mx, my, 5, 0, 0, 0, &m->r2q2scrambling, &s);
  MLRLOOCV(mx, my, &m->q2y, &m->sdep, &m->bias, &m->predicted_y, &m->pred_residuals, &s);

  PrintMLR(m);

  DelMLRModel(&m);
  DelMatrix(&mx);
  DelMatrix(&my);
}


void test1()
{
  matrix *mx, *my;
  MLRMODEL *m;
  ssignal s = SIGSCIENTIFICRUN;
  puts("Test1: Simple Calculation MLR Model with LOOCV and Model Consistency");

  NewMatrix(&mx, 5, 1);
  NewMatrix(&my, 5, 1);

  setMatrixValue(mx, 0, 0, 1.8);
  setMatrixValue(mx, 1, 0, 3.2);
  setMatrixValue(mx, 2, 0, 5.6);
  setMatrixValue(mx, 3, 0, 9.43);
  setMatrixValue(mx, 4, 0, 13.7);

  setMatrixValue(my, 0, 0, 2.1);
  setMatrixValue(my, 1, 0, 1.8);
  setMatrixValue(my, 2, 0, 6.4);
  setMatrixValue(my, 3, 0, 8.2);
  setMatrixValue(my, 4, 0, 14);


  puts("X");
  PrintMatrix(mx);
  puts("Y");
  PrintMatrix(my);
  NewMLRModel(&m);
  MLR(mx, my, m, &s);

  MLRRandomGroupsCV(mx, my, 2, 20, &m->q2y, &m->sdep, &m->bias, &m->predicted_y, &m->pred_residuals, &s);
  /*MLRLOOCV(mx, my, &m->q2y, &m->sdep, &m->bias, &m->predicted_y, &m->pred_residuals, &s);*/
  MLRYScrambling(mx, my, 5, 0, 0, 0, &m->r2q2scrambling, &s);

  PrintMLR(m);

  DelMLRModel(&m);
  DelMatrix(&mx);
  DelMatrix(&my);
}

int main()
{
  test1();
  test2();
  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "array.h"
#include "upca.h"

void test4()
{
  array *t;
  UPCAMODEL *m;

  puts(">>>>>>> Test 4: Compute Multi Way PCA");

  NewArray(&t, 2);

  NewArrayMatrix(&t, 0, 4, 2);
  NewArrayMatrix(&t, 1, 4, 2);

  setArrayValue(t, 0, 0, 0, 7);  setArrayValue(t, 0, 0, 1, 9);
  setArrayValue(t, 0, 1, 0, 5);  setArrayValue(t, 0, 1, 1, 4);
  setArrayValue(t, 0, 2, 0, 2);  setArrayValue(t, 0, 2, 1, 10);
  setArrayValue(t, 0, 3, 0, 4);  setArrayValue(t, 0, 3, 1, 16);

  setArrayValue(t, 1, 0, 0, 7);  setArrayValue(t, 1, 0, 1, 9);
  setArrayValue(t, 1, 1, 0, 5);  setArrayValue(t, 1, 1, 1, 4);
  setArrayValue(t, 1, 2, 0, 4);  setArrayValue(t, 1, 2, 1, 10);
  setArrayValue(t, 1, 3, 0, 6);  setArrayValue(t, 1, 3, 1, 16);

  puts("Array");
  PrintArray(t);

  NewUPCAModel(&m);

  UPCA(t, 2, 1, m, NULL);

  DelUPCAModel(&m);
  DelArray(&t);
}

/*Test from:
 *
 * PCA Nipals algorithm for Multi-Way Principal Component Analysis
 * Journal of Chemometrics, Vol 1, 41-56 (1987)
 * Svante Woold, Paul Geladi and Kim Esbensen
 */
void test3()
{
  array *t;

  UPCAMODEL *m;

  array *tpred;
  matrix *predscore;

  puts(">>>>>>> Test 3: Compute Multi Way PCA and Prediction");

  NewArray(&t, 2);

  NewArrayMatrix(&t, 0, 3, 2);
  NewArrayMatrix(&t, 1, 3, 2);

  setArrayValue(t, 0, 0, 0, 0.424264);  setArrayValue(t, 0, 0, 1, 0.565685);
  setArrayValue(t, 0, 1, 0, 0.565685);  setArrayValue(t, 0, 1, 1, 0.424264);
  setArrayValue(t, 0, 2, 0, 0.707101);  setArrayValue(t, 0, 2, 1, 0.707101);

  setArrayValue(t, 1, 0, 0, 0.565685);  setArrayValue(t, 1, 0, 1, 0.424264);
  setArrayValue(t, 1, 1, 0, 0.424264);  setArrayValue(t, 1, 1, 1, 0.565685);
  setArrayValue(t, 1, 2, 0, 0.707101);  setArrayValue(t, 1, 2, 1, 0.707101);


  NewArray(&tpred, 2);
  NewArrayMatrix(&tpred, 0, 1, 2);
  NewArrayMatrix(&tpred, 1, 1, 2);

  setArrayValue(tpred, 0, 0, 0, 0.5);  setArrayValue(tpred, 0, 0, 1, 0.6);

  setArrayValue(tpred, 1, 0, 0, 0.6);  setArrayValue(tpred, 1, 0, 1, 0.4);


  puts("Array");
  PrintArray(t);


  NewUPCAModel(&m);

  UPCA(t, 2, 1, m, NULL);

  PrintUPCAModel(m);


  puts("Array to pred");
  PrintArray(tpred);

  initMatrix(&predscore);
  UPCAScorePredictor(tpred, m, 9999, &predscore);

  puts("Predicted Scores");
  PrintMatrix(predscore);

  DelMatrix(&predscore);
  DelArray(&tpred);

  DelUPCAModel(&m);
  DelArray(&t);
}


/*
 * Test with Random numbers
 */
void test2()
{
  size_t i, j, k;
  array *t;

  UPCAMODEL *m;

  puts(">>>>>>> Test 2: Compute Multi Way PCA with random variables");

  NewArray(&t, 5);

  for(i = 0; i < t->order; i++){
    NewArrayMatrix(&t, i, 30, 10);
  }

  srand(time(0));
  for(i = 0; i < t->order; i++){
    for(j = 0; j < t->m[i]->row; j++){
      for(k = 0; k < t->m[i]->col; k++)
          setArrayValue(t, i, j, k, rand()/((double)(RAND_MAX)+1));
    }
  }

  /*
  puts("Array");
  PrintArray(t);
  */

  NewUPCAModel(&m);
  ssignal run = SIGSCIENTIFICRUN;
  UPCA(t, 3, 1, m, &run);

  PrintUPCAModel(m);

  DelUPCAModel(&m);
  DelArray(&t);
}

/*Test from:
 *
 * PCA Nipals algorithm for Multi-Way Principal Component Analysis
 * Journal of Chemometrics, Vol 1, 41-56 (1987)
 * Svante Woold, Paul Geladi and Kim Esbensen
 */
void test1()
{
  array *t;
  UPCAMODEL *m;

  puts(">>>>>>> Test 1: Compute Multi Way PCA");

  NewArray(&t, 2);

  NewArrayMatrix(&t, 0, 3, 2);
  NewArrayMatrix(&t, 1, 3, 2);

  setArrayValue(t, 0, 0, 0, 0.424264);  setArrayValue(t, 0, 0, 1, 0.565685);
  setArrayValue(t, 0, 1, 0, 0.565685);  setArrayValue(t, 0, 1, 1, 0.424264);
  setArrayValue(t, 0, 2, 0, 0.707101);  setArrayValue(t, 0, 2, 1, 0.707101);

  setArrayValue(t, 1, 0, 0, 0.565685);  setArrayValue(t, 1, 0, 1, 0.424264);
  setArrayValue(t, 1, 1, 0, 0.424264);  setArrayValue(t, 1, 1, 1, 0.565685);
  setArrayValue(t, 1, 2, 0, 0.707101);  setArrayValue(t, 1, 2, 1, 0.707101);

  puts("Array");
  PrintArray(t);

  NewUPCAModel(&m);

  UPCA(t, 2, 1, m, NULL);

  PrintUPCAModel(m);

  DelUPCAModel(&m);
  DelArray(&t);
}

int main(void)
{
  test1();
  test2();
  test3();
  test4();
  return 0;
}

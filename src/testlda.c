/* testlda.c
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

#include "lda.h"
#include "numeric.h"

void Test4()
{
  size_t i, j;
  matrix *x/*, *p, *probability, *predfeatures*/;
  uivector *y/*, *classpred*/;
  LDAMODEL *lda;

  puts("Test 4: Compute LDA and Validate with Random Group Cross Validation");
  srand(150*128);
  NewMatrix(&x, 150, 128);
  NewUIVector(&y, 150);
  for(i = 0; i < 50; i++){
    for(j = 0; j < x->col; j++){
      setMatrixValue(x, i, j, randDouble(-2, 1));
    }
    y->data[i] = 0;
  }

  for(i = 50; i < 100; i++){
    for(j = 0; j < x->col; j++){
      setMatrixValue(x, i, j, randDouble(-1, 3));
    }
    y->data[i] = 1;
  }

  for(i = 100; i < 150; i++){
    for(j = 0; j < x->col; j++){
      setMatrixValue(x, i, j, randDouble(2, 5));
    }
    y->data[i] = 2;
  }

//   puts("X");PrintMatrix(x);

  NewLDAModel(&lda);
  LDA(x, y, lda);

  LDARandomGroupsCV(x, y, 5, 20, &lda->sens, &lda->spec, &lda->ppv, &lda->npv, &lda->acc, 4, NULL);
  /*LDALOOCV(x, y, &lda->sens, &lda->spec, &lda->ppv, &lda->npv, &lda->acc, 1, NULL);*/

  PrintLDAModel(lda);
  DelLDAModel(&lda);
  DelUIVector(&y);
  DelMatrix(&x);

}
void Test3()
{
  matrix *x, *pred, *probability, *mnpdf, *predfeatures;
  uivector *y, *classpred;
  LDAMODEL *lda;

  puts("Test3: Create LDA Model and Predict one external object");
  NewMatrix(&x, 7, 2);
  NewUIVector(&y, 7);

  x->data[0][0] = 2.95; x->data[0][1] = 6.63; y->data[0] = 0;
  x->data[1][0] = 2.53; x->data[1][1] = 7.79; y->data[1] = 0;
  x->data[2][0] = 3.57; x->data[2][1] = 5.65; y->data[2] = 0;
  x->data[3][0] = 3.16; x->data[3][1] = 5.47; y->data[3] = 0;

  x->data[4][0] = 2.58; x->data[4][1] = 4.46; y->data[4] = 1;
  x->data[5][0] = 2.16; x->data[5][1] = 6.22; y->data[5] = 1;
  x->data[6][0] = 3.27; x->data[6][1] = 3.52; y->data[6] = 1;

  NewMatrix(&pred, 8, 2);
  pred->data[0][0] = 2.95; pred->data[0][1] = 6.63;
  pred->data[1][0] = 2.53; pred->data[1][1] = 7.79;
  pred->data[2][0] = 3.57; pred->data[2][1] = 5.65;
  pred->data[3][0] = 3.16; pred->data[3][1] = 5.47;

  pred->data[4][0] = 2.58; pred->data[4][1] = 4.46;
  pred->data[5][0] = 2.16; pred->data[5][1] = 6.22;
  pred->data[6][0] = 3.27; pred->data[6][1] = 3.52;
  pred->data[7][0] = 2.81; pred->data[7][1] = 5.46;

  PrintMatrix(x);
  PrintUIVector(y);

  NewLDAModel(&lda);
  LDA(x, y, lda);
  PrintLDAModel(lda);

  initMatrix(&predfeatures);
  initMatrix(&probability);
  initUIVector(&classpred);
  initMatrix(&mnpdf);

  LDAPrediction(pred, lda, &predfeatures, &probability, &mnpdf, &classpred);
  puts("Probability");
  PrintMatrix(probability);

  puts("Class Prediction");
  PrintUIVector(classpred);

  puts("Predict New Features");
  PrintMatrix(predfeatures);

  puts("Predicted Multivariate Normal Profile Distribution");
  PrintMatrix(mnpdf);

  PrintUIVector(lda->classid);

  DelMatrix(&mnpdf);
  DelMatrix(&predfeatures);
  DelUIVector(&classpred);
  DelMatrix(&probability);
  DelMatrix(&pred);
  DelLDAModel(&lda);
  DelMatrix(&x);
  DelUIVector(&y);
}

void Test2()
{
  matrix *x, *pred, *probability, *mnpdf, *predfeatures;
  uivector *y, *classpred;
  LDAMODEL *lda;

  puts("Test2: Create LDA Model and Predict itself");

  NewMatrix(&x, 10, 2);
  NewUIVector(&y, 10);
  x->data[0][0] = 4; x->data[0][1] = 2; y->data[0] = 0;
  x->data[1][0] = 2; x->data[1][1] = 4; y->data[1] = 0;
  x->data[2][0] = 2; x->data[2][1] = 3; y->data[2] = 0;
  x->data[3][0] = 3; x->data[3][1] = 6; y->data[3] = 1;
  x->data[4][0] = 4; x->data[4][1] = 4; y->data[4] = 1;

  x->data[5][0] = 9; x->data[5][1] = 10; y->data[5] = 1;
  x->data[6][0] = 6; x->data[6][1] = 8; y->data[6] = 1;
  x->data[7][0] = 9; x->data[7][1] = 5; y->data[7] = 1;
  x->data[8][0] = 8; x->data[8][1] = 7; y->data[8] = 1;
  x->data[9][0] = 10; x->data[9][1] = 8; y->data[9] = 1;

  NewMatrix(&pred, 2, 2);
  pred->data[0][0] = 8; pred->data[0][1] = 6;
  pred->data[1][0] = 5.5; pred->data[1][1] = 5.5;

  PrintMatrix(x);
  PrintUIVector(y);

  NewLDAModel(&lda);
  LDA(x, y, lda);
  PrintLDAModel(lda);

  initMatrix(&predfeatures);
  initMatrix(&probability);
  initUIVector(&classpred);
  initMatrix(&mnpdf);

  LDAPrediction(pred, lda, &predfeatures, &probability, &mnpdf, &classpred);
  puts("Probability");
  PrintMatrix(probability);

  puts("Class Prediction");
  PrintUIVector(classpred);

  puts("Predict New Features");
  PrintMatrix(predfeatures);

  puts("Predicted Multivariate Normal Profile Distribution");
  PrintMatrix(mnpdf);

  DelMatrix(&predfeatures);
  DelMatrix(&pred);
  DelMatrix(&probability);
  DelUIVector(&classpred);
  DelLDAModel(&lda);
  DelMatrix(&x);
  DelUIVector(&y);
}

void Test1()
{
  matrix *x;
  uivector *y;
  LDAMODEL *lda;
  puts("Test2: Create LDA Model");
  NewMatrix(&x, 11, 2);
  NewUIVector(&y, 11);
  x->data[0][0] = 1; x->data[0][1] = 2; y->data[0] = 0;
  x->data[1][0] = 2; x->data[1][1] = 3; y->data[1] = 0;
  x->data[2][0] = 3; x->data[2][1] = 3; y->data[2] = 0;
  x->data[3][0] = 4; x->data[3][1] = 5; y->data[3] = 0;
  x->data[4][0] = 5; x->data[4][1] = 5; y->data[4] = 0;
  x->data[5][0] = 1; x->data[5][1] = 0; y->data[5] = 1;
  x->data[6][0] = 2; x->data[6][1] = 1; y->data[6] = 1;
  x->data[7][0] = 3; x->data[7][1] = 1; y->data[7] = 1;
  x->data[8][0] = 3; x->data[8][1] = 2; y->data[8] = 1;
  x->data[9][0] = 5; x->data[9][1] = 3; y->data[9] = 1;
  x->data[10][0] = 6; x->data[10][1] = 5; y->data[10] = 1;

  PrintMatrix(x);

  NewLDAModel(&lda);

  LDA(x, y, lda);

  PrintLDAModel(lda);

  DelLDAModel(&lda);
  DelMatrix(&x);
  DelUIVector(&y);
}

int main()
{
//   Test1();
    Test2();
//  Test3();
//   Test4();
  return 0;
}

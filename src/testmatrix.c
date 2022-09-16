/* testmatrix.c
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
#include <math.h>

#include "matrix.h"
#include "vector.h"
#include "numeric.h"
#include "algebra.h"


void Test34()
{

  matrix *x;
  dvector *v, *r;

  NewMatrix(&x, 7, 4);
  NewDVector(&v, 7);

  x->data[0][0] = 0.521; x->data[0][1] = 1.221; x->data[0][2] = 0.000; x->data[0][3] = -0.548;
  x->data[1][0] = -0.263; x->data[1][1] = -0.273; x->data[1][2] = 0.000; x->data[1][3] = 0.021;
  x->data[2][0] = -1.690; x->data[2][1] = -3.083; x->data[2][2] = 0.000; x->data[2][3] = 1.125;
  x->data[3][0] = -0.103; x->data[3][1] = -0.067; x->data[3][2] = 0.000; x->data[3][3] = -0.021;
  x->data[4][0] = 0.930; x->data[4][1] = 1.370; x->data[4][2] = 0.000; x->data[4][3] = -0.376;
  x->data[5][0] = -0.596; x->data[5][1] =  -0.990; x->data[5][2] = 0.000; x->data[5][3] = 0.325;
  x->data[6][0] = 1.201; x->data[6][1] = 1.822; x->data[6][2] = 0.000; x->data[6][3] = -0.525;

  v->data[0] = -3.2144;
  v->data[1] = 13.4345;
  v->data[2] = 1.2783;
  v->data[3] = 4.4551;
  v->data[4] = 47.1042;
  v->data[5] = -22.1152;
  v->data[6] = -40.9425;

  NewDVector(&r, 4);
  DVectorMatrixDotProduct(x, v, r);
  PrintDVector(r);
  DelDVector(&r);
  DelMatrix(&x);
  DelDVector(&v);
}


void Test33(){
  matrix *m;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = 3;
  m->data[0][1] = 5;
  m->data[0][2] = 6;

  m->data[1][0] = 9;
  m->data[1][1] = 2;
  m->data[1][2] = -1;

  m->data[2][0] = 5;
  m->data[2][1] = -8;
  m->data[2][2] = 1;

  double det = MatrixDeterminant(m);
  if(det == -580){
    printf("Matrix determinant %f   :OK\n", det);
  }
}

void Test32(){
  /*Large matrix*/
  matrix *mx;
  matrix *ds;
  int nrow = 356123;
  int ncol = 33;
  NewMatrix(&mx, nrow, ncol);
  int i, j;
  srand(nrow+ncol);
  for(i = 0; i < nrow; i++){
    for(j = 0; j < ncol; j++){
      mx->data[i][j] = randDouble(-1, 1);
    }
  }
  initMatrix(&ds);
  puts("Calculating Matrix Statistics Descriptive");
  MatrixColDescStat(mx, &ds);
  PrintMatrix(ds);
  DelMatrix(&ds);
  DelMatrix(&mx);
}

void Test31()
{
  matrix *m, *m_inv;
  puts("Test 31: Matrix Pseudoinversion");
  NewMatrix(&m, 6, 5);

  m->data[0][0] = 8.79; m->data[0][1] = 9.93; m->data[0][2] = 9.83; m->data[0][3] = 5.45; m->data[0][4] = 3.16;
  m->data[1][0] = 6.11; m->data[1][1] = 6.91; m->data[1][2] = 5.04; m->data[1][3] = -0.27; m->data[1][4] = 7.98;
  m->data[2][0] = -9.15; m->data[2][1] = -7.93; m->data[2][2] = 4.86; m->data[2][3] = 4.85; m->data[2][4] = 3.01;
  m->data[3][0] = 9.57; m->data[3][1] = 1.64; m->data[3][2] = 8.83; m->data[3][3] = 0.74; m->data[3][4] = 5.80;
  m->data[4][0] = -3.49; m->data[4][1] = 4.02; m->data[4][2] = 9.80; m->data[4][3] =  10.00; m->data[4][4] = 4.27;
  m->data[5][0] =  9.84; m->data[5][1] = 0.15; m->data[5][2] = -8.99; m->data[5][3] = -6.02; m->data[5][4] = -5.31;

  initMatrix(&m_inv);
  MatrixPseudoinversion(m,  &m_inv);
  PrintMatrix(m);
  PrintMatrix(m_inv);
  DelMatrix(&m);
  DelMatrix(&m_inv);
}

void Test30()
{
  matrix *mxin, *mxout;
  puts("Test 30: SVN Scaling Pretreatment");
  NewMatrix(&mxin, 3, 3);
  mxin->data[0][0] = 1; mxin->data[0][1] = 2; mxin->data[0][2] = 3;
  mxin->data[1][0] = 4; mxin->data[1][1] = 5; mxin->data[1][2] = 0;
  mxin->data[2][0] = 6; mxin->data[2][1] = 7; mxin->data[2][2] = 8;

  PrintMatrix(mxin);
  initMatrix(&mxout);
  MatrixSVNScaling(mxin,  &mxout);
  PrintMatrix(mxout);
  DelMatrix(&mxin);
  DelMatrix(&mxout);
}

void Test29()
{
  matrix *m, *t;
  puts("Test 29: Transpose a matrix of 700000 rows and 128 columns");
  NewMatrix(&m, 700000, 128);
  NewMatrix(&t, 128, 700000);
  MatrixTranspose(m, t);
  MatrixCopy(t, &m);
  DelMatrix(&t);
  DelMatrix(&m);
}

void Test28BIS()
{
  matrix *m, *U, *S, *V_T;
  puts("Test 28: Computing Singular Value Decomposition");

  /*NewMatrix(&m, 3, 3);
  m->data[0][0] = 15.629; m->data[0][1] = -6.399; m->data[0][2] = -4.843;
  m->data[1][0] = 5.305; m->data[1][1] = 14.513; m->data[1][2] = -3.903;
  m->data[2][0] = 2.822; m->data[2][1] = -1.792; m->data[2][2] = -4.535;*/

  NewMatrix(&m, 2, 2);
  m->data[0][0] = 7; m->data[0][1] = 2;
  m->data[1][0] = 3; m->data[1][1] = 4;

  puts("INIT MATRIX");
  PrintMatrix(m);

  initMatrix(&U);
  initMatrix(&S);
  initMatrix(&V_T);

  SVD(m, &U, &S, &V_T);

  puts("U"); PrintMatrix(U);
  puts("S"); PrintMatrix(S);
  puts("V_T"); PrintMatrix(V_T);

  DelMatrix(&U);
  DelMatrix(&S);
  DelMatrix(&V_T);
  DelMatrix(&m);
}

void Test28()
{
  matrix *m, *U, *S, *V_T;
  puts("Test 28: Computing Singular Value Decomposition");

  NewMatrix(&m, 5, 5);
  m->data[0][0] = 2; m->data[0][1] = 0; m->data[0][2] = 8; m->data[0][3] = 6; m->data[0][4] = 0;
  m->data[1][0] = 1; m->data[1][1] = 6; m->data[1][2] = 0; m->data[1][3] = 1; m->data[1][4] = 7;
  m->data[2][0] = 5; m->data[2][1] = 0; m->data[2][2] = 7; m->data[2][3] = 4; m->data[2][4] = 0;
  m->data[3][0] = 7; m->data[3][1] = 0; m->data[3][2] = 8; m->data[3][3] = 5; m->data[3][4] = 0;
  m->data[3][0] = 0; m->data[3][1] = 10; m->data[4][2] = 0; m->data[4][3] = 0; m->data[4][4] = 7;

  puts("INIT MATRIX");
  PrintMatrix(m);

  initMatrix(&U);
  initMatrix(&S);
  initMatrix(&V_T);

  SVD(m, &U, &S, &V_T);

  puts("U"); PrintMatrix(U);
  puts("S"); PrintMatrix(S);
  puts("V_T"); PrintMatrix(V_T);

  DelMatrix(&U);
  DelMatrix(&S);
  DelMatrix(&V_T);
  DelMatrix(&m);
}

void Test27()
{
  size_t i, j, n;
  dvector *eval;
  matrix *A, *evect;
  n = 100;
  NewMatrix(&A, n, n);
  srand(n);
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      setMatrixValue(A, i, j, randDouble(-1, 1));
    }
  }

  puts("Test 27: Calculate Eigenvectors Eigenvalue for a random matrix");
  /*
  puts("Matrix");
  PrintMatrix(A);*/
  initDVector(&eval);
  initMatrix(&evect);
  EVectEval(A, &eval, &evect);

  /*puts("Eigenvalue"); PrintDVector(eval);
  puts("Eigenvectors"); PrintMatrix(evect);*/

  DelDVector(&eval);
  DelMatrix(&evect);
  DelMatrix(&A);
}

void Test26()
{
  matrix *m;
  NewMatrix(&m, 3, 3);

  m->data[0][0] = -149; m->data[0][1] = -50; m->data[0][2] = -154;
  m->data[1][0] = 537; m->data[1][1] = 180; m->data[1][2] = 546;
  m->data[2][0] = -27; m->data[2][1] = -9; m->data[2][2] = -25;

  puts("Test 26: Householder Reduction for a matrix");
  puts("Matrix");
  PrintMatrix(m);
  HouseholderReduction(m);
  puts("Householder Hessenberg form");
  PrintMatrix(m);
  DelMatrix(&m);
}

void Test25()
{
  dvector *a, *r;
  matrix *b;
  NewDVector(&a, 3);
  NewMatrix(&b, 3, 3);
  a->data[0] = 1; a->data[1] = 2; a->data[2] = 3;

  b->data[0][0] = 1; b->data[0][1] = 0; b->data[0][2] = 0;
  b->data[1][0] = 1; b->data[1][1] = 2; b->data[1][2] = 1;
  b->data[2][0] = 0; b->data[2][1] = 1; b->data[2][2] = 3;


  puts("Test 25: Compute vector T matrix Division");
  initDVector(&r);
  DVectorTransposedMatrixDivision(a, b, r);
  puts("Vector T");
  PrintDVector(a);
  puts("Matrix");
  PrintMatrix(b);
  puts("Result");
  PrintDVector(r);
  DelDVector(&a);
  DelDVector(&r);
  DelMatrix(&b);
}

void Test24()
{
  dvector *eval;
  matrix *A, *evect;
  NewMatrix(&A, 8, 8);

  A->data[0][7] = 1;
  A->data[1][0] = 1;
  A->data[2][1] = 1;
  A->data[3][2] = 1;
  A->data[4][3] = 1;
  A->data[5][4] = 1;
  A->data[6][5] = 1;
  A->data[7][6] = 1;

  puts("Test 24: Calculate Eigenvectors Eigenvalue for a sparse matrix");

  puts("Matrix");
  PrintMatrix(A);
  initDVector(&eval);
  initMatrix(&evect);
  EVectEval(A, &eval, &evect);

  puts("Eigenvalue"); PrintDVector(eval);
  puts("Eigenvectors"); PrintMatrix(evect);

  DelDVector(&eval);
  DelMatrix(&evect);
  DelMatrix(&A);
}

void Test23()
{
  dvector *x;
  matrix *h;
  NewDVector(&x, 4);
  x->data[0] = 1;
  x->data[1] = 2;
  x->data[2] = 3;
  x->data[3] = 4;

  puts("Test 23: Build HouseholderMatrix from Vector");
  initMatrix(&h);
  HouseholderMatrix(x, h);

  PrintMatrix(h);

  DelMatrix(&h);
  DelDVector(&x);
}

void Test22()
{
  matrix *m, *L, *U, *mbis;
  NewMatrix(&m, 3, 3);
  puts("Test 22: LU Decomposition");

  m->data[0][0] = 8; m->data[0][1] = 2; m->data[0][2] = 9;
  m->data[1][0] = 4; m->data[1][1] = 9; m->data[1][2] = 4;
  m->data[2][0] = 6; m->data[2][1] = 7; m->data[2][2] = 9;

  initMatrix(&L);
  initMatrix(&U);

  LUDecomposition(m, &L, &U);

  puts("L"); PrintMatrix(L);
  puts("U"); PrintMatrix(U);

  NewMatrix(&mbis, m->row, m->col);
  MatrixDotProduct(L, U, mbis);

  PrintMatrix(m);
  PrintMatrix(mbis);

  DelMatrix(&mbis);
  DelMatrix(&L);
  DelMatrix(&U);
  DelMatrix(&m);
}

void Test21()
{
  matrix *m, *U, *S, *V_T, *US, *mbis;
  puts("Test 21: Computing Singular Value Decomposition");

  NewMatrix(&m, 4, 2);
  m->data[0][0] = 2; m->data[0][1] = 4;
  m->data[1][0] = 1; m->data[1][1] = 3;
  m->data[2][0] = 0; m->data[2][1] = 0;
  m->data[3][0] = 0; m->data[3][1] = 0;

  puts("INIT MATRIX");
  PrintMatrix(m);

  initMatrix(&U);
  initMatrix(&S);
  initMatrix(&V_T);

  SVD(m, &U, &S, &V_T);

  puts("U"); PrintMatrix(U);
  puts("S"); PrintMatrix(S);
  puts("V_T"); PrintMatrix(V_T);

  NewMatrix(&US, m->row, m->col);
  NewMatrix(&mbis, m->row, m->col);
  MatrixDotProduct(U, S, US);
  MatrixDotProduct(US, V_T, mbis);

  puts("U*S*V^T"); PrintMatrix(mbis);

  DelMatrix(&US);
  DelMatrix(&mbis);
  DelMatrix(&U);
  DelMatrix(&S);
  DelMatrix(&V_T);
  DelMatrix(&m);
}

void Test20()
{
  matrix *m, *Q, *R, *m1bis;
  puts("Test 20: QR Decomposition");

  NewMatrix(&m, 5, 3);

  m->data[0][0] = 12; m->data[0][1] = -51;  m->data[0][2] = 4;
  m->data[1][0] = 6; m->data[1][1] = 167;  m->data[1][2] = -68;
  m->data[2][0] = -4; m->data[2][1] = 24;  m->data[2][2] = -41;
  m->data[3][0] = -1; m->data[3][1] = 1;  m->data[3][2] = 0;
  m->data[4][0] = 2; m->data[4][1] = 0;  m->data[4][2] = 3;


  PrintMatrix(m);

  initMatrix(&Q);
  initMatrix(&R);
  QRDecomposition(m, &Q, &R);

  puts("Q"); PrintMatrix(Q);
  puts("R"); PrintMatrix(R);

  NewMatrix(&m1bis, Q->row, R->col);
  MatrixDotProduct(Q, R, m1bis);

  puts("m1bis"); PrintMatrix(m1bis);


  DelMatrix(&m1bis);
  DelMatrix(&R);
  DelMatrix(&Q);
  DelMatrix(&m);
}
void Test19()
{
  matrix *m, *evect, *mbis, *Q, *R;
  dvector *eval;
  puts("Test 19 Eigen Vector Eigen Value LAPACK Method and QR Decomposition");


  NewMatrix(&m, 3, 3);
  m->data[0][0] = 1; m->data[0][1] = -3; m->data[0][2] = 3;
  m->data[1][0] = 3; m->data[1][1] = -5; m->data[1][2] = 3;
  m->data[2][0] = 6; m->data[2][1] = -6; m->data[2][2] = 4;

  NewMatrix(&mbis, m->row, m->col);
  initMatrix(&Q);
  initMatrix(&R);

  QRDecomposition(m, &Q, &R);

  puts("Q"); PrintMatrix(Q);
  puts("R"); PrintMatrix(R);

  puts("Original Matrix");
  PrintMatrix(m);
  puts("Recomputed Matrix");
  MatrixDotProduct(Q, R, mbis);
  PrintMatrix(mbis);

  initDVector(&eval);
  initMatrix(&evect);
  EVectEval(m, &eval, &evect);

  puts("Eigenvalues"); PrintDVector(eval);
  puts("Eigenvectors"); PrintMatrix(evect);


  DelMatrix(&mbis);
  DelMatrix(&Q);
  DelMatrix(&R);
  DelDVector(&eval);
  DelMatrix(&evect);
  DelMatrix(&m);
}

void Test18()
{
  matrix *m, *covm;
  puts("Test Calculation Covariance Matrix");

  NewMatrix(&m, 5, 2);

  m->data[0][0] = 1; m->data[0][1] = 2;
  m->data[1][0] = 2; m->data[1][1] = 3;
  m->data[2][0] = 3; m->data[2][1] = 3;
  m->data[3][0] = 4; m->data[3][1] = 5;
  m->data[4][0] = 5; m->data[4][1] = 5;

  initMatrix(&covm);

  MatrixCovariance(m, &covm);

  puts("Matrix to for covariance matrix");
  PrintMatrix(m);

  puts("Covariance Matrix");
  PrintMatrix(covm);

  DelMatrix(&covm);
  DelMatrix(&m);
}

void Test17()
{
  size_t i, j, row, col;
  matrix *m;

  NewMatrix(&m, 10, 4);
  for(i = 0; i < 10; i++){
    for(j = 0; j < 4; j++){
      srand(time(0)+i+j+rand());
      setMatrixValue(m, i, j, randDouble(1, 24));
    }
  }

  puts("Matrix...");
  PrintMatrix(m);

  MatrixGetMinValueIndex(m, &row, &col);
  printf("Min value [%zu][%zu]: %f\n", row, col, getMatrixValue(m, row, col));

  MatrixGetMaxValueIndex(m, &row, &col);
  printf("Max value [%zu][%zu]: %f\n", row, col, getMatrixValue(m, row, col));

  DelMatrix(&m);
}

void Test16()
{
  size_t row, col;
  matrix *m, *spearmanncorrelmx;

  NewMatrix(&m, 10, 2);
  setMatrixValue(m, 0, 0, 106);  setMatrixValue(m, 0, 1, 7);
  setMatrixValue(m, 1, 0, 86);  setMatrixValue(m, 1, 1, 0);
  setMatrixValue(m, 2, 0, 100);  setMatrixValue(m, 2, 1, 27);
  setMatrixValue(m, 3, 0, 101);  setMatrixValue(m, 3, 1, 50);
  setMatrixValue(m, 4, 0, 99);  setMatrixValue(m, 4, 1, 28);
  setMatrixValue(m, 5, 0, 103);  setMatrixValue(m, 5, 1, 29);
  setMatrixValue(m, 6, 0, 97);  setMatrixValue(m, 6, 1, 20);
  setMatrixValue(m, 7, 0, 113);  setMatrixValue(m, 7, 1, 12);
  setMatrixValue(m, 8, 0, 112);  setMatrixValue(m, 8, 1, 6);
  setMatrixValue(m, 9, 0, 110);  setMatrixValue(m, 9, 1, 17);

  initMatrix(&spearmanncorrelmx);

  puts("Input matrix...");
  PrintMatrix(m);

  SpearmanCorrelMatrix(m, spearmanncorrelmx);

  puts("Spearmann Correl Matrix");
  PrintMatrix(spearmanncorrelmx);

  MatrixGetMaxValueIndex(spearmanncorrelmx, &row, &col);
  printf("Max value [%zu][%zu]: %f\n", row, col, getMatrixValue(spearmanncorrelmx, row, col));

  MatrixGetMinValueIndex(spearmanncorrelmx, &row, &col);
  printf("Min value [%zu][%zu]: %f\n", row, col, getMatrixValue(spearmanncorrelmx, row, col));

  DelMatrix(&spearmanncorrelmx);
  DelMatrix(&m);

}

void Test15()
{
  matrix *m, *pearsoncorrelmx;

  puts("Test Calculation Pearson Correlation Matrix");
  NewMatrix(&m, 6, 2);

  setMatrixValue(m, 0, 0, 24.778);   setMatrixValue(m, 0, 1, 7.0);
  setMatrixValue(m, 1, 0, 27.063);   setMatrixValue(m, 1, 1, 4.1);
  setMatrixValue(m, 2, 0, 25.588);   setMatrixValue(m, 2, 1, 5.0);
  setMatrixValue(m, 3, 0, 23.107);   setMatrixValue(m, 3, 1, 11.4);
  setMatrixValue(m, 4, 0, 14.219);   setMatrixValue(m, 4, 1, 21.5);
  setMatrixValue(m, 5, 0, 15.168);   setMatrixValue(m, 5, 1, 16.1);

  puts("Input matrix...");
  PrintMatrix(m);

  initMatrix(&pearsoncorrelmx);
  PearsonCorrelMatrix(m, pearsoncorrelmx);

  puts("Pearson Correlation Matrix");
  PrintMatrix(pearsoncorrelmx);

  DelMatrix(&pearsoncorrelmx);
  DelMatrix(&m);
}

void Test14()
{
  size_t i, j, nrow, ncol;
  matrix *m1, *m2;
  dvector *d1, *d2;

  nrow = 10;
  ncol = 2;
  NewMatrix(&m1, nrow, ncol);
  NewMatrix(&m2, nrow, ncol);

  for(i = 0; i < nrow; i++){
    for(j = 0; j < ncol; j++){
      srand(time(0)+i+j+rand());
      setMatrixValue(m1, i, j, randDouble(1, 24));
    }
  }

  initDVector(&d1);
  MatrixColAverage(m1, &d1);


  for(i = 0; i < nrow; i++){
    for(j = 0; j < ncol; j++){
      setMatrixValue(m2, i, j, log(getMatrixValue(m1, i, j)));
    }
  }

  PrintMatrix(m1);
  PrintMatrix(m2);

  for(i = 0; i < nrow; i++){
    for(j = 0; j < ncol; j++){
      printf("%f\t", exp(getMatrixValue(m2, i,j)));
    }
    printf("\n");
  }

  initDVector(&d2);
  MatrixColAverage(m2, &d2);

  for(i = 0; i < ncol; i++){
    printf("%f == %f (%f)\n", getDVectorValue(d1, i), exp(getDVectorValue(d2, i)), getDVectorValue(d2, i));
  }

  DelDVector(&d2);
  DelDVector(&d1);
  DelMatrix(&m2);
  DelMatrix(&m1);
}

void Test13()
{
  size_t i, j, nrow, ncol;
  matrix *m, *covm;
  puts("Test Calculation Covariance Matrix");

  nrow = 168;
  ncol = 128;

  NewMatrix(&m, nrow, ncol);

  for(i = 0; i < nrow; i++){
    for(j = 0; j < ncol; j++){
      srand(time(0)+i+j+rand());
      setMatrixValue(m, i, j, randDouble(-19, 24));
    }
  }

  initMatrix(&covm);

  MatrixCovariance(m, &covm);

  puts("Matrix to for covariance matrix");
  PrintMatrix(m);

  puts("Covariance Matrix");
  PrintMatrix(covm);

  DelMatrix(&covm);
  DelMatrix(&m);
}

void Test12()
{
  matrix *m, *covm;

  puts("Test Calculation Covariance Matrix");
  NewMatrix(&m, 6, 2);

  setMatrixValue(m, 0, 0, 24.778);   setMatrixValue(m, 0, 1, 7.0);
  setMatrixValue(m, 1, 0, 27.063);   setMatrixValue(m, 1, 1, 4.1);
  setMatrixValue(m, 2, 0, 25.588);   setMatrixValue(m, 2, 1, 5.0);
  setMatrixValue(m, 3, 0, 23.107);   setMatrixValue(m, 3, 1, 11.4);
  setMatrixValue(m, 4, 0, 14.219);   setMatrixValue(m, 4, 1, 21.5);
  setMatrixValue(m, 5, 0, 15.168);   setMatrixValue(m, 5, 1, 16.1);

  initMatrix(&covm);
  MatrixCovariance(m, &covm);

  puts("Covariance Matrix");
  PrintMatrix(covm);

  DelMatrix(&covm);
  DelMatrix(&m);
}

void Test11()
{
  size_t i, j, n;
  matrix *m1, *m2;
  puts("Matrix Inversion");

  n = 100;
  NewMatrix(&m1, n, n);

  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      srand(time(0)+i+j+rand());
      setMatrixValue(m1, i, j, randDouble(-19, 24));
    }
  }

  initMatrix(&m2);

  MatrixInversion(m1, &m2);

  puts("Matrix to invert");
  PrintMatrix(m1);

  puts("Matrix inverted");
  PrintMatrix(m2);

  DelMatrix(&m2);
  DelMatrix(&m1);
}

/*
   Test 10 the result of inversion matrix is:

   0.025           0.026          -0.004          -0.007          -0.003
   0.013          -0.018           0.013           0.035           0.006
   0.023          -0.031           0.008          -0.025          -0.009
   0.030           0.029          -0.048          -0.059           0.023
   0.001          -0.027           0.050           0.021           0.013

 */
void Test10()
{
  /*size_t i, j;*/
  matrix *m1, *m2;
  puts("Matrix Inversion");

  NewMatrix(&m1, 5, 5);

  setMatrixValue(m1, 0, 0, 18); setMatrixValue(m1, 0, 1, 17); setMatrixValue(m1, 0, 2, 10); setMatrixValue(m1, 0, 3, 3); setMatrixValue(m1, 0, 4, -2);
  setMatrixValue(m1, 1, 0, 21); setMatrixValue(m1, 1, 1, -12); setMatrixValue(m1, 1, 2, -14); setMatrixValue(m1, 1, 3, -2); setMatrixValue(m1, 1, 4, 5);
  setMatrixValue(m1, 2, 0, 15); setMatrixValue(m1, 2, 1, -18); setMatrixValue(m1, 2, 2, 0); setMatrixValue(m1, 2, 3, -5); setMatrixValue(m1, 2, 4, 21);
  setMatrixValue(m1, 3, 0, 1); setMatrixValue(m1, 3, 1, 21); setMatrixValue(m1, 3, 2, -8); setMatrixValue(m1, 3, 3, -4); setMatrixValue(m1, 3, 4, -8);
  setMatrixValue(m1, 4, 0, -16); setMatrixValue(m1, 4, 1, 9); setMatrixValue(m1, 4, 2, -16); setMatrixValue(m1, 4, 3, 21); setMatrixValue(m1, 4, 4, 19);

  initMatrix(&m2);

  MatrixInversion(m1, &m2);

  puts("Matrix to invert");
  PrintMatrix(m1);

  puts("Matrix inverted");
  PrintMatrix(m2);

  DelMatrix(&m2);
  DelMatrix(&m1);
}

void Test9()
{
  matrix *m1, *m2;
  puts("Matrix Inversion");

  NewMatrix(&m1, 3, 3);
  initMatrix(&m2);

  setMatrixValue(m1, 0, 0, 1.0);   setMatrixValue(m1, 0, 1, 4.0);   setMatrixValue(m1, 0, 2, 8.0);
  setMatrixValue(m1, 1, 0, 2.0);   setMatrixValue(m1, 1, 1, 1.0);   setMatrixValue(m1, 1, 2, 3.0);
  setMatrixValue(m1, 2, 0, 6.0);   setMatrixValue(m1, 2, 1, 9.0);   setMatrixValue(m1, 2, 2, 4.0);

  MatrixLUInversion(m1, &m2);

  puts("Matrix to invert");
  PrintMatrix(m1);

  puts("Matrix inverted");
  PrintMatrix(m2);

  DelMatrix(&m2);
  DelMatrix(&m1);
}

void Test8()
{
  size_t i, j;
  matrix *m;
  puts("Test 8: Use initMatrix and ResizeMatrix");

  initMatrix(&m);

  ResizeMatrix(&m, 10, 4);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      setMatrixValue(m, i, j, i+j);
    }
  }

  puts("Initial Matrix");
  PrintMatrix(m);

  puts("Resizin Matrix to 5 x 2");

  ResizeMatrix(&m, 5, 2);

  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      setMatrixValue(m, i, j, i+j);
    }
  }

  puts("Resized Matrix");
  PrintMatrix(m);

  DelMatrix(&m);
}

void Test7()
{
  dvector *colmean;
  matrix *m1, *m2;

  NewMatrix(&m1, 14, 7);
  NewMatrix(&m2, 14, 7);

  m1->data[0][0] = 4.0000;  m1->data[0][1] = 4.0000;  m1->data[0][2] = 1.0000;  m1->data[0][3] = 84.1400;  m1->data[0][4] = 1.0500;  m1->data[0][5] = 235.1500;  m1->data[0][6] = 357.1500;
  m1->data[1][0] = 5.0000;  m1->data[1][1] = 5.0000;  m1->data[1][2] = 1.0000;  m1->data[1][3] = 79.1000;  m1->data[1][4] = 0.9780;  m1->data[1][5] = 1.5090;  m1->data[1][6] = 231.0000;
  m1->data[2][0] = 4.0000;  m1->data[2][1] = 5.0000;  m1->data[2][2] = 1.0000;  m1->data[2][3] = 67.0900;  m1->data[2][4] = 0.9700;  m1->data[2][5] = 249.0000;  m1->data[2][6] = 403.0000;
  m1->data[3][0] = 4.0000;  m1->data[3][1] = 4.0000;  m1->data[3][2] = 1.0000;  m1->data[3][3] = 68.0700;  m1->data[3][4] = 0.9360;  m1->data[3][5] = 187.3500;  m1->data[3][6] = 304.5500;
  m1->data[4][0] = 3.0000;  m1->data[4][1] = 4.0000;  m1->data[4][2] = 2.0000;  m1->data[4][3] = 68.0800;  m1->data[4][4] = 1.0300;  m1->data[4][5] = 363.0000;  m1->data[4][6] = 529.0000;
  m1->data[5][0] = 9.0000;  m1->data[5][1] = 7.0000;  m1->data[5][2] = 1.0000;  m1->data[5][3] = 129.1600;  m1->data[5][4] = 1.0900;  m1->data[5][5] = 258.0000;  m1->data[5][6] = 510.0000;
  m1->data[6][0] = 10.0000;  m1->data[6][1] = 8.0000;  m1->data[6][2] = 0.0000;  m1->data[6][3] = 128.1600;  m1->data[6][4] = 1.1500;  m1->data[6][5] = 352.0000;  m1->data[6][6] = 491.0000;
  m1->data[7][0] = 6.0000;  m1->data[7][1] = 6.0000;  m1->data[7][2] = 0.0000;  m1->data[7][3] = 78.1118;  m1->data[7][4] = 0.8765;  m1->data[7][5] = 278.6400;  m1->data[7][6] = 353.3000;
  m1->data[8][0] = 16.0000;  m1->data[8][1] = 10.0000;  m1->data[8][2] = 0.0000;  m1->data[8][3] = 202.2550;  m1->data[8][4] = 1.2710;  m1->data[8][5] = 429.1500;  m1->data[8][6] = 666.6500;
  m1->data[9][0] = 6.0000;  m1->data[9][1] = 12.0000;  m1->data[9][2] = 0.0000;  m1->data[9][3] = 84.1600;  m1->data[9][4] = 0.7800;  m1->data[9][5] = 279.0000;  m1->data[9][6] = 354.0000;
  m1->data[10][0] = 4.0000;  m1->data[10][1] = 8.0000;  m1->data[10][2] = 1.0000;  m1->data[10][3] = 72.1100;  m1->data[10][4] = 0.8900;  m1->data[10][5] = 164.5000;  m1->data[10][6] = 339.0000;
  m1->data[11][0] = 4.0000;  m1->data[11][1] = 9.0000;  m1->data[11][2] = 1.0000;  m1->data[11][3] = 71.1100;  m1->data[11][4] = 0.8660;  m1->data[11][5] = 210.0000;  m1->data[11][6] = 360.0000;
  m1->data[12][0] = 5.0000;  m1->data[12][1] = 11.0000;  m1->data[12][2] = 1.0000;  m1->data[12][3] = 85.1500;  m1->data[12][4] = 0.8620;  m1->data[12][5] = 266.0000;  m1->data[12][6] = 379.0000;
  m1->data[13][0] = 5.0000;  m1->data[13][1] = 10.0000;  m1->data[13][2] = 1.0000;  m1->data[13][3] = 86.1300;  m1->data[13][4] = 0.8800;  m1->data[13][5] = 228.0000;  m1->data[13][6] = 361.0000;


  printf("Test 7: Mean Centered Matrix\n");
  printf("Original Matrix\n");
  PrintMatrix(m1);
  printf("Mean Centered\n");
  MeanCenteredMatrix(m1, m2);
  PrintMatrix(m2);


  initDVector(&colmean);

  MatrixColAverage(m2, &colmean);

  puts("Column Average");
  PrintDVector(colmean);

  DelDVector(&colmean);
  DelMatrix(&m1);
  DelMatrix(&m2);
}

void Test6()
{
  matrix *m;
  dvector *v;
  dvector *r;

  NewMatrix(&m, 3, 4);
  NewDVector(&v, 3);
  NewDVector(&r, 4);

  setMatrixValue(m, 0, 0, 3); setMatrixValue(m, 0, 1, 1); setMatrixValue(m, 0, 2, 0); setMatrixValue(m, 0, 3, 3);
  setMatrixValue(m, 1, 0, 2); setMatrixValue(m, 1, 1, 4); setMatrixValue(m, 1, 2, 7); setMatrixValue(m, 1, 3, 0);
  setMatrixValue(m, 2, 0, -1); setMatrixValue(m, 2, 1, 3); setMatrixValue(m, 2, 2, 3); setMatrixValue(m, 2, 3, 4);

  setDVectorValue(v, 0, 1); setDVectorValue(v, 1, 2); setDVectorValue(v, 2, 3);

  DVectorMatrixDotProduct(m, v, r);

  printf("Test 6: Product Matrix x Col Vector \n");
  printf("The Matrix\n");
  PrintMatrix(m);
  printf("The Column Vector\n");
  PrintDVector(v);
  printf("The result must be: 4, 18, 23, 15\n");

  printf("The computed result is\n");
  PrintDVector(r);

  DelDVector(&r);
  DelDVector(&v);
  DelMatrix(&m);
}


void Test5_bis()
{
  matrix *m;
  dvector *v;
  dvector *r;
  dvector *r_bis;
  size_t i, j, rowsize = 25000, colsize = 5000;
  clock_t begin, end;
  double time_spent;

  NewMatrix(&m, rowsize, colsize);
  NewDVector(&v, colsize);
  NewDVector(&r, rowsize);
  NewDVector(&r_bis, rowsize);

  srand(rowsize+colsize);
  for(i = 0; i < rowsize; i++){
    for(j = 0; j < colsize; j++){
      m->data[i][j] = rand() % 100;
    }
  }

  for(j = 0; j < colsize; j++){
    v->data[j] = rand() % 25;
  }


  begin = clock();
  MatrixDVectorDotProduct(m, v, r);
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("time_spent(MatrixDVectorDotProduct) %f\n", time_spent);

  begin = clock();
  MT_MatrixDVectorDotProduct(m, v, r_bis);
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("time_spent(MultiThreadMatrixDVectorDotProduct) %f\n", time_spent);

  printf("Test 5 BIS: Product Matrix x Row Vector \n");
  puts("Check the result...");
  for(i = 0; i < r->size; i++){
    double res = r->data[i] - r_bis->data[i];
    if(FLOAT_EQ(res, 0.00, 1e-2))
      continue;
    else{
      printf("%f != %f\n", r->data[i], r_bis->data[i]);
    }
  }
  puts("Done!");

  DelDVector(&r_bis);
  DelDVector(&r);
  DelDVector(&v);
  DelMatrix(&m);
}


void Test5()
{
  matrix *m;
  dvector *v;
  dvector *r;
  dvector *r_bis;

  NewMatrix(&m, 2, 4);
  NewDVector(&v, 4);
  NewDVector(&r, 2);
  NewDVector(&r_bis, 2);

  setMatrixValue(m, 0, 0, -3); setMatrixValue(m, 0, 1, 0); setMatrixValue(m, 0, 2, 3); setMatrixValue(m, 0, 3, 2);
  setMatrixValue(m, 1, 0, 1); setMatrixValue(m, 1, 1, 7); setMatrixValue(m, 1, 2, -1); setMatrixValue(m, 1, 3, 9);

  setDVectorValue(v, 0, 2); setDVectorValue(v, 1, -3); setDVectorValue(v, 2, 4); setDVectorValue(v, 3, -1);

  MatrixDVectorDotProduct(m, v, r);

  MT_MatrixDVectorDotProduct(m, v, r_bis);

  printf("Test 5: Product Matrix x Row Vector \n");
  printf("The Matrix\n");
  PrintMatrix(m);
  printf("The Row Vector\n");
  PrintDVector(v);
  printf("The result must be: 4 and -32\n");

  printf("The computed result is\n");
  PrintDVector(r);

  printf("The computed result is\n");
  PrintDVector(r_bis);

  DelDVector(&r_bis);
  DelDVector(&r);
  DelDVector(&v);
  DelMatrix(&m);
}

void Test4()
{
  size_t i, j;
  matrix *m;
  matrix *t;
  NewMatrix(&m, 10, 15);
  NewMatrix(&t, 15, 10);

  puts("Test 4");
  srand(time(0));
  for(i = 0; i < 10; i++){
    for(j = 0; j < 15; j++){
      setMatrixValue(m, i, j, i+j);
    }
  }
  printf("Original Matrix\n");
  PrintMatrix(m);

  MatrixTranspose(m, t);

  printf("Transposed matrix\n");
  PrintMatrix(t);

  DelMatrix(&m);
  DelMatrix(&t);
}

void Test3()
{
  long int i, j;
  matrix *m;
  dvector *col;

  srand(time(0));

  printf("Test 3 Appending Column\n");

  NewMatrix(&m, 10, 15);

  for(i = 0; i < 10; i++){
    for(j = 0; j < 15; j++){
      m->data[i][j] = i+j;
    }
  }

  NewDVector(&col, 14);
  for(i = 0; i < col->size; i++){
    col->data[i] = rand() % 50;
  }

  MatrixAppendCol(&m, col);

  DelDVector(&col);

  printf("Finale Output\n");
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      printf("%.4f\t", m->data[i][j]);
    }
    printf("\n");
  }

  DelMatrix(&m);
}

void Test2()
{
  long int i, j;
  matrix *m;
  dvector *row;

  srand(time(0));

  printf("Test 2: Appending Row\n");

  NewMatrix(&m, 10, 15);

  for(i = 0; i < 10; i++){
    for(j = 0; j < 15; j++){
      m->data[i][j] = i+j;
    }
  }

  NewDVector(&row, 15);
  for(i = 0; i < row->size; i++){
    row->data[i] = rand() % 50;
  }


  MatrixAppendRow(&m, row);

  DelDVector(&row);

  printf("Finale Output\n");
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      printf("%.4f\t", m->data[i][j]);
    }
    printf("\n");
  }

  DelMatrix(&m);
}

void Test1()
{
  double a = sqrt(-1);
  double b = 1/2.;
  double c = 1/0.;
  long int i, j;
  matrix *m;

  printf("Test 1: Create a Matrix\n");
  NewMatrix(&m, 10, 15);
  for(i = 0; i < 10; i++){
    for(j = 0; j < 15; j++){
      m->data[i][j] = i+j;
    }
  }

  printf("1/0. isnan %d -> %f\n", _isnan_(a), a);
  printf("1/2. isnan %d -> %f\n", _isnan_(b), b);
  printf("1/2. isinf %d -> %f\n", _isinf_(c), c);
  printf("1/2. isinf %d -> %f\n", _isinf_(b), b);

  printf("Finale Output\n");
  PrintMatrix(m);

  DelMatrix(&m);
}

int main(void)
{
  /*Test1();
  Test2();
  Test3();
  Test4();
  Test5();
  Test5_bis();*/
  /*Test6();
  Test7();
  Test8();
  Test9();*/
  /*Test10();
  Test11();
  Test12();
  Test13();
  Test14();
  Test15();
  Test16();*/
  /*Test17();*/
  Test18();
  /*Test19();*/
  /*Test20();
  Test21();*/
  /*Test22();
  Test23();
  Test24();*/
  /*Test25();
  Test26();
  Test27();
  Test28();*/
  Test28BIS();
  /*Test29();
  Test30();
  Test31();
  Test32();
  Test33();*/
  Test34();
  return 0;
}

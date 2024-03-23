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

/*
MISSING TEST! Matrix2IntFactorsMatrix
*/

void Test62()
{
  puts("Test 62: setMatrixValue/getMatrixValue");
  matrix *m;
  NewMatrix(&m, 2,2);
  MatrixSet(m, 1);
  setMatrixValue(m, 0, 0, 0.0/0.0);
  setMatrixValue(m, 1, 1, 1.0/0.0);
  setMatrixValue(m, 3,3, 1); // Error
  getMatrixValue(m, 3,3); // Error
  if(!FLOAT_EQ(m->data[0][0], MISSING, 1e-1) ||
     !FLOAT_EQ(m->data[1][1], MISSING, 1e-1)){
    puts("Error Test62");
    abort();
  }
  else{
    puts("setMatrixValue: OK");
  }
  DelMatrix(&m);
}

void Test61()
{
  /*Large matrix*/
  puts("Test 61: MatrixDotProduct");
  size_t i;
  size_t j;
  matrix *a;
  matrix *b;
  matrix *r1;
  matrix *r2;
  int nrow = 500;
  int ncol = 500;
  NewMatrix(&a, nrow, ncol);
  NewMatrix(&b, nrow, ncol);
  NewMatrix(&r1, nrow, ncol);
  NewMatrix(&r2, nrow, ncol);
  srand(time(0));
  for(i = 0; i < nrow; i++){
    for(j = 0; j < ncol; j++){
      a->data[i][j] = randDouble(-1, 1);
    }
  }
  for(i = 0; i < nrow; i++){
    for(j = 0; j < ncol; j++){
      b->data[i][j] = randDouble(-1, 1);
    }
  }
  MatrixDotProduct_(a, b, r1);
  MatrixDotProduct_LOOP_UNROLLING(a, b, r2);
 for(i = 0; i < r1->row; i++){
    for(j = 0; j < r1->col; j++){
      if(FLOAT_EQ(r1->data[i][j], r2->data[i][j], 1e-12)){
        continue;
      }
      else{
        puts("Error Test61");
        abort();
      }
    }
  }
  puts("MatrixDotProduct_ and MatrixDotProduct_LOOP_UNROLLING: OK");
  DelMatrix(&a);
  DelMatrix(&b);
  DelMatrix(&r1);
  DelMatrix(&r2);
}

void Test60()
{
  puts("Test 60: RowColOuterProduct");
  dvector *a;
  dvector *b;
  matrix *m;
  size_t i;
  size_t j;
  NewDVector(&a, 3);
  NewDVector(&b, 4);
  DVectorSet(a, 2.f);
  DVectorSet(b, 4.f);
  NewMatrix(&m, 3, 4);
  RowColOuterProduct(a, b, m);
  for(i = 0; i < 3; i++){
    for(j = 0; j < 4; j++){
      if(FLOAT_EQ(m->data[i][j], 8.f, 1e-1)){
        continue;
      }
      else{
        puts("Error Test60");
        abort();
      }
    }
  }
  puts("RowColOuterProduct: OK");
  DelDVector(&a);
  DelDVector(&b);
  DelMatrix(&m);
}

void Test59()
{
  puts("Test 56: MatrixSort");
  matrix *m;
  size_t i;
  NewMatrix(&m, 3, 2);
  m->data[0][1] = 3;
  m->data[1][1] = 2;
  m->data[2][1] = 1;
  MatrixSort(m, 1);
  for(i = 0; i < m->row-1; i++){
    if(m->data[i][1] < m->data[i+1][1]){
      continue;
    }
    else{
      puts("Error Test59");
      abort();
    }
  }
  puts("MatrixSort: OK");
  DelMatrix(&m);
}

void Test58()
{
  puts("Test 58: MatrixRowCenterScaling");
  matrix *m;
  matrix *mv;
  NewMatrix(&m, 3, 2);
  m->data[0][0] = -1; m->data[0][1] = 2;
  m->data[1][0] = 3; m->data[1][1] = -2;
  m->data[2][0] = 5; m->data[2][1] = 7;
  initMatrix(&mv);
  MatrixRowCenterScaling(m, mv);
  if(!FLOAT_EQ(mv->data[0][0], -1, 1e-2) || 
     !FLOAT_EQ(mv->data[0][1],  2, 1e-2) ||
     !FLOAT_EQ(mv->data[1][0], 3, 1e-2) || 
     !FLOAT_EQ(mv->data[1][1], -2, 1e-2) ||
     !FLOAT_EQ(mv->data[2][0], 0.417, 1e-2) || 
     !FLOAT_EQ(mv->data[2][1], 0.583, 1e-2)){
    puts("Error Test58");
    abort();
  }
  puts("MatrixRowCenterScaling: OK");
  DelMatrix(&mv);
  DelMatrix(&m);
} 

void Test57()
{
  puts("Test 57: MatrixRowAverage");
  matrix *m;
  dvector *v;
  NewMatrix(&m, 3, 2);
  m->data[0][0] = -1; m->data[0][1] = 2;
  m->data[1][0] = 3; m->data[1][1] = -2;
  m->data[2][0] = 5; m->data[2][1] = 7;
  initDVector(&v);
  MatrixRowAverage(m, v);
  if(!FLOAT_EQ(v->data[0], 0.5, 1e-2) || 
     !FLOAT_EQ(v->data[1], 0.5, 1e-2) ||
     !FLOAT_EQ(v->data[2], 6, 1e-2)){
    puts("Error Test57");
    abort();
  }
  puts("MatrixRowAverage: OK");
  DelDVector(&v);
  DelMatrix(&m);
} 

void Test56()
{
  puts("Test 56: MatrixReverseSort");
  matrix *m;
  size_t i;
  NewMatrix(&m, 3, 2);
  m->data[0][1] = 1;
  m->data[1][1] = 2;
  m->data[2][1] = 3;
  MatrixReverseSort(m, 1);
  for(i = 0; i < m->row-1; i++){
    if(m->data[i][1] > m->data[i+1][1]){
      continue;
    }
    else{
      puts("Error Test56");
      abort();
    }
  }
  puts("MatrixReverseSort: OK");
  DelMatrix(&m);
}

void Test55()
{
  puts("Test 55: MatrixNorm");
  matrix *m;
  matrix *mn;
  NewMatrix(&m, 3, 2);
  m->data[0][0] = -1; m->data[0][1] = 2;
  m->data[1][0] = 3; m->data[1][1] = -2;
  m->data[2][0] = 5; m->data[2][1] = 7;
  NewMatrix(&mn, 3, 2);
  MatrixNorm(m, mn);
  if(!FLOAT_EQ(mn->data[0][0], -0.104, 1e-2) || 
     !FLOAT_EQ(mn->data[0][1],  0.209, 1e-2) ||
     !FLOAT_EQ(mn->data[1][0], 0.313, 1e-2) || 
     !FLOAT_EQ(mn->data[1][1], -0.209, 1e-2) ||
     !FLOAT_EQ(mn->data[2][0], 0.521, 1e-2) || 
     !FLOAT_EQ(mn->data[2][1], 0.730, 1e-2)){
    puts("Error Test55");
    abort();
  }
  DelMatrix(&m);
  DelMatrix(&mn);
  puts("MatrixNorm: OK");
}

void Test54()
{
  puts("Test 54: MatrixMoorePenrosePseudoinverse");
  matrix *m;
  matrix *inv;
  NewMatrix(&m, 3, 2);
  m->data[0][0] = -1; m->data[0][1] = 2;
  m->data[1][0] = 3; m->data[1][1] = -2;
  m->data[2][0] = 5; m->data[2][1] = 7;
  initMatrix(&inv);
  MatrixMoorePenrosePseudoinverse(m, inv);
  if(!FLOAT_EQ(inv->data[0][0], -0.088, EPSILON) ||
     !FLOAT_EQ(inv->data[0][1], 0.178, EPSILON) ||
     !FLOAT_EQ(inv->data[0][2], 0.076, EPSILON) ||
     !FLOAT_EQ(inv->data[1][0], 0.077, EPSILON) ||
     !FLOAT_EQ(inv->data[1][1], -0.119, EPSILON) ||
     !FLOAT_EQ(inv->data[1][2], 0.087, EPSILON)){
    puts("Error Test54");
    abort();
  }
  DelMatrix(&m);
  DelMatrix(&inv);
  puts("MatrixMoorePenrosePseudoinverse: OK");
}

void Test53()
{
  puts("Test 53: MatrixInitRandomInt");
  matrix *m;
  size_t i;
  size_t j;
  NewMatrix(&m, 3, 3);
  MatrixInitRandomInt(m, -100, 100);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if((int)m->data[i][j] >= -100 && (int)m->data[i][j] <= 100 && !FLOAT_EQ(m->data[i][j], 0.f, 1e-8)){
        continue;
      }
      else{
        printf("%f \n", m->data[i][j]);
        puts("Error Test53");
        abort();
      }
    }
  }
  DelMatrix(&m);
  puts("MatrixInitRandomInt: OK");
}

void Test52()
{
  puts("Test 52: MatrixInitRandomFloat");
  matrix *m;
  size_t i;
  size_t j;
  NewMatrix(&m, 3, 3);
  MatrixInitRandomFloat(m, -2.f, 2.f);

  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(m->data[i][j] > -2.f && m->data[i][j] < 2.f && !FLOAT_EQ(m->data[i][j], 0.f, 1e-8)){
        continue;
      }
      else{
        puts("Error Test52");
        abort();
      }
    }
  }
  DelMatrix(&m);
  puts("MatrixInitRandomFloat: OK");
}

void Test51()
{
  puts("Test 51: MatrixDeleteRowAt");
  matrix *m;
  size_t i;
  size_t j;
  NewMatrix(&m, 3, 3);
  m->data[1][0] = 1; m->data[1][1] = 1; m->data[1][2] = 1;
  MatrixDeleteRowAt(m, 1);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(!FLOAT_EQ(m->data[i][j], 1.f, 1e-1)){
        continue;
      }
      else{
        puts("Error Test51");
        abort();
      }
    }
  }
  DelMatrix(&m);
  puts("MatrixDeleteRowAt: OK");
}


void Test50()
{
  puts("Test 50: MatrixDeleteColAt");
  matrix *m;
  size_t i;
  size_t j;
  NewMatrix(&m, 3, 3);
  m->data[0][1] = 1; m->data[1][1] = 1; m->data[2][1] = 1;
  MatrixDeleteColAt(m, 1);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(!FLOAT_EQ(m->data[i][j], 1.f, 1e-1)){
        continue;
      }
      else{
        puts("Error Test50");
        abort();
      }
    }
  }
  DelMatrix(&m);
  puts("MatrixDeleteColAt: OK");
}

void Test49()
{
  puts("Test 48: MatrixCheck");
  matrix *m;
  NewMatrix(&m, 3, 3);
  m->data[1][1] = NAN;
  m->data[0][2] = NAN;
  m->data[2][1] = NAN;
  MatrixCheck(m);
   if(!FLOAT_EQ(m->data[1][1], MISSING, 1e-1) ||
      !FLOAT_EQ(m->data[0][2], MISSING, 1e-1) ||
      !FLOAT_EQ(m->data[2][1], MISSING, 1e-1)){
    puts("Error Test49");
    abort();
  }
  DelMatrix(&m);
  puts("MatrixCheck: OK");
}

void Test48()
{
  puts("Test 48: MatrixAppendUIRow");
  matrix *m;
  uivector *v;
  size_t j;
  NewMatrix(&m, 3, 3);
  NewUIVector(&v, 3);
  v->data[0] = 0; v->data[1] = 1; v->data[2] = 2; 
  MatrixAppendUIRow(m, v);
  for(j = 0; j < m->col; j++){
    if(FLOAT_EQ(m->data[m->row-1][j], v->data[j], EPSILON)){
      continue;
    }
    else{
      puts("Error Test48");
      abort();
    }
  }
  DelUIVector(&v);
  DelMatrix(&m);
  puts("MatrixAppendUIRow: OK");
}

void Test47()
{
  puts("Test 47: MatrixSet");
  matrix *m;
  size_t i;
  size_t j;
  NewMatrix(&m, 3, 3);
  MatrixSet(m, -1);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(!FLOAT_EQ(m->data[i][j], -1.f, 1e-3)){
        puts("Error Test47");
        abort();
      }
      else
        continue;
    }
  }
  puts("MatrixSet: OK");
  DelMatrix(&m);
}

void Test46()
{
  puts("Test 46: MatrixColumnMinMax");
  matrix *m;
  double min;
  double max;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = 1; m->data[0][1] = 0; m->data[0][2] = 3;
  m->data[1][0] = 11.; m->data[1][1] = 5; m->data[1][2] = 2;
  m->data[2][0] = 6; m->data[2][1] = 12; m->data[2][2] = 5;
  
  MatrixColumnMinMax(m, 0, &min, &max);
  if(!FLOAT_EQ(min, 1.f, 1e-4) ||
     !FLOAT_EQ(max, 11.f, 1e-4)){
    puts("Error Test46");
    abort();
  }
  MatrixColumnMinMax(m, 1, &min, &max);
  if(!FLOAT_EQ(min, 0.f, 1e-4) ||
     !FLOAT_EQ(max, 12.f, 1e-4)){
    puts("Error Test46");
    abort();
  }
  MatrixColumnMinMax(m, 1, &min, &max);
  if(!FLOAT_EQ(min, 0.f, 1e-4) ||
     !FLOAT_EQ(max, 12.f, 1e-4)){
    puts("Error Test46");
    abort();
  }
  DelMatrix(&m);
  puts("MatrixColumnMinMax: OK");
}

void Test45()
{
  puts("Test 45: MatrixColVar");
  matrix *m;
  dvector *v;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = 1; m->data[0][1] = 0; m->data[0][2] = 3;
  m->data[1][0] = 11.; m->data[1][1] = 5; m->data[1][2] = 2;
  m->data[2][0] = 6; m->data[2][1] = 12; m->data[2][2] = 5;
  initDVector(&v);
  MatrixColVar(m, v);
  if(!FLOAT_EQ(v->data[0], 25.0000, 1e-4) ||
     !FLOAT_EQ(v->data[1], 36.3333, 1e-4) ||
     !FLOAT_EQ(v->data[2], 2.3333, 1e-4)){
    puts("Error Test45");
    abort();
  }
  DelMatrix(&m);
  DelDVector(&v);
  puts("MatrixColVar: OK");
}

void Test44()
{
  puts("Test 44: MatrixColSDEV");
  matrix *m;
  dvector *v;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = 1; m->data[0][1] = 0; m->data[0][2] = 3;
  m->data[1][0] = 11.; m->data[1][1] = 5; m->data[1][2] = 2;
  m->data[2][0] = 6; m->data[2][1] = 12; m->data[2][2] = 5;
  initDVector(&v);
  MatrixColSDEV(m, v);
  if(!FLOAT_EQ(v->data[0], 5.0000, 1e-4) ||
     !FLOAT_EQ(v->data[1], 6.0277, 1e-4) ||
     !FLOAT_EQ(v->data[2], 1.5275, 1e-4)){
    puts("Error Test44");
    abort();
    }
  DelMatrix(&m);
  DelDVector(&v);
  puts("MatrixColSDEV: OK");
}

void Test43()
{
  puts("Test 43: MatrixColRMS");
  matrix *m;
  dvector *v;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = 1; m->data[0][1] = 0; m->data[0][2] = 3;
  m->data[1][0] = 11.; m->data[1][1] = 5; m->data[1][2] = 2;
  m->data[2][0] = 6; m->data[2][1] = 12; m->data[2][2] = 5;
  initDVector(&v);
  MatrixColRMS(m, v);
  if(!FLOAT_EQ(v->data[0], 7.2572, 1e-3) ||
     !FLOAT_EQ(v->data[1], 7.5056, 1e-3) ||
     !FLOAT_EQ(v->data[2], 3.5590, 1e-3)){
    puts("Error Test43");
    abort();
  }
  DelMatrix(&m);
  DelDVector(&v);
  puts("MatrixColRMS: OK");
}

void Test42()
{
  puts("Test 42: Matrix2SquareMatrix");
  matrix *m;
  matrix *o;
  size_t i;
  size_t j;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = 1; m->data[0][1] = 0; m->data[0][2] = 3;
  m->data[1][0] = 11.; m->data[1][1] = 5; m->data[1][2] = 2;
  m->data[2][0] = 6; m->data[2][1] = 12; m->data[2][2] = 5;
  NewMatrix(&o, 3, 3);
  Matrix2SquareMatrix(m, o);
  for(i = 0; i < m->row; i++){
    for( j = 0; j < m->col; j++){
      if(FLOAT_EQ(o->data[i][j], square(m->data[i][j]), 1e-1))
        continue;
      else{
        puts("Error Test42");
        abort();
      }
    }
  }
  DelMatrix(&m);
  DelMatrix(&o);
  puts("Matrix2SquareMatrix: OK");
}

void Test41()
{
  puts("Test 40: Matrix2SQRTMatrix");
  matrix *m;
  matrix *o;
  size_t i;
  size_t j;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = 1; m->data[0][1] = 0; m->data[0][2] = 3;
  m->data[1][0] = 11.; m->data[1][1] = 5; m->data[1][2] = 2;
  m->data[2][0] = 6; m->data[2][1] = 12; m->data[2][2] = 5;
  NewMatrix(&o, 3, 3);
  Matrix2SQRTMatrix(m, o);
  for(i = 0; i < m->row; i++){
    for( j = 0; j < m->col; j++){
      if(FLOAT_EQ(o->data[i][j], sqrt(m->data[i][j]), 1e-1))
        continue;
      else{
        puts("Error Test41");
        abort();
      }
    }
  }
  DelMatrix(&m);
  DelMatrix(&o);
  puts("Matrix2SQRTMatrix: OK");
}

void Test40()
{
  puts("Test 40: Matrix2LogMatrix");
  matrix *m;
  matrix *o;
  size_t i;
  size_t j;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = 1; m->data[0][1] = 0; m->data[0][2] = 3;
  m->data[1][0] = 11.; m->data[1][1] = 5; m->data[1][2] = 2;
  m->data[2][0] = 6; m->data[2][1] = 12; m->data[2][2] = 5;
  NewMatrix(&o, 3, 3);
  Matrix2LogMatrix(m, o);
  for(i = 0; i < m->row; i++){
    for( j = 0; j < m->col; j++){
      if(FLOAT_EQ(o->data[i][j], log10(m->data[i][j]+1), 1e-1))
        continue;
      else{
        puts("Error Test40");
        abort();
      }
    }
  }
  DelMatrix(&m);
  DelMatrix(&o);
  puts("Matrix2LogMatrix: OK");
}

void Test39()
{
  puts("Test 39: Matrix2ABSMatrix");
  matrix *m;
  matrix *o;
  size_t i;
  size_t j;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = -1; m->data[0][1] = 0; m->data[0][2] = -3;
  m->data[1][0] = -11.; m->data[1][1] = -5; m->data[1][2] = -2;
  m->data[2][0] = -6; m->data[2][1] = -12; m->data[2][2] = -5;
  NewMatrix(&o, 3, 3);
  Matrix2ABSMatrix(m, o);
  for(i = 0; i < m->row; i++){
    for( j = 0; j < m->col; j++){
      if(FLOAT_EQ(o->data[i][j], fabs(m->data[i][j]), 1e-1))
        continue;
      else{
        puts("Error Test39");
        abort();
      }
    }
  }
  DelMatrix(&m);
  DelMatrix(&o);
  puts("Matrix2ABSMatrix: OK");
}
void Test38()
{
  puts("Test38: GenIdentityMatrix");
  matrix *m;
  size_t i;
  size_t j;
  NewMatrix(&m, 10, 10);
  GenIdentityMatrix(m);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(i != j){
        if(FLOAT_EQ(m->data[i][j], 0.f, 1e-2)){
          continue;
        }
        else{
          puts("Error Test38");
          abort();
        }
      }
      else{
        if(FLOAT_EQ(m->data[i][j], 1.f, 1e-2)){
          continue;
        }
        else{
          puts("Error Test38");
          abort();
        }
      }
    }
  }
  DelMatrix(&m);
  puts("GenIdentityMatrix: OK");
}

void Test37()
{
  puts("Test37: FindNan");
  matrix *m;
  NewMatrix(&m, 3, 3);
  m->data[1][1] = NAN;
  m->data[0][2] = NAN;
  m->data[2][1] = NAN;
  FindNan(m);
  DelMatrix(&m);
  puts("FindNan: OK");
}

void Test36()
{
  puts("Test36: DVectorTrasposedDVectorDotProduct");
  dvector *v1;
  dvector *v2;
  matrix *m;
  size_t i;
  size_t j;
  NewDVector(&v1, 5);
  NewDVector(&v2, 10);
  DVectorSet(v1, 2);
  DVectorSet(v2, 3);
  initMatrix(&m);
  DVectorTrasposedDVectorDotProduct(v1, v2, m);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(FLOAT_EQ(m->data[i][j], 6.f, 1e-2)){
        continue;
      }
      else{
        puts("Error Test36");
        abort();
      }
    }
  }
  DelDVector(&v1);
  DelDVector(&v2);
  DelMatrix(&m);
  puts("DVectorTrasposedDVectorDotProduct: OK");
}
void Test35()
{
  puts("Test 35: ValInMatrix");
  matrix *m;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = 1.23; m->data[0][1] = 0.412; m->data[0][2] = 3.123;
  m->data[1][0] = 11.95; m->data[1][1] = 5.14; m->data[1][2] = 2.3112313234234123;
  m->data[2][0] = 6.32; m->data[2][1] = 12.42; m->data[2][2] = -5.423;
  
  if(ValInMatrix(m, 2.3112313234234123) == 1 &&
     ValInMatrix(m, 2.3212313234234123) == 0){
    puts("ValInMatrix: OK");
  }
  else{
    puts("Error Test35");
    abort();
  }
  DelMatrix(&m);
}

void Test34()
{
  puts("Test 34: Calculate MatrixTrace");
  matrix *m;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = 1; m->data[0][1] = 0; m->data[0][2] = 3;
  m->data[1][0] = 11; m->data[1][1] = 5; m->data[1][2] = 2;
  m->data[2][0] = 6; m->data[2][1] = 12; m->data[2][2] = -5;
  if(!FLOAT_EQ(1., MatrixTrace(m), 1e-3)){
    puts("Error Test34");
    abort();
  }
  DelMatrix(&m);
  puts("Calculate MatrixTrace: OK.\n");
}

void Test33()
{
  puts("Test 33: Calculate Matrixnorm");
  matrix *m;
  NewMatrix(&m, 3, 3);
  m->data[0][0] = -4; m->data[0][1] = -3; m->data[0][2] = -2;
  m->data[1][0] = -1; m->data[1][1] = 0; m->data[1][2] = 1;
  m->data[2][0] = 2; m->data[2][1] = 3; m->data[2][2] = 4;
  if(!FLOAT_EQ(7.745966692414834, Matrixnorm(m), 1e-13)){
    puts("Error Test33");
    abort();
  }
  DelMatrix(&m);
  puts("Calculate Matrixnorm: OK.\n");
}

void Test32()
{
  puts("Test 32: Calculate DVectorMatrixDotProduct");
  matrix *x;
  dvector *v;
  dvector *r1;
  dvector *r2;
  size_t i;

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

  NewDVector(&r1, 4);
  DVectorMatrixDotProduct(x, v, r1);

  NewDVector(&r2, 4);
  MT_DVectorMatrixDotProduct(x, v, r2);

  for(i = 0; i < r1->size; i++){
    if(FLOAT_EQ(r1->data[i], r2->data[i], 1e-8)){
      continue;
    }
    else{
      puts("Error Test32");
      abort();
    }
  }
  DelDVector(&r1);
  DelDVector(&r2);
  DelMatrix(&x);
  DelDVector(&v);
  puts("Calculate DVectorMatrixDotProduct: OK");
}


void Test31(){
  puts("Test 31: Calculate the matrix determinant");
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
  DelMatrix(&m);
}

void Test30(){
  /*Large matrix*/
  puts("Test 30: Calculate descriptive statistics from a matrix");
  srand_(123456);
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
  MatrixColDescStat(mx, ds);
  DelMatrix(&ds);
  DelMatrix(&mx);
  puts("Calculating Matrix Statistics Descriptive: OK");
}

void Test29()
{
  matrix *m, *m_inv;
  puts("Test 29: Matrix Pseudoinversion");
  NewMatrix(&m, 5, 5);

  m->data[0][0] = 8.79; m->data[0][1] = 9.93; m->data[0][2] = 9.83; m->data[0][3] = 5.45; m->data[0][4] = 3.16;
  m->data[1][0] = 6.11; m->data[1][1] = 6.91; m->data[1][2] = 5.04; m->data[1][3] = -0.27; m->data[1][4] = 7.98;
  m->data[2][0] = -9.15; m->data[2][1] = -7.93; m->data[2][2] = 4.86; m->data[2][3] = 4.85; m->data[2][4] = 3.01;
  m->data[3][0] = 9.57; m->data[3][1] = 1.64; m->data[3][2] = 8.83; m->data[3][3] = 0.74; m->data[3][4] = 5.80;
  m->data[4][0] = -3.49; m->data[4][1] = 4.02; m->data[4][2] = 9.80; m->data[4][3] =  10.00; m->data[4][4] = 4.27;

  initMatrix(&m_inv);
  MatrixPseudoinversion(m,  m_inv);
  PrintMatrix(m);
  PrintMatrix(m_inv);
  DelMatrix(&m);
  DelMatrix(&m_inv);
}

void Test28()
{
  matrix *mxin, *mxout;
  puts("Test 28: SVN Scaling Pretreatment");
  NewMatrix(&mxin, 3, 3);
  mxin->data[0][0] = 1; mxin->data[0][1] = 2; mxin->data[0][2] = 3;
  mxin->data[1][0] = 4; mxin->data[1][1] = 5; mxin->data[1][2] = 0;
  mxin->data[2][0] = 6; mxin->data[2][1] = 7; mxin->data[2][2] = 8;

  PrintMatrix(mxin);
  initMatrix(&mxout);
  MatrixSVNScaling(mxin,  mxout);
  PrintMatrix(mxout);
  DelMatrix(&mxin);
  DelMatrix(&mxout);
}

void Test27()
{
  matrix *m, *t;
  puts("Test 27: Transpose a matrix of 7000 rows and 128 columns");
  NewMatrix(&m, 7000, 128);
  NewMatrix(&t, 128, 7000);
  MatrixTranspose(m, t);
  MatrixCopy(t, &m);
  DelMatrix(&t);
  DelMatrix(&m);
}
void Test26_tris()
{
  matrix *m, *U, *S, *V_T;

  puts("Test 26 tris: Computing Singular Value Decomposition Lapack method");

  NewMatrix(&m, 6, 4);
  m->data[0][0] = 7.52; m->data[0][1] = -1.10; m->data[0][2] = -7.95; m->data[0][3] = 1.08;
  m->data[1][0] = -0.76; m->data[1][1] = 0.62; m->data[1][2] = 9.34; m->data[1][3] = -7.10;
  m->data[2][0] = 5.13; m->data[2][1] = 6.62; m->data[2][2] = -5.66 ; m->data[2][3] = 0.87;
  m->data[3][0] = -4.75; m->data[3][1] = 8.52; m->data[3][2] = 5.75; m->data[3][3] = 5.30;
  m->data[4][0] = 1.33; m->data[4][1] = 4.91; m->data[4][2] = -5.49 ; m->data[4][3] = -3.52;
  m->data[5][0] = -2.40; m->data[5][1] = -6.77; m->data[5][2] = 2.34 ; m->data[5][3] = 3.95;

  initMatrix(&U);
  initMatrix(&S);
  initMatrix(&V_T);

  SVDlapack(m, U, S, V_T);
  PrintMatrix(U);
  PrintMatrix(S);
  PrintMatrix(V_T);

  DelMatrix(&U);
  DelMatrix(&S);
  DelMatrix(&V_T);
  DelMatrix(&m);
}

void Test26_bis()
{
  size_t i, j, k;
  matrix *m, *U, *S, *V_T;
  matrix *Ur, *Sr, *V_Tr; /* Results */
  puts("Test 26 bis: Computing Singular Value Decomposition Lapack method");

  NewMatrix(&m, 2, 4);
  k = 1;
  for(i = 0; i < 2; i++){
    for(j = 0; j < 4; j++){
      m->data[i][j] = k;
      k++;
    }
  }

  initMatrix(&U);
  initMatrix(&S);
  initMatrix(&V_T);

  SVDlapack(m, U, S, V_T);

  NewMatrix(&Ur, 2, 2);
  NewMatrix(&Sr, 2, 2);
  NewMatrix(&V_Tr, 2, 4);
  Ur->data[0][0] = -0.376168; Ur->data[0][1] = -0.926551;
  Ur->data[1][0] = -0.926551; Ur->data[1][1] =0.376168;
  Sr->data[0][0] = 14.227407; Sr->data[1][1] = 1.257330;
  V_Tr->data[0][0] = -0.352062; V_Tr->data[0][1] = -0.443626; V_Tr->data[0][2] = -0.535190; V_Tr->data[0][3] = -0.626754;
  V_Tr->data[1][0] = 0.758981; V_Tr->data[1][1] = 0.321242; V_Tr->data[1][1] = -0.116498; V_Tr->data[1][1] = -0.554238;

  int ok = 1;
  for(int i = 0; i < 2; i++){
    for(int j = 0; j < 2; j++){
      if(FLOAT_EQ(U->data[i][j], Ur->data[i][j], 1e-2) &&
         FLOAT_EQ(S->data[i][j], Sr->data[i][j], 1e-2) &&
         FLOAT_EQ(V_T->data[i][j], V_Tr->data[i][j], 1e-2)){
         continue;
      }
      else{
        ok = 0;
        break;
      }
    }
  }

  if(ok == 1)
    puts("SVD Lapack: OK");
  else{
    printf("Problem on SVD Lapack decomposition\n");
    puts("Current answer");
    puts("U"); PrintMatrix(U);
    puts("S"); PrintMatrix(S);
    puts("V_T"); PrintMatrix(V_T);
    puts("Correct answer");
    puts("U"); PrintMatrix(Ur);
    puts("S"); PrintMatrix(Sr);
    puts("V_T"); PrintMatrix(V_T);
    puts("Error Test26_bis");
  }

  DelMatrix(&Ur);
  DelMatrix(&Sr);
  DelMatrix(&V_Tr);
  DelMatrix(&U);
  DelMatrix(&S);
  DelMatrix(&V_T);
  DelMatrix(&m);
}

void Test26()
{
  matrix *m, *U, *S, *V_T;
  matrix *Ur, *Sr, *V_Tr; /* Results */
  puts("Test 26: Computing Singular Value Decomposition Lapack method");

  NewMatrix(&m, 2, 2);
  m->data[0][0] = 7; m->data[0][1] = 2;
  m->data[1][0] = 3; m->data[1][1] = 4;

  initMatrix(&U);
  initMatrix(&S);
  initMatrix(&V_T);

  SVDlapack(m, U, S, V_T);

  NewMatrix(&Ur, 2, 2);
  NewMatrix(&Sr, 2, 2);
  NewMatrix(&V_Tr, 2, 2);
  Ur->data[0][0] = -0.847; Ur->data[0][1] = -0.532;
  Ur->data[1][0] = -0.532; Ur->data[1][1] = 0.847;
  Sr->data[0][0] = 8.438; Sr->data[1][1] = 0.f;
  Sr->data[0][0] = 8.438; Sr->data[0][1] = 0.f;
  Sr->data[1][0] = 0.f; Sr->data[1][1] = 2.607;
  V_Tr->data[0][0] = -0.892; V_Tr->data[0][1] = -0.453;
  V_Tr->data[1][0] = -0.453; V_Tr->data[1][1] = 0.892;

  int ok = 1;
  for(int i = 0; i < 2; i++){
    for(int j = 0; j < 2; j++){
      if(FLOAT_EQ(U->data[i][j], Ur->data[i][j], 1e-2) &&
         FLOAT_EQ(S->data[i][j], Sr->data[i][j], 1e-2) &&
         FLOAT_EQ(V_T->data[i][j], V_Tr->data[i][j], 1e-2)){
         continue;
      }
      else{
        ok = 0;
        break;
      }
    }
  }

  if(ok == 1)
    puts("SVD Lapack: OK");
  else{
    printf("Problem on SVD Lapack decomposition\n");
    puts("Current answer");
    puts("U"); PrintMatrix(U);
    puts("S"); PrintMatrix(S);
    puts("V_T"); PrintMatrix(V_T);
    puts("Correct answer");
    puts("U"); PrintMatrix(Ur);
    puts("S"); PrintMatrix(Sr);
    puts("V_T"); PrintMatrix(V_T);
    puts("Error Test26");
    abort();
  }

  DelMatrix(&Ur);
  DelMatrix(&Sr);
  DelMatrix(&V_Tr);
  DelMatrix(&U);
  DelMatrix(&S);
  DelMatrix(&V_T);
  DelMatrix(&m);
}

void Test25()
{
  matrix *m, *U, *S, *V_T;
  matrix *Ur, *Sr, *V_Tr; /* Results */
  puts("Test 25: Computing Singular Value Decomposition internal method");

  NewMatrix(&m, 2, 2);
  m->data[0][0] = 7; m->data[0][1] = 2;
  m->data[1][0] = 3; m->data[1][1] = 4;

  initMatrix(&U);
  initMatrix(&S);
  initMatrix(&V_T);

  SVDlapack(m, U, S, V_T);

  NewMatrix(&Ur, 2, 2);
  NewMatrix(&Sr, 2, 2);
  NewMatrix(&V_Tr, 2, 2);
  Ur->data[0][0] = -0.847; Ur->data[0][1] = -0.532;
  Ur->data[1][0] = -0.532; Ur->data[1][1] = 0.847;
  Sr->data[0][0] = 8.438; Sr->data[0][1] = 0.f;
  Sr->data[1][0] = 0.f; Sr->data[1][1] = 2.607;
  V_Tr->data[0][0] = -0.892; V_Tr->data[0][1] = -0.453;
  V_Tr->data[1][0] = -0.453; V_Tr->data[1][1] = 0.892;

  int ok = 1;
  for(int i = 0; i < 2; i++){
    for(int j = 0; j < 2; j++){
      if(FLOAT_EQ(U->data[i][j], Ur->data[i][j], 1e-2) &&
         FLOAT_EQ(S->data[i][j], Sr->data[i][j], 1e-2) &&
         FLOAT_EQ(V_T->data[i][j], V_Tr->data[i][j], 1e-2)){
         continue;
      }
      else{
        ok = 0;
        break;
      }
    }
  }

  if(ok == 1)
    puts("SVD Internal method: OK");
  else{
    printf("Problem on SVD decomposition\n");
    puts("Current answer");
    puts("U"); PrintMatrix(U);
    puts("S"); PrintMatrix(S);
    puts("V_T"); PrintMatrix(V_T);
    puts("Correct answer");
    puts("U"); PrintMatrix(Ur);
    puts("S"); PrintMatrix(Sr);
    puts("V_T"); PrintMatrix(V_T);
    puts("Error Test25");
    abort();
  }

  DelMatrix(&Ur);
  DelMatrix(&Sr);
  DelMatrix(&V_Tr);
  DelMatrix(&U);
  DelMatrix(&S);
  DelMatrix(&V_T);
  DelMatrix(&m);
}

void Test24()
{
  dvector *eval, *eval_answer;
  matrix *A, *evect, *evect_answer;
  NewMatrix(&A, 3, 3);
  A->data[0][0] = 7; A->data[0][1] = 2; A->data[0][2] = 5;
  A->data[1][0] = 3; A->data[1][1] = 4; A->data[1][2] = 8;
  A->data[2][0] = 7; A->data[2][1] = 3; A->data[2][2] = 9;

  puts("Test 24: Calculate Eigenvectors Eigenvalue for a random matrix");
  /*
  puts("Matrix");
  PrintMatrix(A);*/
  initDVector(&eval);
  initMatrix(&evect);
  EVectEval(A, eval, evect);

  NewDVector(&eval_answer, 3);
  NewMatrix(&evect_answer, 3, 3);
  eval_answer->data[0] = 16.393; eval_answer->data[1] = 2.425;  eval_answer->data[2] = 1.183;

  evect_answer->data[0][0] = -0.630 ; evect_answer->data[0][1] = -0.853; evect_answer->data[0][2] = -0.700;
  evect_answer->data[1][0] = -0.277; evect_answer->data[1][1] = 0.116; evect_answer->data[1][2] = -0.225;
  evect_answer->data[2][0] = -0.726; evect_answer->data[2][1] = 0.508; evect_answer->data[2][2] = 0.678;

  int ok = 1;
  for(int i = 0; i < 3; i++){
    if(FLOAT_EQ(eval_answer->data[i], eval->data[i], 1e-2)){
      for(int j = 0; j < 3; j++){
        if(FLOAT_EQ(evect_answer->data[i][j], evect->data[i][j], 1e-2)){
          continue;
        }
        else{
          ok = 0;
          break;
        }
      }
    }
    else{
      ok = 0;
      break;
    }
  }

  if(ok == 1)
    puts("Eigenvectors/eigenvalues decomposition method: OK");
  else{
    printf("Problem on eigenvectors/eigenvalues decomposition\n");
    puts("Current answer");
    puts("eigenvalues"); PrintDVector(eval);
    puts("eigenvectors"); PrintMatrix(evect);
    puts("Correct answer");
    puts("eigenvalues"); PrintDVector(eval_answer);
    puts("eigenvectors"); PrintMatrix(evect_answer);
    puts("Error Test24");
    abort();
  }

  DelMatrix(&evect_answer);
  DelDVector(&eval_answer);
  /*puts("Eigenvalue"); PrintDVector(eval);
  puts("Eigenvectors"); PrintMatrix(evect);*/

  DelDVector(&eval);
  DelMatrix(&evect);
  DelMatrix(&A);
}


void Test23()
{
  dvector *a, *r;
  matrix *b;
  NewDVector(&a, 3);
  NewMatrix(&b, 3, 3);
  a->data[0] = 1; a->data[1] = 2; a->data[2] = 3;

  b->data[0][0] = 1; b->data[0][1] = 0; b->data[0][2] = 0;
  b->data[1][0] = 1; b->data[1][1] = 2; b->data[1][2] = 1;
  b->data[2][0] = 0; b->data[2][1] = 1; b->data[2][2] = 3;


  puts("Test 23: Compute vector T matrix Division");
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

void Test22()
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

  puts("Test 22: Calculate Eigenvectors Eigenvalue for a sparse matrix");

  puts("Matrix");
  PrintMatrix(A);
  initDVector(&eval);
  initMatrix(&evect);
  EVectEval(A, eval, evect);

  puts("Eigenvalue"); PrintDVector(eval);
  puts("Eigenvectors"); PrintMatrix(evect);

  DelDVector(&eval);
  DelMatrix(&evect);
  DelMatrix(&A);
}


void Test21()
{
  matrix *m, *U, *S, *V_T, *US, *mbis;
  puts("Test 21: Computing Singular Value Decomposition");

  NewMatrix(&m, 2, 2);
  m->data[0][0] = 2; m->data[0][1] = 4;
  m->data[1][0] = 1; m->data[1][1] = 3;

  puts("INIT MATRIX");
  PrintMatrix(m);

  initMatrix(&U);
  initMatrix(&S);
  initMatrix(&V_T);

  SVD(m, U, S, V_T);

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
  QRDecomposition(m, Q, R);

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

  QRDecomposition(m, Q, R);

  puts("Q"); PrintMatrix(Q);
  puts("R"); PrintMatrix(R);

  puts("Original Matrix");
  PrintMatrix(m);
  puts("Recomputed Matrix");
  MatrixDotProduct(Q, R, mbis);
  PrintMatrix(mbis);

  initDVector(&eval);
  initMatrix(&evect);
  EVectEval(m, eval, evect);

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
  size_t i, j;
  matrix *m, *covm, *res;
  puts("Test Calculation Covariance Matrix");

  NewMatrix(&m, 5, 2);
  m->data[0][0] = 1; m->data[0][1] = 2;
  m->data[1][0] = 2; m->data[1][1] = 3;
  m->data[2][0] = 3; m->data[2][1] = 3;
  m->data[3][0] = 4; m->data[3][1] = 5;
  m->data[4][0] = 5; m->data[4][1] = 5;

  NewMatrix(&res, 2, 2);
  res->data[0][0] = 2.5; res->data[0][1] = 2.;
  res->data[1][0] = 2.; res->data[1][1] = 1.8;

  initMatrix(&covm);

  MatrixCovariance(m, covm);

  puts("Matrix to for covariance matrix");
  PrintMatrix(m);

  puts("Covariance Matrix");
  PrintMatrix(covm);

  for(i = 0; i < res->row; i++){
    for(j = 0; j < res->row; j++){
      if(!FLOAT_EQ(res->data[i][j], covm->data[i][j], 1e-1)){
        puts("Error Test18 - Calculation covariance matrix");
        abort();
      }
      else{
        continue;
      }
    }
  }

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
      srand(i+j);
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
      srand(i+j);
      setMatrixValue(m1, i, j, randDouble(1, 24));
    }
  }

  initDVector(&d1);
  MatrixColAverage(m1, d1);


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
  MatrixColAverage(m2, d2);

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
      srand(i+j);
      setMatrixValue(m, i, j, randDouble(-19, 24));
    }
  }

  initMatrix(&covm);

  MatrixCovariance(m, covm);

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
  MatrixCovariance(m, covm);

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
      srand(i+j);
      setMatrixValue(m1, i, j, randDouble(-19, 24));
    }
  }

  initMatrix(&m2);

  MatrixInversion(m1, m2);

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

  MatrixInversion(m1, m2);

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

  MatrixLUInversion(m1, m2);

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

  ResizeMatrix(m, 10, 4);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      setMatrixValue(m, i, j, i+j);
    }
  }

  puts("Initial Matrix");
  PrintMatrix(m);

  puts("Resizin Matrix to 5 x 2");

  ResizeMatrix(m, 5, 2);

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

  MatrixColAverage(m2, colmean);

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
  size_t i, j, rowsize = 2500, colsize = 500;
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
  srand(12);
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

  srand(12);

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

  MatrixAppendCol(m, col);

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

  srand(12);

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


  MatrixAppendRow(m, row);

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
  Test1();
  Test2();
  Test3();
  Test4();
  Test5();
  Test5_bis(); // SLOW
  Test6();
  Test7();
  Test8();
  Test9();
  Test10();
  Test11();
  Test12();
  Test13();
  Test14();
  Test15();
  Test16();
  Test17();
  Test18();
  Test19();
  Test20();
  Test21();
  Test22();
  Test23();
  Test24();
  Test25();
  Test26();
  Test26_bis();
  // Test26_tris();
  Test27(); // SLOW
  Test28();
  Test29();
  Test30();
  Test31();
  Test32();
  Test33();
  Test34();
  Test35();
  Test36();
  Test37();
  Test38();
  Test39();
  Test40();
  Test41();
  Test42();
  Test43();
  Test44();
  Test45();
  Test46();
  Test47();
  Test48();
  Test49();
  Test50();
  Test51();
  Test52();
  Test53();
  Test54();
  Test55();
  Test56();
  Test57();
  Test58();
  Test59();
  Test60();
  Test61();
  Test62();
  return 0;
}

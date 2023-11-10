
/* testmetrics.c
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
#include <math.h>
#include "metricspace.h"
#include "numeric.h"
#include "matrix.h"
#include "statistic.h"

void test13()
{
  puts("Test13: ManhattanDistanceCondensed");
  size_t i;
  size_t j;
  size_t indx;
  matrix *m;
  dvector *dist_condensed;
  matrix *dist_st;

  NewMatrix(&m, 4, 3);
  m->data[0][0] = 1; m->data[0][1] = 2; m->data[0][2] = 3;
  m->data[1][0] = 4; m->data[1][1] = 5; m->data[1][2] = 6;
  m->data[2][0] = 7; m->data[2][1] = 8; m->data[2][2] = 11;
  m->data[3][0] = 9; m->data[3][1] = 10; m->data[3][2] = 12;

  initDVector(&dist_condensed);
  ManhattanDistanceCondensed(m, dist_condensed, 2);

  initMatrix(&dist_st);
  ManhattanDistance_ST(m, m, dist_st);

  for(i = 0; i < dist_st->row; i++){
    double val;
    for(j = 0; j < dist_st->col; j++){
      if(i == j){
        val = 0.f;
      }
      else{
        indx = square_to_condensed_index(i, j, m->row);
        val = dist_condensed->data[indx];
      }

      if(FLOAT_EQ(val, dist_st->data[i][j], EPSILON)){
        continue;
      }
      else{
        puts("Test 15 Failed!");
        abort();
      }
    }
  }
  puts("ManhattanDistanceCondensed: OK");
  DelDVector(&dist_condensed);
  DelMatrix(&dist_st);
  DelMatrix(&m);
}

void test12()
{
  puts("Test12: SquaredEuclideanDistanceCondensed");
  size_t i;
  size_t j;
  size_t indx;
  matrix *m;
  dvector *dist_condensed;
  matrix *dist_st;
  matrix *dist_mt;

  NewMatrix(&m, 4, 3);
  m->data[0][0] = 1; m->data[0][1] = 2; m->data[0][2] = 3;
  m->data[1][0] = 4; m->data[1][1] = 5; m->data[1][2] = 6;
  m->data[2][0] = 7; m->data[2][1] = 8; m->data[2][2] = 11;
  m->data[3][0] = 9; m->data[3][1] = 10; m->data[3][2] = 12;

  initMatrix(& dist_mt);
  CalculateDistance(m, m, dist_mt, 4, SQUARE_EUCLIDEAN);

  initDVector(&dist_condensed);
  SquaredEuclideanDistanceCondensed(m, dist_condensed, 2);

  initMatrix(&dist_st);
  SquaredEuclideanDistance_ST(m, m, dist_st);

  for(i = 0; i < dist_st->row; i++){
    double val;
    for(j = 0; j < dist_st->col; j++){
      if(i == j){
        val = 0.f;
      }
      else{
        indx = square_to_condensed_index(i, j, m->row);
        val = dist_condensed->data[indx];
      }

      if(FLOAT_EQ(val, dist_st->data[i][j], EPSILON) && 
        FLOAT_EQ(val, dist_mt->data[i][j], EPSILON)){
        continue;
      }
      else{
        puts("Test 14 Failed!");
        abort();
      }
    }
  }
  puts("SquaredEuclideanDistanceCondensed: OK");
  DelDVector(&dist_condensed);
  DelMatrix(&dist_st);
  DelMatrix(&dist_mt);
  DelMatrix(&m);
}

void test11()
{
  puts("Test11: EuclideanDistanceCondensed");
  size_t i;
  size_t j;
  size_t indx;
  matrix *m;
  dvector *dist_condensed;
  matrix *dist_st;

  NewMatrix(&m, 4, 3);
  m->data[0][0] = 1; m->data[0][1] = 2; m->data[0][2] = 3;
  m->data[1][0] = 4; m->data[1][1] = 5; m->data[1][2] = 6;
  m->data[2][0] = 7; m->data[2][1] = 8; m->data[2][2] = 11;
  m->data[3][0] = 9; m->data[3][1] = 10; m->data[3][2] = 12;

  initDVector(&dist_condensed);
  EuclideanDistanceCondensed(m, dist_condensed, 2);

  initMatrix(&dist_st);
  EuclideanDistance_ST(m, m, dist_st);

  for(i = 0; i < dist_st->row; i++){
    double val;
    for(j = 0; j < dist_st->col; j++){
      if(i == j){
        val = 0.f;
      }
      else{
        indx = square_to_condensed_index(i, j, m->row);
        val = dist_condensed->data[indx];
      }

      if(FLOAT_EQ(val, dist_st->data[i][j], EPSILON)){
        continue;
      }
      else{
        puts("Test 13 Failed!");
        abort();
      }
    }
  }
  puts("EuclideanDistanceCondensed: OK");
  DelDVector(&dist_condensed);
  DelMatrix(&dist_st);
  DelMatrix(&m);
}

void test10()
{
  puts("Test10: ManhattanDistance");
  size_t i;
  size_t j;
  matrix *m1;
  matrix *m2;
  matrix *dist_mt;
  matrix *dist_st;

  NewMatrix(&m1, 4, 3);
  m1->data[0][0] = 1; m1->data[0][1] = 2; m1->data[0][2] = 3;
  m1->data[1][0] = 4; m1->data[1][1] = 5; m1->data[1][2] = 6;
  m1->data[2][0] = 7; m1->data[2][1] = 8; m1->data[2][2] = 11;
  m1->data[3][0] = 9; m1->data[3][1] = 10; m1->data[3][2] = 12;

  NewMatrix(&m2, 6, 3);
  m2->data[0][0] = 1; m2->data[0][1] = 2; m2->data[0][2] = 3;
  m2->data[1][0] = 1; m2->data[1][1] = 8; m2->data[1][2] = 1;
  m2->data[2][0] = 4; m2->data[2][1] = 0; m2->data[2][2] = 7;
  m2->data[3][0] = 5; m2->data[3][1] = 1; m2->data[3][2] = 6;
  m2->data[4][0] = 3; m2->data[4][1] = 3; m2->data[4][2] = 4;
  m2->data[5][0] = 5; m2->data[5][1] = 2; m2->data[5][2] = 7;


  initMatrix(&dist_mt);
  CalculateDistance(m2, m1, dist_mt, 2, MANHATTAN);

  initMatrix(&dist_st);
  ManhattanDistance_ST(m2, m1, dist_st);

  for(i = 0; i < dist_st->row; i++){
    for(j = 0; j < dist_st->col; j++){
      if(FLOAT_EQ(dist_mt->data[i][j], dist_st->data[i][j], EPSILON)){
        continue;
      }
      else{
        puts("Test 12 Failed!");
        abort();
      }
    }
  }
  puts("ManhattanDistance: OK");
  DelMatrix(&dist_mt);
  DelMatrix(&dist_st);
  DelMatrix(&m1);
  DelMatrix(&m2);
}

void test9()
{
  puts("Test9: CosineDistance");
  size_t i;
  size_t j;
  matrix *m1;
  matrix *m2;
  matrix *dist_mt;
  matrix *dist_st;

  NewMatrix(&m1, 4, 3);
  m1->data[0][0] = 1; m1->data[0][1] = 2; m1->data[0][2] = 3;
  m1->data[1][0] = 4; m1->data[1][1] = 5; m1->data[1][2] = 6;
  m1->data[2][0] = 7; m1->data[2][1] = 8; m1->data[2][2] = 11;
  m1->data[3][0] = 9; m1->data[3][1] = 10; m1->data[3][2] = 12;

  NewMatrix(&m2, 6, 3);
  m2->data[0][0] = 1; m2->data[0][1] = 2; m2->data[0][2] = 3;
  m2->data[1][0] = 1; m2->data[1][1] = 8; m2->data[1][2] = 1;
  m2->data[2][0] = 4; m2->data[2][1] = 0; m2->data[2][2] = 7;
  m2->data[3][0] = 5; m2->data[3][1] = 1; m2->data[3][2] = 6;
  m2->data[4][0] = 3; m2->data[4][1] = 3; m2->data[4][2] = 4;
  m2->data[5][0] = 5; m2->data[5][1] = 2; m2->data[5][2] = 7;


  initMatrix(&dist_mt);
  CalculateDistance(m2, m1, dist_mt, 2, COSINE);

  initMatrix(&dist_st);
  CosineDistance_ST(m2, m1, dist_st);

  for(i = 0; i < dist_st->row; i++){
    for(j = 0; j < dist_st->col; j++){
      if(FLOAT_EQ(dist_mt->data[i][j], dist_st->data[i][j], EPSILON)){
        continue;
      }
      else{
        puts("Test11 Failed!");
        abort();
      }
    }
  }
  puts("CosineDistance: OK");
  DelMatrix(&dist_mt);
  DelMatrix(&dist_st);
  DelMatrix(&m1);
  DelMatrix(&m2);
}

void test8()
{
  puts("Test8: SquaredEuclideanDistance");
  size_t i;
  size_t j;
  matrix *m1;
  matrix *m2;
  matrix *dist_mt;
  matrix *dist_st;

  NewMatrix(&m1, 4, 3);
  m1->data[0][0] = 1; m1->data[0][1] = 2; m1->data[0][2] = 3;
  m1->data[1][0] = 4; m1->data[1][1] = 5; m1->data[1][2] = 6;
  m1->data[2][0] = 7; m1->data[2][1] = 8; m1->data[2][2] = 11;
  m1->data[3][0] = 9; m1->data[3][1] = 10; m1->data[3][2] = 12;

  NewMatrix(&m2, 6, 3);
  m2->data[0][0] = 1; m2->data[0][1] = 2; m2->data[0][2] = 3;
  m2->data[1][0] = 1; m2->data[1][1] = 8; m2->data[1][2] = 1;
  m2->data[2][0] = 4; m2->data[2][1] = 0; m2->data[2][2] = 7;
  m2->data[3][0] = 5; m2->data[3][1] = 1; m2->data[3][2] = 6;
  m2->data[4][0] = 3; m2->data[4][1] = 3; m2->data[4][2] = 4;
  m2->data[5][0] = 5; m2->data[5][1] = 2; m2->data[5][2] = 7;


  initMatrix(&dist_mt);
  CalculateDistance(m2, m1, dist_mt, 2, SQUARE_EUCLIDEAN);

  initMatrix(&dist_st);
  SquaredEuclideanDistance_ST(m2, m1, dist_st);

  for(i = 0; i < dist_st->row; i++){
    for(j = 0; j < dist_st->col; j++){
      if(FLOAT_EQ(dist_mt->data[i][j], dist_st->data[i][j], EPSILON)){
        continue;
      }
      else{
        puts("Test10 Failed!");
        abort();
      }
    }
  }
  puts("SquaredEuclideanDistance: OK");
  DelMatrix(&dist_mt);
  DelMatrix(&dist_st);
  DelMatrix(&m1);
  DelMatrix(&m2);
}

void test7()
{
  puts("Test7: MatrixMatrixDistance");
  /* Euclidean distance between two matrix test */
  matrix *m1;
  matrix *m2;
  double mxmxdist;

  NewMatrix(&m1, 4, 3);
  m1->data[0][0] = 1; m1->data[0][1] = 2; m1->data[0][2] = 3;
  m1->data[1][0] = 4; m1->data[1][1] = 5; m1->data[1][2] = 6;
  m1->data[2][0] = 7; m1->data[2][1] = 8; m1->data[2][2] = 11;
  m1->data[3][0] = 9; m1->data[3][1] = 10; m1->data[3][2] = 12;

  NewMatrix(&m2, 6, 3);
  m2->data[0][0] = 1; m2->data[0][1] = 2; m2->data[0][2] = 3;
  m2->data[1][0] = 1; m2->data[1][1] = 8; m2->data[1][2] = 1;
  m2->data[2][0] = 4; m2->data[2][1] = 0; m2->data[2][2] = 7;
  m2->data[3][0] = 5; m2->data[3][1] = 1; m2->data[3][2] = 6;
  m2->data[4][0] = 3; m2->data[4][1] = 3; m2->data[4][2] = 4;
  m2->data[5][0] = 5; m2->data[5][1] = 2; m2->data[5][2] = 7;

  mxmxdist = MatrixMatrixDistance(m1, m2);
  printf("%f\n", mxmxdist);
  if(!FLOAT_EQ(mxmxdist, 4.980518, 1e-4))
    abort();
  puts("MatrixMatrixDistance: OK");
  DelMatrix(&m1);
  DelMatrix(&m2);
}
void test6()
{
  puts("Test6: CosineDistance");
  size_t i, j, maxrow, maxcol;
  matrix *m, *dist;

  maxrow = 600;
  maxcol = 10;
  srand_(maxrow);
  NewMatrix(&m, maxrow, maxcol);

  srand(maxrow+maxcol);
  for(i = 0; i < maxrow; i++){
    for(j = 0; j < maxcol; j++){
      setMatrixValue(m, i, j, rand() % 100);
    }
  }

  initMatrix(&dist);
  CalculateDistance(m, m, dist, 4, COSINE);

  puts("CosineDistance: OK");
  DelMatrix(&m);
  DelMatrix(&dist);
}

void test5()
{
  puts("Test5: SquaredEuclideanDistance");
  size_t i, j, maxrow, maxcol;
  matrix *m, *dist;

  maxrow = 600;
  maxcol = 10;
  srand_(maxrow);
  NewMatrix(&m, maxrow, maxcol);

  srand(maxrow+maxcol);
  for(i = 0; i < maxrow; i++){
    for(j = 0; j < maxcol; j++){
      setMatrixValue(m, i, j, rand() % 100);
    }
  }

  initMatrix(&dist);
  CalculateDistance(m, m, dist, 4, SQUARE_EUCLIDEAN);
  puts("SquaredEuclideanDistance: OK");
  DelMatrix(&m);
  DelMatrix(&dist);
}

void test4()
{
  puts("Test4: CosineDistanceCondensed");
  size_t i, j, indx;
  matrix *m;
  dvector *dist;
  NewMatrix(&m, 4, 4);
  m->data[0][0] = 0; m->data[0][1] = 1;
  m->data[1][0] = 1; m->data[1][1] = 1;
  m->data[2][0] = 3; m->data[2][1] = 5;
  m->data[3][0] = 15; m->data[3][1] = 5;
  initDVector(&dist);
  CosineDistanceCondensed(m, dist, 2);
  for(i = 0; i < 4; i++){
    for(j = 0; j < 4; j++){
      if(i == j){
        printf(" 0.00");
      }
      else{
        indx = square_to_condensed_index(i, j, m->row);
        printf(" %.2f", (dist->data[indx]));
      }
    }
    printf("\n");
  }
  puts("CosineDistanceCondensed: OK");
  DelDVector(&dist);
  DelMatrix(&m);
}

void test3()
{
  puts("Test3: CovarianceDistanceMap.");
  matrix *mi, *mo;

  NewMatrix(&mi, 100, 2);
  NewMatrix(&mo, 100, 2);

  mi->data[0][0] = -4.503286; mi->data[0][1] = -5.138264;
  mi->data[1][0] = -4.352311; mi->data[1][1] = -3.476970;
  mi->data[2][0] = -5.234153; mi->data[2][1] = -5.234137;
  mi->data[3][0] = -3.420787; mi->data[3][1] = -4.232565;
  mi->data[4][0] = -5.469474; mi->data[4][1] = -4.457440;
  mi->data[5][0] = -5.463418; mi->data[5][1] = -5.465730;
  mi->data[6][0] = -4.758038; mi->data[6][1] = -6.913280;
  mi->data[7][0] = -6.724918; mi->data[7][1] = -5.562288;
  mi->data[8][0] = -6.012831; mi->data[8][1] = -4.685753;
  mi->data[9][0] = -5.908024; mi->data[9][1] = -6.412304;
  mi->data[10][0] = -3.534351; mi->data[10][1] = -5.225776;
  mi->data[11][0] = -4.932472; mi->data[11][1] = -6.424748;
  mi->data[12][0] = -5.544383; mi->data[12][1] = -4.889077;
  mi->data[13][0] = -6.150994; mi->data[13][1] = -4.624302;
  mi->data[14][0] = -5.600639; mi->data[14][1] = -5.291694;
  mi->data[15][0] = -5.601707; mi->data[15][1] = -3.147722;
  mi->data[16][0] = -5.013497; mi->data[16][1] = -6.057711;
  mi->data[17][0] = -4.177455; mi->data[17][1] = -6.220844;
  mi->data[18][0] = -4.791136; mi->data[18][1] = -6.959670;
  mi->data[19][0] = -6.328186; mi->data[19][1] = -4.803139;
  mi->data[20][0] = -4.261533; mi->data[20][1] = -4.828632;
  mi->data[21][0] = -5.115648; mi->data[21][1] = -5.301104;
  mi->data[22][0] = -6.478522; mi->data[22][1] = -5.719844;
  mi->data[23][0] = -5.460639; mi->data[23][1] = -3.942878;
  mi->data[24][0] = -4.656382; mi->data[24][1] = -6.763040;
  mi->data[25][0] = -4.675916; mi->data[25][1] = -5.385082;
  mi->data[26][0] = -5.676922; mi->data[26][1] = -4.388324;
  mi->data[27][0] = -3.969000; mi->data[27][1] = -4.068720;
  mi->data[28][0] = -5.839218; mi->data[28][1] = -5.309212;
  mi->data[29][0] = -4.668737; mi->data[29][1] = -4.024455;
  mi->data[30][0] = -5.479174; mi->data[30][1] = -5.185659;
  mi->data[31][0] = -6.106335; mi->data[31][1] = -6.196207;
  mi->data[32][0] = -4.187474; mi->data[32][1] = -3.643760;
  mi->data[33][0] = -5.072010; mi->data[33][1] = -3.996467;
  mi->data[34][0] = -4.638364; mi->data[34][1] = -5.645120;
  mi->data[35][0] = -4.638604; mi->data[35][1] = -3.461963;
  mi->data[36][0] = -5.035826; mi->data[36][1] = -3.435356;
  mi->data[37][0] = -7.619745; mi->data[37][1] = -4.178097;
  mi->data[38][0] = -4.912953; mi->data[38][1] = -5.299007;
  mi->data[39][0] = -4.908239; mi->data[39][1] = -6.987569;
  mi->data[40][0] = -5.219672; mi->data[40][1] = -4.642887;
  mi->data[41][0] = -3.522106; mi->data[41][1] = -5.518270;
  mi->data[42][0] = -5.808494; mi->data[42][1] = -5.501757;
  mi->data[43][0] = -4.084598; mi->data[43][1] = -4.671249;
  mi->data[44][0] = -5.529760; mi->data[44][1] = -4.486733;
  mi->data[45][0] = -4.902922; mi->data[45][1] = -4.031355;
  mi->data[46][0] = -5.702053; mi->data[46][1] = -5.327662;
  mi->data[47][0] = -5.392108; mi->data[47][1] = -6.463515;
  mi->data[48][0] = -4.703880; mi->data[48][1] = -4.738945;
  mi->data[49][0] = -4.994887; mi->data[49][1] = -5.234587;
  mi->data[50][0] = 1.584629; mi->data[50][1] = 4.579355;
  mi->data[51][0] = 2.657285; mi->data[51][1] = 4.197723;
  mi->data[52][0] = 2.838714; mi->data[52][1] = 5.404051;
  mi->data[53][0] = 4.886186; mi->data[53][1] = 5.174578;
  mi->data[54][0] = 3.257550; mi->data[54][1] = 4.925554;
  mi->data[55][0] = 1.081229; mi->data[55][1] = 4.973486;
  mi->data[56][0] = 3.060230; mi->data[56][1] = 7.463242;
  mi->data[57][0] = 2.807639; mi->data[57][1] = 5.301547;
  mi->data[58][0] = 2.965288; mi->data[58][1] = 3.831322;
  mi->data[59][0] = 4.142823; mi->data[59][1] = 5.751933;
  mi->data[60][0] = 3.791032; mi->data[60][1] = 4.090613;
  mi->data[61][0] = 4.402794; mi->data[61][1] = 3.598149;
  mi->data[62][0] = 3.586857; mi->data[62][1] = 7.190456;
  mi->data[63][0] = 2.009464; mi->data[63][1] = 4.433702;
  mi->data[64][0] = 3.099651; mi->data[64][1] = 4.496524;
  mi->data[65][0] = 1.449337; mi->data[65][1] = 5.068563;
  mi->data[66][0] = 1.937696; mi->data[66][1] = 5.473592;
  mi->data[67][0] = 2.080576; mi->data[67][1] = 6.549934;
  mi->data[68][0] = 2.216747; mi->data[68][1] = 4.677938;
  mi->data[69][0] = 3.813517; mi->data[69][1] = 3.769136;
  mi->data[70][0] = 3.227460; mi->data[70][1] = 6.307143;
  mi->data[71][0] = 1.392517; mi->data[71][1] = 5.184634;
  mi->data[72][0] = 3.259883; mi->data[72][1] = 5.781823;
  mi->data[73][0] = 1.763049; mi->data[73][1] = 3.679543;
  mi->data[74][0] = 3.521942; mi->data[74][1] = 5.296985;
  mi->data[75][0] = 3.250493; mi->data[75][1] = 5.346448;
  mi->data[76][0] = 2.319975; mi->data[76][1] = 5.232254;
  mi->data[77][0] = 3.293072; mi->data[77][1] = 4.285649;
  mi->data[78][0] = 4.865775; mi->data[78][1] = 5.473833;
  mi->data[79][0] = 1.808697; mi->data[79][1] = 5.656554;
  mi->data[80][0] = 2.025318; mi->data[80][1] = 5.787085;
  mi->data[81][0] = 4.158596; mi->data[81][1] = 4.179318;
  mi->data[82][0] = 3.963376; mi->data[82][1] = 5.412781;
  mi->data[83][0] = 3.822060; mi->data[83][1] = 6.896793;
  mi->data[84][0] = 2.754612; mi->data[84][1] = 4.246264;
  mi->data[85][0] = 2.110486; mi->data[85][1] = 4.184190;
  mi->data[86][0] = 2.922898; mi->data[86][1] = 5.341152;
  mi->data[87][0] = 3.276691; mi->data[87][1] = 5.827183;
  mi->data[88][0] = 3.013002; mi->data[88][1] = 6.453534;
  mi->data[89][0] = 2.735343; mi->data[89][1] = 7.720169;
  mi->data[90][0] = 3.625667; mi->data[90][1] = 4.142842;
  mi->data[91][0] = 1.929108; mi->data[91][1] = 5.482472;
  mi->data[92][0] = 2.776537; mi->data[92][1] = 5.714000;
  mi->data[93][0] = 3.473238; mi->data[93][1] = 4.927171;
  mi->data[94][0] = 2.153206; mi->data[94][1] = 3.485153;
  mi->data[95][0] = 2.553485; mi->data[95][1] = 5.856399;
  mi->data[96][0] = 3.214094; mi->data[96][1] = 3.754261;
  mi->data[97][0] = 3.173181; mi->data[97][1] = 5.385317;
  mi->data[98][0] = 2.116143; mi->data[98][1] = 5.153725;
  mi->data[99][0] = 3.058209; mi->data[99][1] = 3.857030;
  CovarianceDistanceMap(mi, mo);
  DelMatrix(&mi);
  DelMatrix(&mo);
}


void test2()
{
  puts("Test2: EuclideanDistance");
  size_t i, j, maxrow, maxcol;
  matrix *m, *dist;

  maxrow = 100;
  maxcol = 50;
  srand_(maxrow);
  NewMatrix(&m, maxrow, maxcol);

  srand(maxrow+maxcol);
  for(i = 0; i < maxrow; i++){
    for(j = 0; j < maxcol; j++){
      setMatrixValue(m, i, j, rand() % 100);
    }
  }

  initMatrix(&dist);
  CalculateDistance(m, m, dist, 4, EUCLIDEAN);
  DelMatrix(&m);
  DelMatrix(&dist);
}

void test1()
{
  puts("Test1: EuclideanDistance");
  size_t i;
  size_t j;
  matrix *m1;
  matrix *m2;
  matrix *dist_mt;
  matrix *dist_st;

  NewMatrix(&m1, 4, 3);
  m1->data[0][0] = 1; m1->data[0][1] = 2; m1->data[0][2] = 3;
  m1->data[1][0] = 4; m1->data[1][1] = 5; m1->data[1][2] = 6;
  m1->data[2][0] = 7; m1->data[2][1] = 8; m1->data[2][2] = 11;
  m1->data[3][0] = 9; m1->data[3][1] = 10; m1->data[3][2] = 12;

  NewMatrix(&m2, 6, 3);
  m2->data[0][0] = 1; m2->data[0][1] = 2; m2->data[0][2] = 3;
  m2->data[1][0] = 1; m2->data[1][1] = 8; m2->data[1][2] = 1;
  m2->data[2][0] = 4; m2->data[2][1] = 0; m2->data[2][2] = 7;
  m2->data[3][0] = 5; m2->data[3][1] = 1; m2->data[3][2] = 6;
  m2->data[4][0] = 3; m2->data[4][1] = 3; m2->data[4][2] = 4;
  m2->data[5][0] = 5; m2->data[5][1] = 2; m2->data[5][2] = 7;


  initMatrix(&dist_mt);
  CalculateDistance(m2, m1, dist_mt, 4, EUCLIDEAN);

  initMatrix(&dist_st);
  EuclideanDistance_ST(m2, m1, dist_st);

  for(i = 0; i < dist_st->row; i++){
    for(j = 0; j < dist_st->col; j++){
      if(FLOAT_EQ(dist_mt->data[i][j], dist_st->data[i][j], EPSILON)){
        continue;
      }
      else{
        puts("Test1 Failed!");
        abort();
      }
    }
  }
  puts("EuclideanDistance: OK");
  DelMatrix(&dist_mt);
  DelMatrix(&dist_st);
  DelMatrix(&m1);
  DelMatrix(&m2);
}

int main(void)
{
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
  test8();
  test9();
  test10();
  test11();
  test12();
  test13();
}

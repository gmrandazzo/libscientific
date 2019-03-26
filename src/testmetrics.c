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


void test6()
{
  size_t i, j, maxrow, maxcol;
  matrix *m;
  dvector *dist;

  maxrow = 60000;
  maxcol = 784;

  NewMatrix(&m, maxrow, maxcol);

  srand(maxrow+maxcol);
  for(i = 0; i < maxrow; i++){
    for(j = 0; j < maxcol; j++){
      setMatrixValue(m, i, j, rand() % 100);
    }
  }

  initDVector(&dist);
  MahalanobisDistance(m, NULL, NULL, &dist);

  DelMatrix(&m);
  DelDVector(&dist);
}

void test5()
{
  size_t i, j;
  matrix *m, *c, *edst;
  dvector *mdst;

  NewMatrix(&m, 100, 2);

  m->data[0][0] = -4.503286; m->data[0][1] = -5.138264;
  m->data[1][0] = -4.352311; m->data[1][1] = -3.476970;
  m->data[2][0] = -5.234153; m->data[2][1] = -5.234137;
  m->data[3][0] = -3.420787; m->data[3][1] = -4.232565;
  m->data[4][0] = -5.469474; m->data[4][1] = -4.457440;
  m->data[5][0] = -5.463418; m->data[5][1] = -5.465730;
  m->data[6][0] = -4.758038; m->data[6][1] = -6.913280;
  m->data[7][0] = -6.724918; m->data[7][1] = -5.562288;
  m->data[8][0] = -6.012831; m->data[8][1] = -4.685753;
  m->data[9][0] = -5.908024; m->data[9][1] = -6.412304;
  m->data[10][0] = -3.534351; m->data[10][1] = -5.225776;
  m->data[11][0] = -4.932472; m->data[11][1] = -6.424748;
  m->data[12][0] = -5.544383; m->data[12][1] = -4.889077;
  m->data[13][0] = -6.150994; m->data[13][1] = -4.624302;
  m->data[14][0] = -5.600639; m->data[14][1] = -5.291694;
  m->data[15][0] = -5.601707; m->data[15][1] = -3.147722;
  m->data[16][0] = -5.013497; m->data[16][1] = -6.057711;
  m->data[17][0] = -4.177455; m->data[17][1] = -6.220844;
  m->data[18][0] = -4.791136; m->data[18][1] = -6.959670;
  m->data[19][0] = -6.328186; m->data[19][1] = -4.803139;
  m->data[20][0] = -4.261533; m->data[20][1] = -4.828632;
  m->data[21][0] = -5.115648; m->data[21][1] = -5.301104;
  m->data[22][0] = -6.478522; m->data[22][1] = -5.719844;
  m->data[23][0] = -5.460639; m->data[23][1] = -3.942878;
  m->data[24][0] = -4.656382; m->data[24][1] = -6.763040;
  m->data[25][0] = -4.675916; m->data[25][1] = -5.385082;
  m->data[26][0] = -5.676922; m->data[26][1] = -4.388324;
  m->data[27][0] = -3.969000; m->data[27][1] = -4.068720;
  m->data[28][0] = -5.839218; m->data[28][1] = -5.309212;
  m->data[29][0] = -4.668737; m->data[29][1] = -4.024455;
  m->data[30][0] = -5.479174; m->data[30][1] = -5.185659;
  m->data[31][0] = -6.106335; m->data[31][1] = -6.196207;
  m->data[32][0] = -4.187474; m->data[32][1] = -3.643760;
  m->data[33][0] = -5.072010; m->data[33][1] = -3.996467;
  m->data[34][0] = -4.638364; m->data[34][1] = -5.645120;
  m->data[35][0] = -4.638604; m->data[35][1] = -3.461963;
  m->data[36][0] = -5.035826; m->data[36][1] = -3.435356;
  m->data[37][0] = -7.619745; m->data[37][1] = -4.178097;
  m->data[38][0] = -4.912953; m->data[38][1] = -5.299007;
  m->data[39][0] = -4.908239; m->data[39][1] = -6.987569;
  m->data[40][0] = -5.219672; m->data[40][1] = -4.642887;
  m->data[41][0] = -3.522106; m->data[41][1] = -5.518270;
  m->data[42][0] = -5.808494; m->data[42][1] = -5.501757;
  m->data[43][0] = -4.084598; m->data[43][1] = -4.671249;
  m->data[44][0] = -5.529760; m->data[44][1] = -4.486733;
  m->data[45][0] = -4.902922; m->data[45][1] = -4.031355;
  m->data[46][0] = -5.702053; m->data[46][1] = -5.327662;
  m->data[47][0] = -5.392108; m->data[47][1] = -6.463515;
  m->data[48][0] = -4.703880; m->data[48][1] = -4.738945;
  m->data[49][0] = -4.994887; m->data[49][1] = -5.234587;
  m->data[50][0] = 1.584629; m->data[50][1] = 4.579355;
  m->data[51][0] = 2.657285; m->data[51][1] = 4.197723;
  m->data[52][0] = 2.838714; m->data[52][1] = 5.404051;
  m->data[53][0] = 4.886186; m->data[53][1] = 5.174578;
  m->data[54][0] = 3.257550; m->data[54][1] = 4.925554;
  m->data[55][0] = 1.081229; m->data[55][1] = 4.973486;
  m->data[56][0] = 3.060230; m->data[56][1] = 7.463242;
  m->data[57][0] = 2.807639; m->data[57][1] = 5.301547;
  m->data[58][0] = 2.965288; m->data[58][1] = 3.831322;
  m->data[59][0] = 4.142823; m->data[59][1] = 5.751933;
  m->data[60][0] = 3.791032; m->data[60][1] = 4.090613;
  m->data[61][0] = 4.402794; m->data[61][1] = 3.598149;
  m->data[62][0] = 3.586857; m->data[62][1] = 7.190456;
  m->data[63][0] = 2.009464; m->data[63][1] = 4.433702;
  m->data[64][0] = 3.099651; m->data[64][1] = 4.496524;
  m->data[65][0] = 1.449337; m->data[65][1] = 5.068563;
  m->data[66][0] = 1.937696; m->data[66][1] = 5.473592;
  m->data[67][0] = 2.080576; m->data[67][1] = 6.549934;
  m->data[68][0] = 2.216747; m->data[68][1] = 4.677938;
  m->data[69][0] = 3.813517; m->data[69][1] = 3.769136;
  m->data[70][0] = 3.227460; m->data[70][1] = 6.307143;
  m->data[71][0] = 1.392517; m->data[71][1] = 5.184634;
  m->data[72][0] = 3.259883; m->data[72][1] = 5.781823;
  m->data[73][0] = 1.763049; m->data[73][1] = 3.679543;
  m->data[74][0] = 3.521942; m->data[74][1] = 5.296985;
  m->data[75][0] = 3.250493; m->data[75][1] = 5.346448;
  m->data[76][0] = 2.319975; m->data[76][1] = 5.232254;
  m->data[77][0] = 3.293072; m->data[77][1] = 4.285649;
  m->data[78][0] = 4.865775; m->data[78][1] = 5.473833;
  m->data[79][0] = 1.808697; m->data[79][1] = 5.656554;
  m->data[80][0] = 2.025318; m->data[80][1] = 5.787085;
  m->data[81][0] = 4.158596; m->data[81][1] = 4.179318;
  m->data[82][0] = 3.963376; m->data[82][1] = 5.412781;
  m->data[83][0] = 3.822060; m->data[83][1] = 6.896793;
  m->data[84][0] = 2.754612; m->data[84][1] = 4.246264;
  m->data[85][0] = 2.110486; m->data[85][1] = 4.184190;
  m->data[86][0] = 2.922898; m->data[86][1] = 5.341152;
  m->data[87][0] = 3.276691; m->data[87][1] = 5.827183;
  m->data[88][0] = 3.013002; m->data[88][1] = 6.453534;
  m->data[89][0] = 2.735343; m->data[89][1] = 7.720169;
  m->data[90][0] = 3.625667; m->data[90][1] = 4.142842;
  m->data[91][0] = 1.929108; m->data[91][1] = 5.482472;
  m->data[92][0] = 2.776537; m->data[92][1] = 5.714000;
  m->data[93][0] = 3.473238; m->data[93][1] = 4.927171;
  m->data[94][0] = 2.153206; m->data[94][1] = 3.485153;
  m->data[95][0] = 2.553485; m->data[95][1] = 5.856399;
  m->data[96][0] = 3.214094; m->data[96][1] = 3.754261;
  m->data[97][0] = 3.173181; m->data[97][1] = 5.385317;
  m->data[98][0] = 2.116143; m->data[98][1] = 5.153725;
  m->data[99][0] = 3.058209; m->data[99][1] = 3.857030;

  initDVector(&mdst);
  dvector *mu;
  matrix *invcov;
  initDVector(&mu);
  initMatrix(&invcov);
  MahalanobisDistance(m, &invcov, &mu, &mdst);

  PrintMatrix(invcov);
  PrintDVector(mu);
  DelMatrix(&invcov);
  DelDVector(&mu);

  NewMatrix(&c, 1, 2);

  for(j = 0; j < m->col; j++){
    for(i = 0; i < m->row; i++){
      c->data[0][j] += m->data[i][j];
    }
    c->data[0][j] /= (double)m->row;
  }



  initMatrix(&edst);
  ManhattanDistance(m, c, &edst, 4);

  for(i = 0; i < m->row; i++){
    printf("%.4f %.4f\n", mdst->data[i], edst->data[0][i]);
  }

  DelMatrix(&edst);
  DelMatrix(&c);
  DelMatrix(&m);
  DelDVector(&mdst);
}

void test4()
{
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

  PrintMatrix(mi);

  CovarianceDistanceMap(mi, &mo);

  PrintMatrix(mo);

  DelMatrix(&mi);
  DelMatrix(&mo);
}


void test3()
{
  size_t i, j, maxrow, maxcol;
  matrix *m, *dist;

  maxrow = 10000;
  maxcol = 50;

  NewMatrix(&m, maxrow, maxcol);

  srand(maxrow+maxcol);
  for(i = 0; i < maxrow; i++){
    for(j = 0; j < maxcol; j++){
      setMatrixValue(m, i, j, rand() % 100);
    }
  }

  initMatrix(&dist);
  EuclideanDistance(m, m, &dist, 4);
  DelMatrix(&m);
  DelMatrix(&dist);
}

void test2()
{
  puts("Test3: ROC and Precision-Recall test.");
  /* ROC curve test */
  dvector *y_score, *y_true;
  matrix *roc, *pr;
  double auc, ap;
  NewDVector(&y_score, 13);
  NewDVector(&y_true, 13);

  y_score->data[0] = 0.1;
  y_score->data[1] = 0.2;
  y_score->data[2] = 0.3;
  y_score->data[3] = 0.4;
  y_score->data[4] = 0.5;
  y_score->data[5] = 0.6;
  y_score->data[6] = 0.7;
  y_score->data[7] = 0.8;
  y_score->data[8] = 0.9;
  y_score->data[9] = 1.0;
  y_score->data[10] = 2.0;
  y_score->data[11] = 3.0;
  y_score->data[12] = 3.4;


  y_true->data[0] = 0;
  y_true->data[1] = 0;
  y_true->data[2] = 0;
  y_true->data[3] = 0;
  y_true->data[4] = 0;
  y_true->data[5] = 1;
  y_true->data[6] = 1;
  y_true->data[7] = 1;
  y_true->data[8] = 1;
  y_true->data[9] = 0;
  y_true->data[10] = 1;
  y_true->data[11] = 1;
  y_true->data[12] = 1;

  initMatrix(&roc);
  ROC(y_true, y_score,  &roc, &auc);
  PrintMatrix(roc);
  if(FLOAT_EQ(auc, 0.904762, 1e-6)){
    printf("AUC OK!\n");
  }else{
    printf("AUC ERROR!\n");
  }
  printf("AUC: %f\n", auc);

  initMatrix(&pr);
  PrecisionRecall(y_true, y_score,  &pr, &ap);
  PrintMatrix(pr);

  if(FLOAT_EQ(ap, 0.900425, 1e-6)){
    printf("AVERAGE PRECISION-RECALL OK!\n");
  }else{
    printf("AVERAGE PRECISION-RECALL ERROR!\n");
  }
  printf("AVERAGE PRECISION-RECALL: %f\n", ap);

  DelMatrix(&pr);
  DelMatrix(&roc);
  DelDVector(&y_true);
  DelDVector(&y_score);
}

void test1()
{
  printf("Test1: Euclidean distance between two matrix.");
  /* Euclidean distance between two matrix test */
  matrix *m;
  matrix *dist;

  NewMatrix(&m, 4, 3);
  m->data[0][0] = 1; m->data[0][1] = 2; m->data[0][2] = 3;
  m->data[1][0] = 4; m->data[1][1] = 5; m->data[1][2] = 6;
  m->data[2][0] = 7; m->data[2][1] = 8; m->data[2][2] = 11;
  m->data[3][0] = 9; m->data[3][1] = 10; m->data[3][2] = 12;

  initMatrix(&dist);
  EuclideanDistance(m, m, &dist, 4);

  puts("Matrix");
  PrintMatrix(m);
  puts("Distance Matrix");
  PrintMatrix(dist);

  DelMatrix(&dist);
  DelMatrix(&m);
}

int main(void)
{
  /**/test1();
  test2();
  test3();
  test4();
  test5();
  test6();
}

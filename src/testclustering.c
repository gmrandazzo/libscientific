/* testclustering.c
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
#include <time.h>
#include "matrix.h"
#include "clustering.h"
#include "metricspace.h"

void test19()
{
  size_t i, j, maxrow, maxcol;
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;

  uivector *selections;

  NewMatrix(&m, 14, 2);


  maxrow = 1000;
  maxcol = 50;

  NewMatrix(&m, maxrow, maxcol);

  srand(maxrow+maxcol);
  for(i = 0; i < maxrow; i++){
    for(j = 0; j < maxcol; j++){
      setMatrixValue(m, i, j, rand() % 100);
    }
  }

  initUIVector(&selections);

  MDC(m, 0, 0, &selections, 8, &run);
  printf("Selected compounds %llu\n", selections->size);
  DelUIVector(&selections);
  DelMatrix(&m);
}

void test18()
{
  matrix *m, *centroids;
  int run = SIGSCIENTIFICRUN;
  uivector *clusters;

  size_t i, j, maxrow, maxcol;
  maxrow = 10000;
  maxcol = 20;

  NewMatrix(&m, maxrow, maxcol);

  srand(maxrow+maxcol);
  for(i = 0; i < maxrow; i++){
    for(j = 0; j < maxcol; j++){
      setMatrixValue(m, i, j, rand() % 100);
    }
  }

  initUIVector(&clusters);
  initMatrix(&centroids);
  puts("KMeans++ Clustering");

  KMeans(m, 1000, 0, &clusters, &centroids, 4, &run);


  DelUIVector(&clusters);
  DelMatrix(&centroids);
  DelMatrix(&m);
}

void test17()
{
  matrix *m, *centroids;
  int run = SIGSCIENTIFICRUN;
  uivector *clusters;

  size_t i, j, maxrow, maxcol;
  maxrow = 100;
  maxcol = 50;

  NewMatrix(&m, maxrow, maxcol);

  srand(maxrow+maxcol);
  for(i = 0; i < maxrow; i++){
    for(j = 0; j < maxcol; j++){
      setMatrixValue(m, i, j, rand() % 100);
    }
  }

  initUIVector(&clusters);
  initMatrix(&centroids);
  puts("Hierarchical Clustering");

  HierarchicalClustering(m, 3, &clusters, &centroids, NULL, 0, 4, &run);

  /*puts("clusters");
  PrintUIVector(clusters);
  */

  DelUIVector(&clusters);
  DelMatrix(&centroids);
  DelMatrix(&m);
}

void test16()
{
  matrix *m;
  hgmbins *bins_id;
  HyperGridModel *hgm;

  NewMatrix(&m, 10, 3);

  setMatrixValue(m, 0, 0, 1.9063449);       setMatrixValue(m, 0, 1, 4.7654846);
  setMatrixValue(m, 1, 0, 4.8621579);       setMatrixValue(m, 1, 1, 8.7060234);
  setMatrixValue(m, 2, 0, 3.9597626);       setMatrixValue(m, 2, 1, 4.1309202);
  setMatrixValue(m, 3, 0, 5.7198639);       setMatrixValue(m, 3, 1, 4.8299515);
  setMatrixValue(m, 4, 0, 5.0944097);       setMatrixValue(m, 4, 1, 5.2799994);
  setMatrixValue(m, 5, 0, 9.5844892);       setMatrixValue(m, 5, 1, 3.5483752);
  setMatrixValue(m, 6, 0, 9.5067052);       setMatrixValue(m, 6, 1, 7.0238615);
  setMatrixValue(m, 7, 0, 4.8633205);       setMatrixValue(m, 7, 1, 7.717007);
  setMatrixValue(m, 8, 0, 4.3323122);       setMatrixValue(m, 8, 1, 9.0220105);
  setMatrixValue(m, 9, 0, 6.955182);        setMatrixValue(m, 9, 1, 7.0026549);

  NewHyperGridMap(&hgm);
  HyperGridMap(m, 4, &bins_id, &hgm);
  printf("Total number of bins : %lf\n", hgm->bsize);
  printf("Bins apparteinance id\n");
  PrintHGMBins(bins_id);
  printf("Grid Map\n");
  PrintMatrix(hgm->gmap);
  DelHyperGridMap(&hgm);
  DelMatrix(&m);
  DelHGMBins(&bins_id);
}

void test15()
{
  matrix *m, *centroids;
  int run = SIGSCIENTIFICRUN;
  uivector *clusters;
  strvector *dendogram;
  NewMatrix(&m, 10, 2);

  setMatrixValue(m, 0, 0, 1.9063449);       setMatrixValue(m, 0, 1, 4.7654846);
  setMatrixValue(m, 1, 0, 4.8621579);       setMatrixValue(m, 1, 1, 8.7060234);
  setMatrixValue(m, 2, 0, 3.9597626);       setMatrixValue(m, 2, 1, 4.1309202);
  setMatrixValue(m, 3, 0, 5.7198639);       setMatrixValue(m, 3, 1, 4.8299515);
  setMatrixValue(m, 4, 0, 5.0944097);       setMatrixValue(m, 4, 1, 5.2799994);
  setMatrixValue(m, 5, 0, 9.5844892);       setMatrixValue(m, 5, 1, 3.5483752);
  setMatrixValue(m, 6, 0, 9.5067052);       setMatrixValue(m, 6, 1, 7.0238615);
  setMatrixValue(m, 7, 0, 4.8633205);       setMatrixValue(m, 7, 1, 7.717007);
  setMatrixValue(m, 8, 0, 4.3323122);       setMatrixValue(m, 8, 1, 9.0220105);
  setMatrixValue(m, 9, 0, 6.955182);        setMatrixValue(m, 9, 1, 7.0026549);


  initUIVector(&clusters);
  initMatrix(&centroids);
  initStrVector(&dendogram);
  puts("Hierarchical Clustering");

  HierarchicalClustering(m, 3, &clusters, &centroids, &dendogram, average_linkage, 4, &run);

  puts("clusters");
  PrintUIVector(clusters);
  puts("Dendogram");
  PrintStrVector(dendogram);

  DelStrVector(&dendogram);
  DelUIVector(&clusters);
  DelMatrix(&centroids);
  DelMatrix(&m);
}

void test14()
{
  matrix *m1, *m2;

  double dist;

  NewMatrix(&m1, 10, 2);
  NewMatrix(&m2, 5, 2);

  setMatrixValue(m1, 0, 0, 2);    setMatrixValue(m1, 0, 1, 2);
  setMatrixValue(m1, 1, 0, 2);    setMatrixValue(m1, 1, 1, 5);
  setMatrixValue(m1, 2, 0, 6);    setMatrixValue(m1, 2, 1, 5);
  setMatrixValue(m1, 3, 0, 7);    setMatrixValue(m1, 3, 1, 3);
  setMatrixValue(m1, 4, 0, 4);    setMatrixValue(m1, 4, 1, 7);
  setMatrixValue(m1, 5, 0, 6);    setMatrixValue(m1, 5, 1, 4);
  setMatrixValue(m1, 6, 0, 5);    setMatrixValue(m1, 6, 1, 3);
  setMatrixValue(m1, 7, 0, 4);    setMatrixValue(m1, 7, 1, 6);
  setMatrixValue(m1, 8, 0, 2);    setMatrixValue(m1, 8, 1, 5);
  setMatrixValue(m1, 9, 0, 1);    setMatrixValue(m1, 9, 1, 3);


  setMatrixValue(m2, 0, 0, 6);    setMatrixValue(m2, 0, 1, 5);
  setMatrixValue(m2, 1, 0, 7);    setMatrixValue(m2, 1, 1, 4);
  setMatrixValue(m2, 2, 0, 8);    setMatrixValue(m2, 2, 1, 7);
  setMatrixValue(m2, 3, 0, 5);    setMatrixValue(m2, 3, 1, 6);
  setMatrixValue(m2, 4, 0, 5);    setMatrixValue(m2, 4, 1, 4);

  dist = MatrixMahalanobisDistance(m1, m2);

  printf("dist %f\n", dist);

  DelMatrix(&m1);
  DelMatrix(&m2);
}

void test13()
{
 matrix *m, *centroids;
 uivector *clusters;

 NewMatrix(&m, 6, 1);
 setMatrixValue(m, 0, 0, 1.2);
 setMatrixValue(m, 1, 0, 5.6);
 setMatrixValue(m, 2, 0, 3.7);
 setMatrixValue(m, 3, 0, 0.6);
 setMatrixValue(m, 4, 0, 0.1);
 setMatrixValue(m, 5, 0, 2.6);


  initUIVector(&clusters);
  initMatrix(&centroids);

  puts("KMeans ++");
  KMeans(m, 2, 1, &clusters, &centroids, 4, NULL);

  puts("Centroids");
  PrintMatrix(centroids);

  puts("Selections");
  PrintUIVector(clusters);

  puts("Data Matrix to compare with selections");
  PrintMatrix(m);

  puts("Results: cluster 1 = 5.6, 3.7   cluster 2 = 1.2, 0.6, 0.1, 2.6");

  DelUIVector(&clusters);
  DelMatrix(&centroids);
  DelMatrix(&m);
}

void test12()
{
  size_t i, j;
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;
  dvector *ssdist;

  NewMatrix(&m, 100, 2);
  for(i = 0; i < 100; i++){
    for(j = 0; j < 2; j++){
      srand(time(0)+i+j+rand());
      setMatrixValue(m, i, j, rand() % 100);
    }
  }

  initDVector(&ssdist);
  puts("KMeans ++ Cross Validation");
  KMeansRandomGroupsCV(m, 15, 1, 3, 10, &ssdist, 4, &run);

  puts("ssdist");
  PrintDVector(ssdist);

  DelDVector(&ssdist);
  DelMatrix(&m);
}

void test11()
{
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;
  dvector *ssdist;

  NewMatrix(&m, 14, 2);

  m->data[0][0] =  84.1400;  m->data[0][1] =  357.1500;
  m->data[1][0] =  79.1000;  m->data[1][1] = 231.0000;
  m->data[2][0] =  67.0900;  m->data[2][1] = 403.0000;
  m->data[3][0] = 68.0700;  m->data[3][1] = 304.5500;
  m->data[4][0] = 68.0800;  m->data[4][1] = 529.0000;
  m->data[5][0] = 129.1600;  m->data[5][1] =  510.0000;
  m->data[6][0] = 128.1600;  m->data[6][1] =  491.0000;
  m->data[7][0] = 78.1118;  m->data[7][1] =  353.3000;
  m->data[8][0] = 202.2550;  m->data[8][1] =  666.6500;
  m->data[9][0] = 84.1600;  m->data[9][1] = 354.0000;
  m->data[10][0] = 72.1100;  m->data[10][1] = 339.0000;
  m->data[11][0] = 71.1100;  m->data[11][1] = 360.0000;
  m->data[12][0] = 85.1500;  m->data[12][1] = 379.0000;
  m->data[13][0] = 86.1300;  m->data[13][1] =  361.0000;

  initDVector(&ssdist);
  puts("KMeans ++ Cross Validation");
  /*KMeansRandomGroupsCV(m, 9, 1, 3, 20, &ssdist, &run);*/

  KMeansJumpMethod(m, 9, 1, &ssdist, 4, &run);

  puts("ssdist");
  PrintDVector(ssdist);

  DelDVector(&ssdist);
  DelMatrix(&m);
}

void test10()
{
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;
  uivector *clusters;
  matrix *centroids;

  NewMatrix(&m, 14, 2);

  m->data[0][0] =  84.1400;  m->data[0][1] =  357.1500;
  m->data[1][0] =  79.1000;  m->data[1][1] = 231.0000;
  m->data[2][0] =  67.0900;  m->data[2][1] = 403.0000;
  m->data[3][0] = 68.0700;  m->data[3][1] = 304.5500;
  m->data[4][0] = 68.0800;  m->data[4][1] = 529.0000;
  m->data[5][0] = 129.1600;  m->data[5][1] =  510.0000;
  m->data[6][0] = 128.1600;  m->data[6][1] =  491.0000;
  m->data[7][0] = 78.1118;  m->data[7][1] =  353.3000;
  m->data[8][0] = 202.2550;  m->data[8][1] =  666.6500;
  m->data[9][0] = 84.1600;  m->data[9][1] = 354.0000;
  m->data[10][0] = 72.1100;  m->data[10][1] = 339.0000;
  m->data[11][0] = 71.1100;  m->data[11][1] = 360.0000;
  m->data[12][0] = 85.1500;  m->data[12][1] = 379.0000;
  m->data[13][0] = 86.1300;  m->data[13][1] =  361.0000;

  initUIVector(&clusters);
  initMatrix(&centroids);

  puts("KMeans ++");
  KMeans(m, 3, 1, &clusters, &centroids, 4, &run);

  puts("Centroids");
  PrintMatrix(centroids);

  puts("Selections");
  PrintUIVector(clusters);

  //PruneResults(m, centroids, 3, 0, clusters);

  /*puts("New Selections");
  PrintUIVector(clusters);*/

  DelUIVector(&clusters);
  DelMatrix(&centroids);
  DelMatrix(&m);
}

void test9()
{
  matrix *m; /* Data matrix */
  uivector *clusters;
  int run = SIGSCIENTIFICRUN;

  NewMatrix(&m, 100, 2);

  setMatrixValue(m, 0, 0, 7.651165);  setMatrixValue(m, 0, 1, 33.374403);
  setMatrixValue(m, 1, 0, 6.218197);  setMatrixValue(m, 1, 1, 69.393923);
  setMatrixValue(m, 2, 0, 9.968406);  setMatrixValue(m, 2, 1, 4.948435);
  setMatrixValue(m, 3, 0, 5.714325);  setMatrixValue(m, 3, 1, 69.416728);
  setMatrixValue(m, 4, 0, 4.924057);  setMatrixValue(m, 4, 1, 12.706405);
  setMatrixValue(m, 5, 0, 4.825453);  setMatrixValue(m, 5, 1, 53.826898);
  setMatrixValue(m, 6, 0, 2.891843);  setMatrixValue(m, 6, 1, 24.575599);
  setMatrixValue(m, 7, 0, 5.960617);  setMatrixValue(m, 7, 1, 42.125457);
  setMatrixValue(m, 8, 0, 7.671809);  setMatrixValue(m, 8, 1, 21.719972);
  setMatrixValue(m, 9, 0, 3.438646);  setMatrixValue(m, 9, 1, 92.106844);
  setMatrixValue(m, 10, 0, 9.475234);  setMatrixValue(m, 10, 1, 4.049699);
  setMatrixValue(m, 11, 0, 8.626795);  setMatrixValue(m, 11, 1, 26.943452);
  setMatrixValue(m, 12, 0, 5.860959);  setMatrixValue(m, 12, 1, 37.790985);
  setMatrixValue(m, 13, 0, 6.124417);  setMatrixValue(m, 13, 1, 19.057411);
  setMatrixValue(m, 14, 0, 1.734119);  setMatrixValue(m, 14, 1, 91.159850);
  setMatrixValue(m, 15, 0, 3.572315);  setMatrixValue(m, 15, 1, 52.605869);
  setMatrixValue(m, 16, 0, 9.087432);  setMatrixValue(m, 16, 1, 20.092922);
  setMatrixValue(m, 17, 0, 4.905813);  setMatrixValue(m, 17, 1, 53.304572);
  setMatrixValue(m, 18, 0, 9.539760);  setMatrixValue(m, 18, 1, 51.435908);
  setMatrixValue(m, 19, 0, 6.959215);  setMatrixValue(m, 19, 1, 54.982246);
  setMatrixValue(m, 20, 0, 8.386170);  setMatrixValue(m, 20, 1, 30.504135);
  setMatrixValue(m, 21, 0, 7.268718);  setMatrixValue(m, 21, 1, 87.612271);
  setMatrixValue(m, 22, 0, 2.423585);  setMatrixValue(m, 22, 1, 31.551512);
  setMatrixValue(m, 23, 0, 9.624445);  setMatrixValue(m, 23, 1, 78.476988);
  setMatrixValue(m, 24, 0, 4.433680);  setMatrixValue(m, 24, 1, 71.037004);
  setMatrixValue(m, 25, 0, 6.798443);  setMatrixValue(m, 25, 1, 25.033095);
  setMatrixValue(m, 26, 0, 9.900906);  setMatrixValue(m, 26, 1, 81.120367);
  setMatrixValue(m, 27, 0, 3.968972);  setMatrixValue(m, 27, 1, 88.536517);
  setMatrixValue(m, 28, 0, 8.022121);  setMatrixValue(m, 28, 1, 65.349167);
  setMatrixValue(m, 29, 0, 8.026969);  setMatrixValue(m, 29, 1, 52.926456);
  setMatrixValue(m, 30, 0, 6.413291);  setMatrixValue(m, 30, 1, 28.124218);
  setMatrixValue(m, 31, 0, 1.879453);  setMatrixValue(m, 31, 1, 75.084904);
  setMatrixValue(m, 32, 0, 9.296220);  setMatrixValue(m, 32, 1, 76.371613);
  setMatrixValue(m, 33, 0, 3.920501);  setMatrixValue(m, 33, 1, 36.414331);
  setMatrixValue(m, 34, 0, 1.092812);  setMatrixValue(m, 34, 1, 56.834326);
  setMatrixValue(m, 35, 0, 8.201150);  setMatrixValue(m, 35, 1, 59.379607);
  setMatrixValue(m, 36, 0, 4.397655);  setMatrixValue(m, 36, 1, 61.071957);
  setMatrixValue(m, 37, 0, 4.220824);  setMatrixValue(m, 37, 1, 70.393470);
  setMatrixValue(m, 38, 0, 9.733481);  setMatrixValue(m, 38, 1, 42.032896);
  setMatrixValue(m, 39, 0, 3.998391);  setMatrixValue(m, 39, 1, 46.566763);
  setMatrixValue(m, 40, 0, 1.784019);  setMatrixValue(m, 40, 1, 46.779213);
  setMatrixValue(m, 41, 0, 7.535995);  setMatrixValue(m, 41, 1, 93.725636);
  setMatrixValue(m, 42, 0, 7.394030);  setMatrixValue(m, 42, 1, 44.959715);
  setMatrixValue(m, 43, 0, 9.304112);  setMatrixValue(m, 43, 1, 12.562119);
  setMatrixValue(m, 44, 0, 5.433701);  setMatrixValue(m, 44, 1, 30.377871);
  setMatrixValue(m, 45, 0, 2.474974);  setMatrixValue(m, 45, 1, 89.772628);
  setMatrixValue(m, 46, 0, 3.583825);  setMatrixValue(m, 46, 1, 75.089247);
  setMatrixValue(m, 47, 0, 5.967499);  setMatrixValue(m, 47, 1, 36.262162);
  setMatrixValue(m, 48, 0, 5.533060);  setMatrixValue(m, 48, 1, 90.168083);
  setMatrixValue(m, 49, 0, 2.588625);  setMatrixValue(m, 49, 1, 89.612832);
  setMatrixValue(m, 50, 0, 9.932260);  setMatrixValue(m, 50, 1, 8.981204);
  setMatrixValue(m, 51, 0, 3.158407);  setMatrixValue(m, 51, 1, 34.050263);
  setMatrixValue(m, 52, 0, 2.131164);  setMatrixValue(m, 52, 1, 43.370761);
  setMatrixValue(m, 53, 0, 5.594023);  setMatrixValue(m, 53, 1, 54.449151);
  setMatrixValue(m, 54, 0, 8.958604);  setMatrixValue(m, 54, 1, 69.407942);
  setMatrixValue(m, 55, 0, 6.024750);  setMatrixValue(m, 55, 1, 46.827483);
  setMatrixValue(m, 56, 0, 6.418956);  setMatrixValue(m, 56, 1, 53.286935);
  setMatrixValue(m, 57, 0, 5.780031);  setMatrixValue(m, 57, 1, 5.008359);
  setMatrixValue(m, 58, 0, 3.422119);  setMatrixValue(m, 58, 1, 83.216003);
  setMatrixValue(m, 59, 0, 9.052767);  setMatrixValue(m, 59, 1, 74.315114);
  setMatrixValue(m, 60, 0, 8.601263);  setMatrixValue(m, 60, 1, 58.054798);
  setMatrixValue(m, 61, 0, 3.998506);  setMatrixValue(m, 61, 1, 88.703496);
  setMatrixValue(m, 62, 0, 7.667506);  setMatrixValue(m, 62, 1, 17.890521);
  setMatrixValue(m, 63, 0, 2.337662);  setMatrixValue(m, 63, 1, 14.210792);
  setMatrixValue(m, 64, 0, 5.156890);  setMatrixValue(m, 64, 1, 80.958857);
  setMatrixValue(m, 65, 0, 2.018825);  setMatrixValue(m, 65, 1, 97.187163);
  setMatrixValue(m, 66, 0, 4.432068);  setMatrixValue(m, 66, 1, 83.327432);
  setMatrixValue(m, 67, 0, 2.114427);  setMatrixValue(m, 67, 1, 7.031047);
  setMatrixValue(m, 68, 0, 5.496374);  setMatrixValue(m, 68, 1, 4.867239);
  setMatrixValue(m, 69, 0, 9.673838);  setMatrixValue(m, 69, 1, 68.282363);
  setMatrixValue(m, 70, 0, 3.409410);  setMatrixValue(m, 70, 1, 33.769704);
  setMatrixValue(m, 71, 0, 2.067344);  setMatrixValue(m, 71, 1, 86.043228);
  setMatrixValue(m, 72, 0, 3.522978);  setMatrixValue(m, 72, 1, 19.634175);
  setMatrixValue(m, 73, 0, 6.631368);  setMatrixValue(m, 73, 1, 6.293639);
  setMatrixValue(m, 74, 0, 8.929563);  setMatrixValue(m, 74, 1, 58.177679);
  setMatrixValue(m, 75, 0, 1.009022);  setMatrixValue(m, 75, 1, 19.041718);
  setMatrixValue(m, 76, 0, 2.616159);  setMatrixValue(m, 76, 1, 84.733896);
  setMatrixValue(m, 77, 0, 5.480356);  setMatrixValue(m, 77, 1, 8.508343);
  setMatrixValue(m, 78, 0, 4.365632);  setMatrixValue(m, 78, 1, 18.537420);
  setMatrixValue(m, 79, 0, 2.725384);  setMatrixValue(m, 79, 1, 75.165470);
  setMatrixValue(m, 80, 0, 9.175719);  setMatrixValue(m, 80, 1, 87.846565);
  setMatrixValue(m, 81, 0, 3.245108);  setMatrixValue(m, 81, 1, 91.361911);
  setMatrixValue(m, 82, 0, 4.456010);  setMatrixValue(m, 82, 1, 69.457073);
  setMatrixValue(m, 83, 0, 2.896024);  setMatrixValue(m, 83, 1, 79.771741);
  setMatrixValue(m, 84, 0, 3.258866);  setMatrixValue(m, 84, 1, 87.442361);
  setMatrixValue(m, 85, 0, 3.022637);  setMatrixValue(m, 85, 1, 69.437017);
  setMatrixValue(m, 86, 0, 5.552419);  setMatrixValue(m, 86, 1, 9.583123);
  setMatrixValue(m, 87, 0, 9.192921);  setMatrixValue(m, 87, 1, 13.282271);
  setMatrixValue(m, 88, 0, 2.593360);  setMatrixValue(m, 88, 1, 59.944635);
  setMatrixValue(m, 89, 0, 6.159689);  setMatrixValue(m, 89, 1, 16.881533);
  setMatrixValue(m, 90, 0, 7.189455);  setMatrixValue(m, 90, 1, 64.119091);
  setMatrixValue(m, 91, 0, 5.082821);  setMatrixValue(m, 91, 1, 84.611977);
  setMatrixValue(m, 92, 0, 6.116967);  setMatrixValue(m, 92, 1, 89.405369);
  setMatrixValue(m, 93, 0, 7.255329);  setMatrixValue(m, 93, 1, 62.257379);
  setMatrixValue(m, 94, 0, 6.905527);  setMatrixValue(m, 94, 1, 73.774295);
  setMatrixValue(m, 95, 0, 6.180253);  setMatrixValue(m, 95, 1, 25.157409);
  setMatrixValue(m, 96, 0, 4.670816);  setMatrixValue(m, 96, 1, 11.266562);
  setMatrixValue(m, 97, 0, 7.109005);  setMatrixValue(m, 97, 1, 11.772769);
  setMatrixValue(m, 98, 0, 3.113531);  setMatrixValue(m, 98, 1, 60.972161);
  setMatrixValue(m, 99, 0, 3.855047);  setMatrixValue(m, 99, 1, 91.930780);

  initUIVector(&clusters);


  KMeans(m, 10, 3, &clusters, NULL, 4, &run);

  puts("Selections");
  PrintUIVector(clusters);

  DelUIVector(&clusters);
  DelMatrix(&m);
}

void test8()
{
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;
  uivector *clusters;
  matrix *centroids;

  NewMatrix(&m, 14, 2);

  m->data[0][0] =  84.1400;  m->data[0][1] =  357.1500;
  m->data[1][0] =  79.1000;  m->data[1][1] = 231.0000;
  m->data[2][0] =  67.0900;  m->data[2][1] = 403.0000;
  m->data[3][0] = 68.0700;  m->data[3][1] = 304.5500;
  m->data[4][0] = 68.0800;  m->data[4][1] = 529.0000;
  m->data[5][0] = 129.1600;  m->data[5][1] =  510.0000;
  m->data[6][0] = 128.1600;  m->data[6][1] =  491.0000;
  m->data[7][0] = 78.1118;  m->data[7][1] =  353.3000;
  m->data[8][0] = 202.2550;  m->data[8][1] =  666.6500;
  m->data[9][0] = 84.1600;  m->data[9][1] = 354.0000;
  m->data[10][0] = 72.1100;  m->data[10][1] = 339.0000;
  m->data[11][0] = 71.1100;  m->data[11][1] = 360.0000;
  m->data[12][0] = 85.1500;  m->data[12][1] = 379.0000;
  m->data[13][0] = 86.1300;  m->data[13][1] =  361.0000;

  initUIVector(&clusters);
  initMatrix(&centroids);

  puts("KMeans MaxDis");
  KMeans(m, 2, 3, &clusters, &centroids, 4, &run);

  puts("Selections");
  PrintUIVector(clusters);

  PrintMatrix(centroids);

  DelUIVector(&clusters);
  DelMatrix(&centroids);
  DelMatrix(&m);
}

void test7()
{
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;
  uivector *clusters;
  matrix *centroids;

  NewMatrix(&m, 14, 2);

  m->data[0][0] =  84.1400;  m->data[0][1] =  357.1500;
  m->data[1][0] =  79.1000;  m->data[1][1] = 231.0000;
  m->data[2][0] =  67.0900;  m->data[2][1] = 403.0000;
  m->data[3][0] = 68.0700;  m->data[3][1] = 304.5500;
  m->data[4][0] = 68.0800;  m->data[4][1] = 529.0000;
  m->data[5][0] = 129.1600;  m->data[5][1] =  510.0000;
  m->data[6][0] = 128.1600;  m->data[6][1] =  491.0000;
  m->data[7][0] = 78.1118;  m->data[7][1] =  353.3000;
  m->data[8][0] = 202.2550;  m->data[8][1] =  666.6500;
  m->data[9][0] = 84.1600;  m->data[9][1] = 354.0000;
  m->data[10][0] = 72.1100;  m->data[10][1] = 339.0000;
  m->data[11][0] = 71.1100;  m->data[11][1] = 360.0000;
  m->data[12][0] = 85.1500;  m->data[12][1] = 379.0000;
  m->data[13][0] = 86.1300;  m->data[13][1] =  361.0000;

  initUIVector(&clusters);
  initMatrix(&centroids);

  puts("KMeans MDC");
  KMeans(m, 2, 2, &clusters, &centroids, 4, &run);

  puts("Selections");
  PrintUIVector(clusters);

  PrintMatrix(centroids);

  DelUIVector(&clusters);
  DelMatrix(&centroids);
  DelMatrix(&m);
}

void test6()
{
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;
  uivector *clusters;
  matrix *centroids;

  NewMatrix(&m, 14, 2);

  m->data[0][0] =  84.1400;  m->data[0][1] =  357.1500;
  m->data[1][0] =  79.1000;  m->data[1][1] = 231.0000;
  m->data[2][0] =  67.0900;  m->data[2][1] = 403.0000;
  m->data[3][0] = 68.0700;  m->data[3][1] = 304.5500;
  m->data[4][0] = 68.0800;  m->data[4][1] = 529.0000;
  m->data[5][0] = 129.1600;  m->data[5][1] =  510.0000;
  m->data[6][0] = 128.1600;  m->data[6][1] =  491.0000;
  m->data[7][0] = 78.1118;  m->data[7][1] =  353.3000;
  m->data[8][0] = 202.2550;  m->data[8][1] =  666.6500;
  m->data[9][0] = 84.1600;  m->data[9][1] = 354.0000;
  m->data[10][0] = 72.1100;  m->data[10][1] = 339.0000;
  m->data[11][0] = 71.1100;  m->data[11][1] = 360.0000;
  m->data[12][0] = 85.1500;  m->data[12][1] = 379.0000;
  m->data[13][0] = 86.1300;  m->data[13][1] =  361.0000;

  initUIVector(&clusters);
  initMatrix(&centroids);

  puts("KMeans RANDOM");
  KMeans(m, 2, 0, &clusters, &centroids, 4, &run);

  puts("Selections");
  PrintUIVector(clusters);

  PrintMatrix(centroids);

  DelUIVector(&clusters);
  DelMatrix(&centroids);
  DelMatrix(&m);
}

void test5()
{
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;
  uivector *clusters;
  matrix *centroids;

  NewMatrix(&m, 14, 2);

  m->data[0][0] =  84.1400;  m->data[0][1] =  357.1500;
  m->data[1][0] =  79.1000;  m->data[1][1] = 231.0000;
  m->data[2][0] =  67.0900;  m->data[2][1] = 403.0000;
  m->data[3][0] = 68.0700;  m->data[3][1] = 304.5500;
  m->data[4][0] = 68.0800;  m->data[4][1] = 529.0000;
  m->data[5][0] = 129.1600;  m->data[5][1] =  510.0000;
  m->data[6][0] = 128.1600;  m->data[6][1] =  491.0000;
  m->data[7][0] = 78.1118;  m->data[7][1] =  353.3000;
  m->data[8][0] = 202.2550;  m->data[8][1] =  666.6500;
  m->data[9][0] = 84.1600;  m->data[9][1] = 354.0000;
  m->data[10][0] = 72.1100;  m->data[10][1] = 339.0000;
  m->data[11][0] = 71.1100;  m->data[11][1] = 360.0000;
  m->data[12][0] = 85.1500;  m->data[12][1] = 379.0000;
  m->data[13][0] = 86.1300;  m->data[13][1] =  361.0000;

  initUIVector(&clusters);
  initMatrix(&centroids);

  puts("KMeans++ TEST");
  KMeans(m, 2, 1, &clusters, &centroids, 4, &run);

  puts("Selections");
  PrintUIVector(clusters);

  puts("Centroids");
  PrintMatrix(centroids);

  DelUIVector(&clusters);
  DelMatrix(&centroids);
  DelMatrix(&m);
}

void test4()
{
  size_t i;
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;
  uivector *selections;

  NewMatrix(&m, 14, 2);

  m->data[0][0] =  84.1400;  m->data[0][1] =  357.1500;
  m->data[1][0] =  79.1000;  m->data[1][1] = 231.0000;
  m->data[2][0] =  67.0900;  m->data[2][1] = 403.0000;
  m->data[3][0] = 68.0700;  m->data[3][1] = 304.5500;
  m->data[4][0] = 68.0800;  m->data[4][1] = 529.0000;
  m->data[5][0] = 129.1600;  m->data[5][1] =  510.0000;
  m->data[6][0] = 128.1600;  m->data[6][1] =  491.0000;
  m->data[7][0] = 78.1118;  m->data[7][1] =  353.3000;
  m->data[8][0] = 202.2550;  m->data[8][1] =  666.6500;
  m->data[9][0] = 84.1600;  m->data[9][1] = 354.0000;
  m->data[10][0] = 72.1100;  m->data[10][1] = 339.0000;
  m->data[11][0] = 71.1100;  m->data[11][1] = 360.0000;
  m->data[12][0] = 85.1500;  m->data[12][1] = 379.0000;
  m->data[13][0] = 86.1300;  m->data[13][1] =  361.0000;

  initUIVector(&selections);

  puts("KMeansppCenters TEST");
  KMeansppCenters(m, 2, &selections, 1, &run);

  puts("Selections");
  PrintUIVector(selections);

  for(i = 0; i < selections->size; i++){
    printf("%f\t%f\n", getMatrixValue(m, getUIVectorValue(selections, i), 0), getMatrixValue(m, getUIVectorValue(selections, i), 1));
  }

  DelUIVector(&selections);
  DelMatrix(&m);
}

void test3()
{
  size_t i;
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;
  uivector *selections;

  NewMatrix(&m, 14, 2);


  m->data[0][0] =  84.1400;  m->data[0][1] =  357.1500;
  m->data[1][0] =  79.1000;  m->data[1][1] = 231.0000;
  m->data[2][0] =  67.0900;  m->data[2][1] = 403.0000;
  m->data[3][0] = 68.0700;  m->data[3][1] = 304.5500;
  m->data[4][0] = 68.0800;  m->data[4][1] = 529.0000;
  m->data[5][0] = 129.1600;  m->data[5][1] =  510.0000;
  m->data[6][0] = 128.1600;  m->data[6][1] =  491.0000;
  m->data[7][0] = 78.1118;  m->data[7][1] =  353.3000;
  m->data[8][0] = 202.2550;  m->data[8][1] =  666.6500;
  m->data[9][0] = 84.1600;  m->data[9][1] = 354.0000;
  m->data[10][0] = 72.1100;  m->data[10][1] = 339.0000;
  m->data[11][0] = 71.1100;  m->data[11][1] = 360.0000;
  m->data[12][0] = 85.1500;  m->data[12][1] = 379.0000;
  m->data[13][0] = 86.1300;  m->data[13][1] =  361.0000;

  PrintMatrix(m);
  initUIVector(&selections);

  MaxDis(m, floor(0.35*m->row), 0, &selections, 4, &run);

  puts("Selections");
  PrintUIVector(selections);

  for(i = 0; i < selections->size; i++){
    printf("%f\t%f\n", getMatrixValue(m, getUIVectorValue(selections, i), 0), getMatrixValue(m, getUIVectorValue(selections, i), 1));
  }

  DelUIVector(&selections);
  DelMatrix(&m);
}

void test2()
{
  size_t i;
  matrix *m; /* Data matrix */
  int run = SIGSCIENTIFICRUN;

  uivector *selections;

  NewMatrix(&m, 14, 2);


  m->data[0][0] =  84.1400;  m->data[0][1] =  357.1500;
  m->data[1][0] =  79.1000;  m->data[1][1] = 231.0000;
  m->data[2][0] =  67.0900;  m->data[2][1] = 403.0000;
  m->data[3][0] = 68.0700;  m->data[3][1] = 304.5500;
  m->data[4][0] = 68.0800;  m->data[4][1] = 529.0000;
  m->data[5][0] = 129.1600;  m->data[5][1] =  510.0000;
  m->data[6][0] = 128.1600;  m->data[6][1] =  491.0000;
  m->data[7][0] = 78.1118;  m->data[7][1] =  353.3000;
  m->data[8][0] = 202.2550;  m->data[8][1] =  666.6500;
  m->data[9][0] = 84.1600;  m->data[9][1] = 354.0000;
  m->data[10][0] = 72.1100;  m->data[10][1] = 339.0000;
  m->data[11][0] = 71.1100;  m->data[11][1] = 360.0000;
  m->data[12][0] = 85.1500;  m->data[12][1] = 379.0000;
  m->data[13][0] = 86.1300;  m->data[13][1] =  361.0000;

  initUIVector(&selections);

  MDC(m, 0, 0, &selections, 4, &run);

  puts("Selections");
  PrintUIVector(selections);

  for(i = 0; i < selections->size; i++){
    printf("%f\t%f\n", getMatrixValue(m, getUIVectorValue(selections, i), 0), getMatrixValue(m, getUIVectorValue(selections, i), 1));
  }

  DelUIVector(&selections);
  DelMatrix(&m);
}


void test1()
{
  size_t i;
  matrix *m; /* Data matrix */
  uivector *selections;
  int run = SIGSCIENTIFICRUN;

  NewMatrix(&m, 100, 2);

  setMatrixValue(m, 0, 0, 7.651165);  setMatrixValue(m, 0, 1, 33.374403);
  setMatrixValue(m, 1, 0, 6.218197);  setMatrixValue(m, 1, 1, 69.393923);
  setMatrixValue(m, 2, 0, 9.968406);  setMatrixValue(m, 2, 1, 4.948435);
  setMatrixValue(m, 3, 0, 5.714325);  setMatrixValue(m, 3, 1, 69.416728);
  setMatrixValue(m, 4, 0, 4.924057);  setMatrixValue(m, 4, 1, 12.706405);
  setMatrixValue(m, 5, 0, 4.825453);  setMatrixValue(m, 5, 1, 53.826898);
  setMatrixValue(m, 6, 0, 2.891843);  setMatrixValue(m, 6, 1, 24.575599);
  setMatrixValue(m, 7, 0, 5.960617);  setMatrixValue(m, 7, 1, 42.125457);
  setMatrixValue(m, 8, 0, 7.671809);  setMatrixValue(m, 8, 1, 21.719972);
  setMatrixValue(m, 9, 0, 3.438646);  setMatrixValue(m, 9, 1, 92.106844);
  setMatrixValue(m, 10, 0, 9.475234);  setMatrixValue(m, 10, 1, 4.049699);
  setMatrixValue(m, 11, 0, 8.626795);  setMatrixValue(m, 11, 1, 26.943452);
  setMatrixValue(m, 12, 0, 5.860959);  setMatrixValue(m, 12, 1, 37.790985);
  setMatrixValue(m, 13, 0, 6.124417);  setMatrixValue(m, 13, 1, 19.057411);
  setMatrixValue(m, 14, 0, 1.734119);  setMatrixValue(m, 14, 1, 91.159850);
  setMatrixValue(m, 15, 0, 3.572315);  setMatrixValue(m, 15, 1, 52.605869);
  setMatrixValue(m, 16, 0, 9.087432);  setMatrixValue(m, 16, 1, 20.092922);
  setMatrixValue(m, 17, 0, 4.905813);  setMatrixValue(m, 17, 1, 53.304572);
  setMatrixValue(m, 18, 0, 9.539760);  setMatrixValue(m, 18, 1, 51.435908);
  setMatrixValue(m, 19, 0, 6.959215);  setMatrixValue(m, 19, 1, 54.982246);
  setMatrixValue(m, 20, 0, 8.386170);  setMatrixValue(m, 20, 1, 30.504135);
  setMatrixValue(m, 21, 0, 7.268718);  setMatrixValue(m, 21, 1, 87.612271);
  setMatrixValue(m, 22, 0, 2.423585);  setMatrixValue(m, 22, 1, 31.551512);
  setMatrixValue(m, 23, 0, 9.624445);  setMatrixValue(m, 23, 1, 78.476988);
  setMatrixValue(m, 24, 0, 4.433680);  setMatrixValue(m, 24, 1, 71.037004);
  setMatrixValue(m, 25, 0, 6.798443);  setMatrixValue(m, 25, 1, 25.033095);
  setMatrixValue(m, 26, 0, 9.900906);  setMatrixValue(m, 26, 1, 81.120367);
  setMatrixValue(m, 27, 0, 3.968972);  setMatrixValue(m, 27, 1, 88.536517);
  setMatrixValue(m, 28, 0, 8.022121);  setMatrixValue(m, 28, 1, 65.349167);
  setMatrixValue(m, 29, 0, 8.026969);  setMatrixValue(m, 29, 1, 52.926456);
  setMatrixValue(m, 30, 0, 6.413291);  setMatrixValue(m, 30, 1, 28.124218);
  setMatrixValue(m, 31, 0, 1.879453);  setMatrixValue(m, 31, 1, 75.084904);
  setMatrixValue(m, 32, 0, 9.296220);  setMatrixValue(m, 32, 1, 76.371613);
  setMatrixValue(m, 33, 0, 3.920501);  setMatrixValue(m, 33, 1, 36.414331);
  setMatrixValue(m, 34, 0, 1.092812);  setMatrixValue(m, 34, 1, 56.834326);
  setMatrixValue(m, 35, 0, 8.201150);  setMatrixValue(m, 35, 1, 59.379607);
  setMatrixValue(m, 36, 0, 4.397655);  setMatrixValue(m, 36, 1, 61.071957);
  setMatrixValue(m, 37, 0, 4.220824);  setMatrixValue(m, 37, 1, 70.393470);
  setMatrixValue(m, 38, 0, 9.733481);  setMatrixValue(m, 38, 1, 42.032896);
  setMatrixValue(m, 39, 0, 3.998391);  setMatrixValue(m, 39, 1, 46.566763);
  setMatrixValue(m, 40, 0, 1.784019);  setMatrixValue(m, 40, 1, 46.779213);
  setMatrixValue(m, 41, 0, 7.535995);  setMatrixValue(m, 41, 1, 93.725636);
  setMatrixValue(m, 42, 0, 7.394030);  setMatrixValue(m, 42, 1, 44.959715);
  setMatrixValue(m, 43, 0, 9.304112);  setMatrixValue(m, 43, 1, 12.562119);
  setMatrixValue(m, 44, 0, 5.433701);  setMatrixValue(m, 44, 1, 30.377871);
  setMatrixValue(m, 45, 0, 2.474974);  setMatrixValue(m, 45, 1, 89.772628);
  setMatrixValue(m, 46, 0, 3.583825);  setMatrixValue(m, 46, 1, 75.089247);
  setMatrixValue(m, 47, 0, 5.967499);  setMatrixValue(m, 47, 1, 36.262162);
  setMatrixValue(m, 48, 0, 5.533060);  setMatrixValue(m, 48, 1, 90.168083);
  setMatrixValue(m, 49, 0, 2.588625);  setMatrixValue(m, 49, 1, 89.612832);
  setMatrixValue(m, 50, 0, 9.932260);  setMatrixValue(m, 50, 1, 8.981204);
  setMatrixValue(m, 51, 0, 3.158407);  setMatrixValue(m, 51, 1, 34.050263);
  setMatrixValue(m, 52, 0, 2.131164);  setMatrixValue(m, 52, 1, 43.370761);
  setMatrixValue(m, 53, 0, 5.594023);  setMatrixValue(m, 53, 1, 54.449151);
  setMatrixValue(m, 54, 0, 8.958604);  setMatrixValue(m, 54, 1, 69.407942);
  setMatrixValue(m, 55, 0, 6.024750);  setMatrixValue(m, 55, 1, 46.827483);
  setMatrixValue(m, 56, 0, 6.418956);  setMatrixValue(m, 56, 1, 53.286935);
  setMatrixValue(m, 57, 0, 5.780031);  setMatrixValue(m, 57, 1, 5.008359);
  setMatrixValue(m, 58, 0, 3.422119);  setMatrixValue(m, 58, 1, 83.216003);
  setMatrixValue(m, 59, 0, 9.052767);  setMatrixValue(m, 59, 1, 74.315114);
  setMatrixValue(m, 60, 0, 8.601263);  setMatrixValue(m, 60, 1, 58.054798);
  setMatrixValue(m, 61, 0, 3.998506);  setMatrixValue(m, 61, 1, 88.703496);
  setMatrixValue(m, 62, 0, 7.667506);  setMatrixValue(m, 62, 1, 17.890521);
  setMatrixValue(m, 63, 0, 2.337662);  setMatrixValue(m, 63, 1, 14.210792);
  setMatrixValue(m, 64, 0, 5.156890);  setMatrixValue(m, 64, 1, 80.958857);
  setMatrixValue(m, 65, 0, 2.018825);  setMatrixValue(m, 65, 1, 97.187163);
  setMatrixValue(m, 66, 0, 4.432068);  setMatrixValue(m, 66, 1, 83.327432);
  setMatrixValue(m, 67, 0, 2.114427);  setMatrixValue(m, 67, 1, 7.031047);
  setMatrixValue(m, 68, 0, 5.496374);  setMatrixValue(m, 68, 1, 4.867239);
  setMatrixValue(m, 69, 0, 9.673838);  setMatrixValue(m, 69, 1, 68.282363);
  setMatrixValue(m, 70, 0, 3.409410);  setMatrixValue(m, 70, 1, 33.769704);
  setMatrixValue(m, 71, 0, 2.067344);  setMatrixValue(m, 71, 1, 86.043228);
  setMatrixValue(m, 72, 0, 3.522978);  setMatrixValue(m, 72, 1, 19.634175);
  setMatrixValue(m, 73, 0, 6.631368);  setMatrixValue(m, 73, 1, 6.293639);
  setMatrixValue(m, 74, 0, 8.929563);  setMatrixValue(m, 74, 1, 58.177679);
  setMatrixValue(m, 75, 0, 1.009022);  setMatrixValue(m, 75, 1, 19.041718);
  setMatrixValue(m, 76, 0, 2.616159);  setMatrixValue(m, 76, 1, 84.733896);
  setMatrixValue(m, 77, 0, 5.480356);  setMatrixValue(m, 77, 1, 8.508343);
  setMatrixValue(m, 78, 0, 4.365632);  setMatrixValue(m, 78, 1, 18.537420);
  setMatrixValue(m, 79, 0, 2.725384);  setMatrixValue(m, 79, 1, 75.165470);
  setMatrixValue(m, 80, 0, 9.175719);  setMatrixValue(m, 80, 1, 87.846565);
  setMatrixValue(m, 81, 0, 3.245108);  setMatrixValue(m, 81, 1, 91.361911);
  setMatrixValue(m, 82, 0, 4.456010);  setMatrixValue(m, 82, 1, 69.457073);
  setMatrixValue(m, 83, 0, 2.896024);  setMatrixValue(m, 83, 1, 79.771741);
  setMatrixValue(m, 84, 0, 3.258866);  setMatrixValue(m, 84, 1, 87.442361);
  setMatrixValue(m, 85, 0, 3.022637);  setMatrixValue(m, 85, 1, 69.437017);
  setMatrixValue(m, 86, 0, 5.552419);  setMatrixValue(m, 86, 1, 9.583123);
  setMatrixValue(m, 87, 0, 9.192921);  setMatrixValue(m, 87, 1, 13.282271);
  setMatrixValue(m, 88, 0, 2.593360);  setMatrixValue(m, 88, 1, 59.944635);
  setMatrixValue(m, 89, 0, 6.159689);  setMatrixValue(m, 89, 1, 16.881533);
  setMatrixValue(m, 90, 0, 7.189455);  setMatrixValue(m, 90, 1, 64.119091);
  setMatrixValue(m, 91, 0, 5.082821);  setMatrixValue(m, 91, 1, 84.611977);
  setMatrixValue(m, 92, 0, 6.116967);  setMatrixValue(m, 92, 1, 89.405369);
  setMatrixValue(m, 93, 0, 7.255329);  setMatrixValue(m, 93, 1, 62.257379);
  setMatrixValue(m, 94, 0, 6.905527);  setMatrixValue(m, 94, 1, 73.774295);
  setMatrixValue(m, 95, 0, 6.180253);  setMatrixValue(m, 95, 1, 25.157409);
  setMatrixValue(m, 96, 0, 4.670816);  setMatrixValue(m, 96, 1, 11.266562);
  setMatrixValue(m, 97, 0, 7.109005);  setMatrixValue(m, 97, 1, 11.772769);
  setMatrixValue(m, 98, 0, 3.113531);  setMatrixValue(m, 98, 1, 60.972161);
  setMatrixValue(m, 99, 0, 3.855047);  setMatrixValue(m, 99, 1, 91.930780);

  initUIVector(&selections);

  /*MDC(m, 0, 0, &selections, 4, &run);*/
  MaxDis(m, floor(0.35*m->row), 0, &selections, 4, &run);

  puts("Selections");

  for(i = 0; i < selections->size; i++){
    printf("%f\t%f\n", getMatrixValue(m, getUIVectorValue(selections, i), 0), getMatrixValue(m, getUIVectorValue(selections, i), 1));
  }

  DelUIVector(&selections);
  DelMatrix(&m);
}

int main(void){
  /* Selection Tests
  test1();*/
  //test2();
  test3();
  /*test4();*/

  /*Clustering Tests*/
  // test5();
  /*test6();
  test7();
  test8();
  test9();
  test10();*/

  //test11();
  /*test12();*/

  /*test13();

  test14();
  test15();*/
   // test16(); HERE LAST TIME
  /*test17();*/
  //test18();
  //test19();
  return 0;
}

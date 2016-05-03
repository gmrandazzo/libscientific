/* clustering.c
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
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "scientificinfo.h"
#include "clustering.h"
#include "metricspace.h"
#include "numeric.h"
#include "matrix.h"
#include "memwrapper.h"

void MDC(matrix* m, size_t n, int metric, uivector** selections, ssignal *s)
{
  size_t i, j, k, l, mdc, nmdc;
  double d, dist, d_a, d_b;
  dvector *vectinfo, *rankvector;
  matrix *tmprank;

  NewDVector(&vectinfo, m->row);
  NewDVector(&rankvector, m->row);
  DVectorSet(vectinfo, 0.f);
  NewMatrix(&tmprank, m->row, 2);

  for(i = 0; i < m->row; i++){
    for(k = 0; k < m->row; k++){
      dist = 0.f;
        if(metric == 0){ /* EUCLIDEAN DISTANCE */
        for(j = 0; j < m->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          dist += square(m->data[i][j] - m->data[k][j]);
        }
        tmprank->data[k][0] = sqrt(dist); tmprank->data[k][1] = k;
      }
      else if(metric == 1){ /*Manhattan DISTANCE*/
        for(j = 0; j < m->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          dist += fabs(m->data[i][j] - m->data[k][j]);
        }
        tmprank->data[k][0] = sqrt(dist); tmprank->data[k][1] = k;
      }
      else{ /* COSINE DISTANCE */
        d_a = d_b = 0.f;
        for(j = 0; j < m->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          dist += (m->data[i][j] * m->data[k][j]);
          d_a += square(m->data[i][j]);
          d_b += square(m->data[k][j]);
        }
        tmprank->data[k][0] = dist/(sqrt(d_a)*sqrt(d_b)); tmprank->data[k][1] = k;
      }
    }

    MatrixSort(tmprank, 0);

    /*
    PrintMatrix(tmprank);
    */
    /*
    * Calculate the reciprocal of the rank
    */
    d = 2;
    for(k = 0; k < tmprank->row; k++){
      if(k == i){
        vectinfo->data[k] += 1;
      }
      else{
        j = tmprank->data[k][1];
        vectinfo->data[j] += 1/d;
        d += 1;
      }
    }
  }

  /*
  puts("vectinifo");
  PrintDVector(vectinfo);
  sleep(2);
  */

  nmdc = 0;
  while(1){
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      /*puts("---------------------- Start Selection --------------------"); */
      /* Find the compound with largest value in vectinfo */

      dist = vectinfo->data[0];
      mdc = 0;

      /*
      printf("init dist[0] %f\n", dist);
      */

      for(i = 0; i < vectinfo->size; i++){
        if(vectinfo->data[i] > dist){
          dist = vectinfo->data[i];
          mdc = i;
        }
        else{
          continue;
        }
      }

      /*
      printf("max dist[%lu] %f\n", (size_t)mdc, dist);
      */

      nmdc++;

      /*now mdc is the MDC */

      UIVectorAppend(selections, mdc);
      /*Recalculate the disances of the MDC to all the other compounds and the reciprocal ranks as aboce subtract these reciprocal ranks from 1 and store in the rank vector R*/

      for(k = 0; k < m->row; k++){
        dist = 0.f;

        if(metric == 0){ /* EUCLIDEAN DISTANCE */
          for(j = 0; j < m->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
            dist += square(m->data[mdc][j] - m->data[k][j]);
          }
          tmprank->data[k][0] = sqrt(dist); tmprank->data[k][1] = k;
        }
        else if(metric == 1){ /*Manhattan DISTANCE*/
          for(j = 0; j < m->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
            dist += fabs(m->data[mdc][j] - m->data[k][j]);
          }
          tmprank->data[k][0] = sqrt(dist); tmprank->data[k][1] = k;
        }
        else{ /* COSINE DISTANCE */
          d_a = d_b = 0.f;
          for(j = 0; j < m->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
            dist += m->data[mdc][j] * m->data[k][j];
            d_a += square(m->data[mdc][j]);
            d_b += square(m->data[k][j]);
          }
          tmprank->data[k][0] = dist/(sqrt(d_a)*sqrt(d_b)); tmprank->data[k][1] = k;
        }

//         setDVectorValue(tmp, k, sqrt(dist));
      }

      /*
      puts("tmpvector wich is rankvector");
      PrintDVector(tmp);
      */

      MatrixSort(tmprank, 0);

      d = 2;
      for(i = 0; i < tmprank->row; i++){
        j = tmprank->data[i][1];
        if(j == mdc){
          rankvector->data[j] = 0.f;/* set to 0 the mdc in order to unselect this..*/
        }
        else{
          rankvector->data[j] = 1 - (1 / d);
          d += 1;
        }
      }

      /*
      puts("rankvector");
      PrintDVector(rankvector);
      */

      /* Multiply values in I by the corresponding values in  R. Store the result in I.*/
      for(i = 0 ; i < vectinfo->size; i++){
        vectinfo->data[i] *= rankvector->data[i];
      }

      /*check that all number in the I Exceded 1 if they do then go to "Find the compound with largest value in vectinfo"
      * else stop....
      */

      /*
      puts("vectinfo * rankvector");
      PrintDVector(vectinfo);
      sleep(1);
      */
      if(n > 0){
        if(nmdc < n){
          continue;
        }
        else{
          break;
        }
      }
      else{
        k = 0; // objects > 1
        l = 0; // extracted objects
        for(i = 0; i < vectinfo->size; i++){
          double val = getDVectorValue(vectinfo, i);
          if(val > 1){
            k++;
          }
          else{
            if(FLOAT_EQ(val, 0.f, EPSILON)){
              l++;
            }
            else{
              continue;
            }
          }
        }

        if(k == vectinfo->size-l){
          continue;
        }
        else{
          break;
        }
      }
    }
  }


  DelMatrix(&tmprank);
  DelDVector(&rankvector);
  DelDVector(&vectinfo);
}

void MaxDis(matrix* m, size_t n, int metric, uivector** selections, ssignal *s)
{
  size_t i, j, l, nobj, ntotobj;
  double dis;
  matrix *m1, *m2, *distances;
  uivector *idselection/*, *discounter*/;
  dvector *tmp, *mindists;

  initMatrix(&m2);

  initUIVector(&idselection);
  MDC(m, 1, metric, &idselection, s);

  /* 1. Initialise Subset by transferring to it a componund */
  if(s != NULL && (*s) == SIGSCIENTIFICRUN){
    UIVectorAppend(selections, getUIVectorValue(idselection, 0));

    tmp = getMatrixRow(m, getUIVectorValue(idselection, 0));
    MatrixAppendRow(&m2, tmp);
    DelDVector(&tmp);

    DelUIVector(&idselection);

    /* 2. Calculate the dissimilarity between each remaining object in Database and the compounds in Subset */
    if(n > m->row){
      ntotobj = m->row;
    }
    else{
      ntotobj = n;
    }

    for(nobj = 1; nobj < ntotobj; nobj++){
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        initMatrix(&m1);
        initUIVector(&idselection);

        for(i = 0; i < m->row; i++){
          if(UIVectorHasValue((*selections), i) == 1){
            tmp = getMatrixRow(m, i);
            MatrixAppendRow(&m1, tmp);
            UIVectorAppend(&idselection, i);
            DelDVector(&tmp);
          }
          else{
            continue;
          }
        }

        initMatrix(&distances);

        if(metric == 0){
          EuclideanDistance(m1, m2, &distances);
        }
        else if(metric == 1){
          ManhattanDistance(m1, m2, &distances);
        }
        else{
          CosineDistance(m1, m2, &distances);
        }

        /*
        puts("--------------------------------------------");
        PrintMatrix(distances);
        */

        /* 3. The next object to be selected is always as distant as possible from already selected molecules */

        NewDVector(&mindists, m1->row);

        /* Select the minumum distances from all distances */
        for(j = 0; j < distances->col; j++){ /*for each molecule remaining in database */
          dis = getMatrixValue(distances, 0, j);
          for(i = 1; i < distances->row; i++){
            if(getMatrixValue(distances, i, j) < dis){
              dis = getMatrixValue(distances, i, j);
            }
            else{
              continue;
            }
          }
          setDVectorValue(mindists, j, dis);
        }


        /*Select the maximum object distant from all minimum distances */

        l = 0;
        for(i = 1; i < mindists->size; i++){
          if(getDVectorValue(mindists, i) > getDVectorValue(mindists, l)){
            l = i;
          }
          else{
            continue;
          }
        }

        /* l is the max min object to select */
        UIVectorAppend(selections, getUIVectorValue(idselection, l));

        tmp = getMatrixRow(m, getUIVectorValue(idselection, l));
        MatrixAppendRow(&m2, tmp);
        DelDVector(&tmp);


        DelDVector(&mindists);
        DelMatrix(&distances);
        DelUIVector(&idselection);
        DelMatrix(&m1);
      }
    }

  }
  else{
    DelUIVector(&idselection);
  }

  DelMatrix(&m2);
}

typedef struct{
  uivector *id; /* each index is the step id... */
} OBJ;


void HyperGridMap(matrix* m, size_t step, size_t n, uivector** selections, ssignal *s)
{
  size_t i, j, k, step_;
  matrix *maxminstep, *membership;
  uivector *population;

  OBJ *objs = xmalloc(sizeof(OBJ)*m->row);

  NewMatrix(&maxminstep, m->col, 3);


  step_ = step;

  /* get the max and min for each column.... */
  for(j = 0; j < m->col; j++){
    MatrixColumnMinMax(m, j, &maxminstep->data[j][0], &maxminstep->data[j][1]);
    setMatrixValue(maxminstep, j, 2, (getMatrixValue(maxminstep, j, 1)-getMatrixValue(maxminstep, j, 0)) / step_);
  }

  /* for each object check what is the bin membership and store in de id vector OJB.id */
  for(i = 0; i < m->row; i++){
    NewUIVector(&objs[i].id, m->col);

    for(j = 0; j < m->col; j++){
      setUIVectorValue(objs[i].id, j, (size_t) floor((getMatrixValue(m, i, j) - getMatrixValue(maxminstep, j, 0)) / getMatrixValue(maxminstep, j, 2)));
    }
  }

  for(i = 0; i < m->row; i++){
    PrintUIVector(objs[i].id); /* the value objs[i].id[j]+1 multiplied for the j step maxminstep[j][2]
    is equal to the min bin....*/
    puts("##############################");
  }

  /* Build the membership square matrix */
  NewMatrix(&membership, m->row, m->row);
  NewUIVector(&population, m->row);

  for(i = 0; i < m->row; i++){
    setMatrixValue(membership, i, i, i+1);
    for(j = i+1; j < m->row; j++){
      for(k = 0; k < objs[i].id->size; k++){
        if(getUIVectorValue(objs[i].id, k) == getUIVectorValue(objs[j].id, k)){
          continue;
        }
        else{
          break;
        }
      }

      if(k == objs[i].id->size){ /* j object in the same beans of i and so mark with 1*/
        setMatrixValue(membership, i, j, 1.f);
        setMatrixValue(membership, j, i, 1.f);
        setUIVectorValue(population, i, getUIVectorValue(population, i)+1);
      }
      else{ /* not in the same beans then mark with 0*/
        setMatrixValue(membership, i, j, 0.f);
        setMatrixValue(membership, j, i, 0.f);
      }
    }
  }

  PrintMatrix(membership);

  DelUIVector(&population);

  DelMatrix(&membership);

  for(i = 0; i < m->row; i++){
    DelUIVector(&objs[i].id);
  }

//   xfree(&objs);

  DelMatrix(&maxminstep);
}


/*
 * David Arthur KMeans++ init centers
 * Q: Can you explain how k-means++ works in more detail?

  k-means++ is just a way of choosing the initial centers for k-means. After that, we run k-means as
  normal. So, suppose we want to choose k initial centers from a point-set (x_1, x_2, ..., x_n). Here
  is the full algorithm to choose a set C of centers:

    1. Choose one point uniformly at random from (x_1, x_2, ..., x_n), and add it to C.
    2. For each point x_i, set D(x_i) to be the distance between x_i and the nearest point in C.
    3. Choose a real number y uniformly at random between 0 and D(x_1)^2 + D(x_2)^2 + ... + D(x_n)^2.
    4. Find the unique integer i so that
         D(x_1)^2 + D(x_2)^2 + ... D(x_i)^2 >= y > D(x_1)^2 + D(x_2)^2 + ... + D(x_(i-1))^2.

         se ho 1000 punti e in questo caso sto esaminando il 38 punti allora
         A =    D(x_1)^2 + D(x_2)^2 + ... D(x_38)^2
         B =  D(x_1)^2 + D(x_2)^2 + ... + D(x_37)^2
         quindi deve essere soddisfatta la relazione
         A >= y
         y > B
         in maniera tale da trovare il punto che si trova nella distriobuzione random più densa possibile.

    5. Add x_i to C.
    6. Repeat Steps 2-5 until we have chosen k centers.

Q: Why bother with randomness? Wouldn't it just be better to choose x_i with maximum D(x_i) instead?

  Absolutely not. This is a bad idea for several reasons:
    1. It performs worse in practice.
    2. Unlike k-means++, this way offers no approximation guarantees.
    3. It virtually guarantees you will pick every outlier as a center. That's bad.
    4. It is not random. If a method is random, you can run it many times and take the best clustering.
       Randomness is a GOOD thing for k-means.
*/

void KMeansppCenters(matrix *m, size_t n, uivector **selections, ssignal *s)
{
  size_t i, j, k;
  double dist, tmp, y, A, B;
  dvector *D, *D_min, *D_square;
  size_t q = n; /*get the number of clusters*/

  /* Step 1 Set Random Function
  srand(time(0));*/
  srand(m->col+m->row+n);
  UIVectorAppend(selections, rand() % m->row); /* random point from 0 to the max data points */

  /* Step 2 */
  while(q > 1){
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      initDVector(&D); /* vettore distanza di tutti i punti. Per ogni punto c'è un valore di distanza */
      for(i = 0; i < m->row; i++){
        initDVector(&D_min);
        for(k = 0; k < (*selections)->size; k++){
          dist = 0.f;
          for(j = 0; j < m->col; j++ ){
            dist += square(getMatrixValue(m, i, j) - getMatrixValue(m, getUIVectorValue((*selections), k), j));
          }
          DVectorAppend(&D_min, sqrt(dist));
        }

        /* get the min value distance of the point x_i from C */
        dist = getDVectorValue(D_min, 0);
        for(k = 1; k < D_min->size; k++){
          if(getDVectorValue(D_min, k) < dist){
            dist = getDVectorValue(D_min, k);
          }
          else{
            continue;
          }
        }

        DVectorAppend(&D, dist);
        DelDVector(&D_min);
      }

      /* Step 3 Calculate the square of distances and store in
      * a vector and in a sum (dist)
      */
      dist = 0.f;
      initDVector(&D_square);
      for(i = 0; i < D->size; i++){
        tmp = square(getDVectorValue(D, i));
        DVectorAppend(&D_square, tmp);
        dist += tmp;
      }

      y = (rand() * dist) / RAND_MAX;

      /* Step 4 */
      A = 0.f;
      B = 0.f;
      /* Step 4: choose two number A and B that satisfy these releaction:
        * A = D(x_1)^2 + D(x_2)^2 + ... D(x_i)^2
        * B = D(x_1)^2 + D(x_2)^2 + ... + D(x_(i-1))^2
        */
      for(i = 0; i < D_square->size; i++){
        for(j = 0; j <= i; j++){
          A += getDVectorValue(D_square, j);
          if(j < i){
            B += getDVectorValue(D_square, j);
          }
          else{
            continue;
          }
        }

        if(A >= y && y > B ){
          if(UIVectorHasValue((*selections), i) == 1){
            UIVectorAppend(selections, i);
            q--;
            break;
          }
          else{
            continue;
          }
        }
      }
      DelDVector(&D_square);
      DelDVector(&D);
    }
  }
}

/*This function rank and get the nmaxobj near or far from centroids */
void PruneResults(matrix *m, matrix *centroids, size_t nmaxobj, int type, uivector* clusters)
{
  size_t i, j, n, k, l;
  double var;
  matrix *submx, *subcentroid, *distmx;
  dvector *tmp;
  uivector *clusters_, *ids;

  NewUIVector(&clusters_, clusters->size);

  UIVectorSet(clusters_, 0);

  for(n = 0; n < centroids->row; n++){
    initMatrix(&subcentroid);

    tmp = getMatrixRow(centroids, n);
    MatrixAppendRow(&subcentroid, tmp);
    DelDVector(&tmp);

    initMatrix(&submx);
    initUIVector(&ids);
    for(i = 0; i < clusters->size; i++){
      if(getUIVectorValue(clusters, i) == n+1){
        tmp = getMatrixRow(m, i);
        MatrixAppendRow(&submx, tmp);
        DelDVector(&tmp);
        UIVectorAppend(&ids, i);
      }
      else{
        continue;
      }
    }

    initMatrix(&distmx);
    EuclideanDistance(subcentroid, submx, &distmx);

    if(type == 0){ /* Near Object */
      for(j = 0; j < nmaxobj; j++){
        l = 0;
        k = getUIVectorValue(ids, 0);
        var = getMatrixValue(distmx, l, 0);
        for(i = 1; i < distmx->row; i++){
          if(getMatrixValue(distmx, i, 0) < var){
            var = getMatrixValue(distmx, i, 0);
            k =  getUIVectorValue(ids, i);
            l = i;
          }
          else{
            continue;
          }
        }

        if(FLOAT_EQ(var, 99999, EPSILON)){
          setMatrixValue(distmx, l, 0, 99999);
        }
        else{
          setMatrixValue(distmx, l, 0, 99999);
          setUIVectorValue(clusters_, k, n+1);
        }

      }
    }
    else{ /*  Faar Object */
      for(j = 0; j < nmaxobj; j++){
        l = 0;
        k = getUIVectorValue(ids, 0);
        var = getMatrixValue(distmx, l, 0);
        for(i = 1; i < distmx->row; i++){
          if(getMatrixValue(distmx, i, 0) > var){
            var = getMatrixValue(distmx, i, 0);
            k =  getUIVectorValue(ids, i);
            l = i;
          }
          else{
            continue;
          }
        }

        if(FLOAT_EQ(var, 0, EPSILON)){
          setMatrixValue(distmx, l, 0, 0);
        }
        else{
          setMatrixValue(distmx, l, 0, 0);
          setUIVectorValue(clusters_, k, n+1);
        }
      }
    }

    DelMatrix(&distmx);
    DelMatrix(&subcentroid);
    DelUIVector(&ids);
    DelMatrix(&submx);
  }

  for(i = 0; i < clusters->size; i++){
    setUIVectorValue(clusters, i, getUIVectorValue(clusters_, i));
  }

  DelUIVector(&clusters_);
}


void RecalcCenctroids(matrix *m, matrix *G, matrix *centroids)
{
  size_t i, j, k, n, l, q;
  double sum;
  matrix *tmp;

  for(i = 0; i < centroids->row; i++){
    for(k = 0; k < m->col; k++){
      sum = 0.f;
      n = 0;
      for(j = 0; j< m->row; j++){
        if((int)getMatrixValue(G, i, j) == 1){
          sum += getMatrixValue(m, j, k);
          n++;
        }
        else{
          continue;
        }
      }

      if(n != 0){
        setMatrixValue(centroids, i, k, (double)sum/(double)n);
      }
      else{
        /* Remove the i centroid and remove the cluster i in the group matrix G */
        initMatrix(&tmp);
        MatrixCopy(centroids, &tmp);
        ResizeMatrix(&centroids, centroids->row-1, centroids->col);
        q = 0;
        for(l = 0; l < tmp->row; l++){
          if(l != i){
            for(j = 0; j < tmp->col; j++){
              setMatrixValue(centroids, q, j, getMatrixValue(tmp, l, j));
            }
            q++;
          }
          else{
            continue;
          }
        }

        ResizeMatrix(&tmp, G->row, G->col);

        for(l = 0; l < G->row; l++){
          for(j = 0; j < G->col; j++){
            setMatrixValue(tmp, l, j, getMatrixValue(G, l, j));
          }
        }

        ResizeMatrix(&G, G->row-1, G->col);

        q = 0;
        for(l = 0; l < tmp->row; l++){
          if(l != i){
            for(j = 0; j < tmp->col; j++){
              setMatrixValue(G, q, j, getMatrixValue(tmp, l, j));
            }
            q++;
          }
          else{
            continue;
          }
        }

        DelMatrix(&tmp);
        i = 0;
        break;
      }
    }
  }
}

/*
 * 1. Random Centroid made by kmeans++ step 1 algorithm
 * 2. Distance object-centroids
 * 3. Group based on minimum distance
 * 4. If no object moved to group stop; else recompute new centroids and go to step 2
*/

void KMeans(matrix* m, size_t nclusters, int initializer, uivector** clusters, matrix **_centroids_, ssignal *s)
{
  size_t i, j, k, iter;
  int iterate = 1;
  uivector *selections;
  matrix *G, *G_old,  *distmx, *centroids_, **centroids;

  NewMatrix(&G, nclusters, m->row); /* These are int matrix */
  NewMatrix(&G_old, nclusters, m->row); /* These are int matrix */

  if(_centroids_ == NULL){
    centroids = &centroids_;
    NewMatrix(centroids, nclusters, m->col);
  }
  else{
    centroids = _centroids_;
    ResizeMatrix(centroids, nclusters, m->col);
  }

  initUIVector(&selections);

  /* Step 1. Select start centroids */
  if(initializer == 0){ /* Random */
    for(i = 0; i < nclusters; i++){
      srand(m->col+m->row+nclusters+i);
      UIVectorAppend(&selections, rand() % m->row);
    }
  }
  else if(initializer == 1){ /* KMeansppCenters */
    KMeansppCenters(m, nclusters, &selections, s);
  }
  else if(initializer == 2){ /* MDC */
    MDC(m, nclusters, 0, &selections, s);
  }
  else{ /*if(initializer == 3){  MaxDis */
    MaxDis(m, nclusters, 0, &selections, s);
  }

  /* else personal centroid configuration */

  for(i = 0; i < selections->size; i++){
    for(j = 0; j < m->col; j++){
      setMatrixValue((*centroids), i, j, getMatrixValue(m, getUIVectorValue(selections, i), j));
    }
  }

  DelUIVector(&selections);


  #ifdef DEBUG
  puts("Start Centroids");
  PrintMatrix((*centroids));
  #endif

  iter = 0;
  while(iterate == 1){
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      /* Step 2. Distance all object from (*centroids) */
      initMatrix(&distmx);
      EuclideanDistance(m, (*centroids), &distmx);

      if((*centroids)->row != G_old->row){ /*If one cluster was excluded from the previoys recalculation of centroid then resize the Group matrix */
        ResizeMatrix(&G, G_old->row, m->row);
      }

      /* Step 3. Group based on minimum distance. G is the Group matrix mark with 1 if is inside the cluster */
      for(j = 0; j < distmx->col; j++ ){
        k = 0;
        for(i = 1; i < distmx->row; i++ ){
          if(getMatrixValue(distmx, i, j) < getMatrixValue(distmx, k, j)){
            k = i;
          }
          else{
            continue;
          }
        }
        setMatrixValue(G, k, j, 1.0);
      }

      DelMatrix(&distmx);

      #ifdef DEBUG
      puts("Group Matrix");
      PrintMatrix(G);
      #endif


      /* Step 4. If no object moved to group stop; else recompute new (*centroids) and go to step 2 */
      if(iter != 0){ /* for first iteration store onle G in G_old */
        for(i = 0; i < G->row; i++){
          for(j = 0; j < G->col; j++)
            if((int)getMatrixValue(G, i, j) != (int)getMatrixValue(G_old, i, j)){
              iterate = 1;
            }
            else{
              iterate = 0;
            }
        }
        iter++;
      }
      else{
        iter++;
      }

      /* store in G G_old for next iteration and reset G*/
      for(i = 0; i < G->row; i++){
        for(j = 0; j < G->col; j++){
          setMatrixValue(G_old, i, j, getMatrixValue(G, i, j));
          setMatrixValue(G, i, j, 0.f);
        }
      }

      #ifdef DEBUG
      puts("Old Centroid coordinate");
      PrintMatrix((*centroids));
      #endif

      /* recompute (*centroids) */
      RecalcCenctroids(m, G_old, (*centroids));

      #ifdef DEBUG
      puts("New Centroid Coordinate");
      PrintMatrix((*centroids));
      #endif
    }
  }

  if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
    DelMatrix(&G);
    DelMatrix(&G_old);
    if(_centroids_ == NULL){
      DelMatrix(centroids);
    }
  }
  else{
    /* End Kmeans and store centroids coordinate and cluster number for each row */
    /* recompute final centroids to insert in the output file */
    RecalcCenctroids(m, G_old, (*centroids));

    UIVectorResize(clusters, m->row);

    for(i = 0; i < G_old->row; i++){
      for(j = 0; j < G_old->col; j++){
        if(getMatrixValue(G_old, i, j) == 1){
          setUIVectorValue((*clusters), j, i+1);
        }
        else{
          continue;
        }
      }
    }

    DelMatrix(&G);
    DelMatrix(&G_old);
    if(_centroids_ == NULL){
      DelMatrix(centroids);
    }
  }
}

void KMeansRandomGroupsCV(matrix* m, size_t maxnclusters, int initializer, size_t groups, size_t iterations, dvector** ssdist, ssignal *s)
{

  size_t i, j, k, g, n, a, iterations_;
  double mindist;
  matrix *gid, *subm, *predm, *centroids, *distances;
  uivector *clusters;

  NewMatrix(&gid, groups, (size_t)ceil(m->row/(double)groups));

  srand(groups*m->row*iterations);

  DVectorResize(ssdist, maxnclusters);

  iterations_ = 0;
  while(iterations_ <  iterations){
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      /* Divide in group  all the Dataset */
      MatrixSet(gid, -1);

      /* step 1 generate the random groups */
      k = 0;
      for(i = 0; i <  gid->row; i++){
        for(j = 0; j <  gid->col; j++){
          do{
            n = (size_t)rand() % (m->row);
          } while(ValInMatrix(gid, n) == 1 && k < (m->row));
          if(k < m->row){
            setMatrixValue(gid, i, j, n);
            k++;
          }
          else
            continue;
        }
      }

      #ifdef DEBUG
      puts("Gid Matrix");
      PrintMatrix(gid);
      #endif

      /*step 2*/
      for(g = 0; g < gid->row; g++){ /*For aeach group */
        /* Estimate how many objects are inside the sub model without the group "g" */
        n = 0;
        for(i = 0; i < gid->row; i++){
          if(i != g){
            for(j = 0; j < gid->col; j++){
              if((int)getMatrixValue(gid, i, j) != -1)
                n++;
              else
                continue;
            }
          }
          else
            continue;
        }

        /*Allocate the submodel*/
        NewMatrix(&subm, n, m->col);

        /* Estimate how many objects are inside the group "g" to predict*/
        n = 0;
        for(j = 0; j < gid->col; j++){
          if((int)getMatrixValue(gid, g, j) != -1)
            n++;
          else
            continue;
        }


        /*Allocate the */
        NewMatrix(&predm, n, m->col);


        /* copy the submodel values */

        for(i = 0, k = 0; i < gid->row; i++){
          if(i != g){
            for(j = 0; j < gid->col; j++){
              a =  (size_t)getMatrixValue(gid, i, j); /* get the row index */
              if(a != -1){
                for(n = 0; n < m->col; n++){
                  setMatrixValue(subm, k, n, getMatrixValue(m, a, n));
                }
                k++;
              }
              else{
                continue;
              }
            }
          }
          else{
            continue;
          }
        }

        /* copy the objects to predict into predictm*/
        for(j = 0, k = 0; j < gid->col; j++){
          a = (size_t)getMatrixValue(gid, g, j);
          if(a != -1){
            for(n = 0; n < m->col; n++){
              setMatrixValue(predm, k, n, getMatrixValue(m, a, n));
            }
            k++;
          }
          else{
            continue;
          }
        }

        /* Kmeans
          * calculate the sum of square of the distance for each object from all the centroids...
          * store in the vector ssdist
        */
        for(j = 1; j <= maxnclusters; j++){
          initMatrix(&centroids);
          initUIVector(&clusters);

          KMeans(subm, j, initializer, &clusters, &centroids, s);

          initMatrix(&distances);
          EuclideanDistance(centroids, predm, &distances);

          #ifdef DEBUG
          puts("Centroids");
          PrintMatrix(centroids);
          puts("Matrix to predict");
          PrintMatrix(predm);
          puts("Distances");
          PrintMatrix(distances);
          #endif

          /*Get the minumum distance for each point ant sum it in ssdist. This is the kmeans clustering*/
          for(i = 0; i < distances->row; i++){ /* for each object */
            mindist = getMatrixValue(distances, i, 0);
            for(k = 1; k < distances->col; k++){ /*for each centroid */
              if(getMatrixValue(distances, i, k) < mindist){
                mindist = getMatrixValue(distances, i, k);
              }
              else{
                continue;
              }
            }
            setDVectorValue((*ssdist), j-1, getDVectorValue((*ssdist), j-1) + mindist);
          }
          DelMatrix(&distances);
          DelUIVector(&clusters);
          DelMatrix(&centroids);
        }


        DelMatrix(&subm);
        DelMatrix(&predm);
      }
    }
    iterations_++;
  }

  if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
    DVectorResize(ssdist, 0);
    DelMatrix(&gid);
  }
  else{
    /* divide all the value of the ssdist for the number of iteration */
    for(i = 0; i < (*ssdist)->size; i++){
      setDVectorValue((*ssdist), i, getDVectorValue((*ssdist), i) / iterations);
    }

    DVectNorm((*ssdist), (*ssdist));
    DelMatrix(&gid);
  }
}

/*

  Catherine A. Sugar and Gareth M. James (2003): Finding the number of Clusters in a Dataset,
  Journal of the American Statistical Association
  98:463, 750-763

  JumpMethod(X):
    Let Y = (p/2)
    Init a list D, of size n+1
    Let D[0] = 0
    For k = 1 ... n:
        Cluster X with k clusters (e.g., with k-means)
        Let d = Distortion of the resulting clustering
        D[k] = d^(-Y)
    Define J(i) = D[i] - D[i-1]
    Return the k between 1 and n that maximizes J(k)
*/
void KMeansJumpMethod(matrix* m, size_t maxnclusters, int initializer, dvector** jumps, ssignal *s)
{
  size_t k, i;
  double dist,  y, min, max;
  matrix *centroids, *xcenter;
  uivector *clusters;
  dvector *d;

  DVectorResize(jumps, maxnclusters);

  NewMatrix(&xcenter, 1, m->col);

  NewDVector(&d, maxnclusters+1);

  y = (double)m->col/2.f;

  for(k = 2; k <= maxnclusters; k++){
    dist = 0.f;
    initUIVector(&clusters);
    initMatrix(&centroids);

    KMeans(m, k, initializer, &clusters, &centroids, s);

    dist = MahalanobisDistance(m, centroids);

    printf("dist %f  %f\n", dist, pow(dist, -y));

    if(FLOAT_EQ(dist, 0.f, EPSILON)){
      setDVectorValue(d, k, 0.f);
    }
    else{
      setDVectorValue(d, k, pow(dist, -y));
    }

    DelUIVector(&clusters);
    DelMatrix(&centroids);
  }

  for(i = 1; i < d->size; i++){
    setDVectorValue((*jumps), i-1, getDVectorValue(d, i)-getDVectorValue(d, i-1));
  }

  DVectorMinMax((*jumps), &min, &max);

  y = max - min;
  for(i = 0; i < (*jumps)->size; i++){
    setDVectorValue((*jumps), i, (getDVectorValue((*jumps), i)- min) / y);
  }

  DelDVector(&d);
  DelMatrix(&xcenter);
}

/*hierarchical clustering with different linkage criterion*/
void HierarchicalClustering(matrix* _m, size_t nclusters, uivector** _clusters, matrix **_centroids_, strvector **dendogram, enum LinkageType linktype, ssignal *s)
{
  size_t i, j, k, l, m, min_i, min_j;
  char buffer[MAXCHARSIZE];
  /* calc the matrix distance */
  matrix *distmx, *distmx_new;
  strvector *pointname, *pointname_new, *clusters, *tokens;
  dvector *clusterdist, *tmp;

  initMatrix(&distmx);
  initStrVector(&pointname);
  initStrVector(&clusters);
  initDVector(&clusterdist);

  if(linktype > 2){
    SquaredEuclideanDistance(_m, _m, &distmx);
  }
  else{
    EuclideanDistance(_m, _m, &distmx);
  }

  for(i = 0; i < _m->row; i++){
    StrVectorAppendInt(&pointname, i);
  }


  m = 1; /* identifier for cluster */

  while(m < _m->row){
    /* Step 2. find the minimum distance value inside the distance matrix */
    min_i = 0;
    min_j = 1;
    for(i = 0; i < distmx->row; i++ ){
      for(j = i+1; j < distmx->col; j++ ){
        if(getMatrixValue(distmx, i, j) < getMatrixValue(distmx, min_i, min_j) &&
          !FLOAT_EQ(getMatrixValue(distmx, i, j), 0.f, EPSILON)){
          min_i = i;
          min_j = j;
        }
        else{
          continue;
        }
      }
    }

    /* Step 3 merge the cluster and increase m. */
    DVectorAppend(&clusterdist, getMatrixValue(distmx, min_i, min_j));
    strcpy(buffer, "");
    strcpy(buffer, getStr(pointname, min_i));
    strcat(buffer, ";");
    strcat(buffer, getStr(pointname, min_j));
    strcat(buffer, ";");
    StrVectorAppend(&clusters, buffer);
//     printf("clusters: %s == %s?\n", getStr(clusters, clusters->size-1), buffer);

    strcpy(buffer, "");
    strcpy(buffer, getStr(pointname, min_i));
    strcat(buffer, ";");
    strcat(buffer, getStr(pointname, min_j));
    setStr(pointname, min_i, buffer);

    strcpy(buffer, "");
    strcpy(buffer, getStr(pointname, min_j));
    strcat(buffer, ";");
    strcat(buffer, getStr(pointname, min_i));
    setStr(pointname, min_j, buffer);
    m++;

    /*Step 4 Remove the min_i and min_j from the matrix and update the distance matrix
    * Erase Row with index "min_i"
    */

    initDVector(&tmp);
    for(i = 0; i < distmx->row; i++){
      if(i == min_i){
        DVectorAppend(&tmp, 0.f);
      }
      else if(i == min_j){
        continue;
      }
      else{
        if(linktype == 0){
          /* Single-Linkage Criterion */
        /* printf("ok: %f\t%f\n", getMatrixValue(distmx, min_i, i), getMatrixValue(distmx, min_j, i)); */
          if( getMatrixValue(distmx, min_i, i) < getMatrixValue(distmx, min_j, i)){
            DVectorAppend(&tmp, getMatrixValue(distmx, min_i, i));
          }
          else{
            DVectorAppend(&tmp, getMatrixValue(distmx, min_j, i));
          }
        }
        else if(linktype == 1){
          /*Complete Linkage*/
          if( getMatrixValue(distmx, min_i, i) > getMatrixValue(distmx, min_j, i)){
            DVectorAppend(&tmp, getMatrixValue(distmx, min_i, i));
          }
          else{
            DVectorAppend(&tmp, getMatrixValue(distmx, min_j, i));
          }
        }
        else if(linktype == 2){
          /*Average Linkage*/
          DVectorAppend(&tmp, sqrt(square(distmx->data[min_i][i]-distmx->data[min_j][i]))/2.);
        }
        else{
          /* Ward-Linkage Criterion */
          DVectorAppend(&tmp, square(distmx->data[min_i][i]-distmx->data[min_j][i]));
        }
      }
    }

    NewMatrix(&distmx_new, distmx->row-1, distmx->col-1);
    NewStrVector(&pointname_new, pointname->size-1);
    if(min_i < min_j){
      /* Erase the row min_j and col min_j */
      k = 0;
      for(i = 0; i < distmx->row; i++){
        if(i != min_j){
          setStr(pointname_new, k, getStr(pointname, i));
          l = 0;
          for(j = 0; j < distmx->col; j++){
            if(j != min_j){
              setMatrixValue(distmx_new, k, l, getMatrixValue(distmx, i, j));
              l++;
            }
            else{
              continue;
            }
          }
          k++;
        }
        else{
          continue;
        }
      }

      /*
      puts("OLD MX");
      PrintMatrix(distmx);
      printf("NEW MX without row and column %u\n", (unsigned int)min_j);
      PrintMatrix( distmx_new);
      */

      /* Fill the column min_i with the row[i] */
      for(i = 0; i < tmp->size; i++){
        setMatrixValue(distmx, i, min_i, getDVectorValue(tmp, i));
      }
      /* Fill the row min_i with the row[i] */
      for( i=0; i < tmp->size; i++ ){
        setMatrixValue(distmx, min_i, i, getDVectorValue(tmp, i));
      }
    }
    else{
      /* Erase the row min_i and col min_i */
      k = 0;
      for(i = 0; i < distmx->row; i++){
        if(i != min_i){
          setStr(pointname_new, k, getStr(pointname, i));
          l = 0;
          for(j = 0; j < distmx->col; j++){
            if(j != min_i){
              setMatrixValue(distmx_new, k, l, getMatrixValue(distmx, i, j));
              l++;
            }
            else{
              continue;
            }
          }
          k++;
        }
        else{
          continue;
        }
      }

      /*
      puts("OLD MX");
      PrintMatrix(distmx);
      printf("NEW MX without row and column %u\n", (unsigned int)min_j);
      PrintMatrix( distmx_new);
      */

      /* Fill the column min_j with the row[i] */
      for(i = 0; i < tmp->size; i++){
        setMatrixValue(distmx, i, min_j, getDVectorValue(tmp, i));
      }
      /* Fill the row min_j with the row[i] */
      for( i=0; i < tmp->size; i++ ){
        setMatrixValue(distmx, min_j, i, getDVectorValue(tmp, i));
      }
    }

    MatrixCopy(distmx_new, &distmx);
    StrVectorResize(&pointname, pointname_new->size);
    for(i = 0; i < pointname_new->size; i++){
      setStr(pointname, i, getStr(pointname_new, i));
    }

    DelMatrix(&distmx_new);
    DelStrVector(&pointname_new);
    DelDVector(&tmp);
  }


  #ifdef DEBUG
  puts("Dendogram Clusters");
  for(i = 0; i < clusters->size; i++ ){
    printf("%s\t%f\n", getStr(clusters, i), getDVectorValue(clusterdist, i));
  }
  #endif

  if(dendogram != NULL){
    for(i = 0; i < clusters->size; i++){
      strcpy(buffer, "");
      strcpy(buffer, getStr(clusters, i));
      char fl2str[50];
      snprintf (fl2str, sizeof(fl2str), "%f", clusterdist->data[i]);
      strcat(buffer, fl2str);
      StrVectorAppend(dendogram, buffer);
    }
  }

  /* Merge Clusters to _clusters: split and select clusters. */

  UIVectorResize(_clusters, _m->row);

  k = 1;

  for(i = clusters->size-nclusters-1; i < clusters->size; i++){
    initStrVector(&tokens);
//     printf("Cluster %s\n", getStr(clusters, i));
    SplitString(getStr(clusters, i), ";", &tokens);
    m = 0;
    for(j = 0; j < tokens->size; j++){
      l = atoi(getStr(tokens, j));
      if(getUIVectorValue((*_clusters), l) == 0){
        setUIVectorValue((*_clusters), l, k);
        m = 1;
      }
      else{
        continue;
      }
    }
    DelStrVector(&tokens);

    /* if are foundet new element added to the cluster than encrease the cluster level. */
    if(m == 1){
      k++;
    }
    else{
      continue;
    }
  }

  DelMatrix(&distmx);
  DelStrVector(&pointname);
  DelStrVector(&clusters);
  DelDVector(&clusterdist);
}
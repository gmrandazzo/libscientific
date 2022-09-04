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
#include <pthread.h>

//#ifdef DEBUG
#include <time.h>
//#endif

#include "scientificinfo.h"
#include "clustering.h"
#include "metricspace.h"
#include "numeric.h"
#include "matrix.h"
#include "memwrapper.h"

typedef struct{
  matrix *m;
  matrix *tmprank;
  size_t from, to;
  size_t mdc, metric;
} mdc_th_args;

void *MDCWorker(void *arg_)
{
  size_t k, j;
  double dist;
  mdc_th_args *arg = (mdc_th_args*) arg_;

  for(k = arg->from; k < arg->to; k++){
    dist = 0.f;
    // EUCLIDEAN DISTANCE
    if(arg->metric == 0){
      for(j = 0; j < arg->m->col; j++){
        dist += square(arg->m->data[arg->mdc][j] - arg->m->data[k][j]);
      }
      dist = sqrt(dist);
    }
    // MANHATAN DISTANCE
    else if(arg->metric == 1){
      for(j = 0; j < arg->m->col; j++){
        dist += fabs(arg->m->data[arg->mdc][j] - arg->m->data[k][j]);
      }
    }
    // COSINE DISTANCE
    else{
      double d_a, d_b;
      d_a = d_b = 0.f;
      for(j = 0; j < arg->m->col; j++){
        dist += arg->m->data[arg->mdc][j] * arg->m->data[k][j];
        d_a += square(arg->m->data[arg->mdc][j]);
        d_b += square(arg->m->data[k][j]);
      }
      dist = dist/(sqrt(d_a)*sqrt(d_b));
    }
    arg->tmprank->data[k][0] = dist;
    arg->tmprank->data[k][1] = k;
  }
  return 0;
}

/*
 *
 */
void MDC(matrix* m, size_t n, int metric, uivector** selections, size_t nthreads, ssignal *s)
{
  size_t i, j, k, l, mdc, nmdc, th;
  double d, dist;
  dvector *vectinfo, *rankvector;
  matrix *tmprank;

  pthread_t *threads;
  mdc_th_args *args;

  NewDVector(&vectinfo, m->row);
  NewDVector(&rankvector, m->row);
  DVectorSet(vectinfo, 0.f);
  NewMatrix(&tmprank, m->row, 2);


  threads = xmalloc(sizeof(pthread_t)*nthreads);
  args = xmalloc(sizeof(mdc_th_args)*nthreads);
  /* reference some memory data */
  for(th = 0; th < nthreads; th++){
    args[th].m = m;
    args[th].tmprank = tmprank;
    args[th].metric = metric;
  }

  /* MULTI THREAD IMPLEMENTATION */
  matrix *dm;
  NewMatrix(&dm, m->row, m->row);
  if(metric == 0){
    EuclideanDistance(m, m, &dm, nthreads);
  }
  else if(metric == 1){
    ManhattanDistance(m, m, &dm, nthreads);
  }
  else{
    CosineDistance(m, m, &dm, nthreads);
  }

  for(i = 0; i < m->row; i++){
    for(k = 0; k < m->row; k++){
      tmprank->data[k][0] = dm->data[k][i];
      tmprank->data[k][1] = k;
    }

    MatrixSort(tmprank, 0);

    //Calculate the reciprocal of the rank
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
  DelMatrix(&dm);

  /* SINGLE THREAD IMPLEMENTATION

  for(i = 0; i < m->row; i++){
    for(k = 0; k < m->row; k++){
      dist = 0.f;
      // EUCLIDEAN DISTANCE
      if(metric == 0){
        for(j = 0; j < m->col; j++){
          dist += square(m->data[i][j] - m->data[k][j]);
        }
        tmprank->data[k][0] = sqrt(dist); tmprank->data[k][1] = k;
      }
      //Manhattan DISTANCE
      else if(metric == 1){
        for(j = 0; j < m->col; j++){
          dist += fabs(m->data[i][j] - m->data[k][j]);
        }
        tmprank->data[k][0] = sqrt(dist); tmprank->data[k][1] = k;
      }
      //  COSINE DISTANCE
      else{
        double d_a, d_b;
        d_a = d_b = 0.f;
        for(j = 0; j < m->col; j++){
          dist += (m->data[i][j] * m->data[k][j]);
          d_a += square(m->data[i][j]);
          d_b += square(m->data[k][j]);
        }
        tmprank->data[k][0] = dist/(sqrt(d_a)*sqrt(d_b)); tmprank->data[k][1] = k;
      }
    }

    MatrixSort(tmprank, 0);

    // PrintMatrix(tmprank);


    // Calculate the reciprocal of the rank
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
  */


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
      /* MULTITHREAD IMPLEMENTATION */
      size_t step = (size_t)ceil((double)m->row/(double)nthreads);
      size_t from = 0, to = step;
      for(th = 0; th < nthreads; th++){
        args[th].from = from;
        args[th].to = to;
        args[th].mdc = mdc;
        pthread_create(&threads[th], NULL, MDCWorker, (void*) &args[th]);

        from = to;
        if(from+step > m->row){
          to = m->row;
        }
        else{
          to+=step;
        }
      }

      for(th = 0; th < nthreads; th++){
        pthread_join(threads[th], NULL);
      }

      /* SINGLE TRHEAD IMPLEMENTATION
      for(k = 0; k < m->row; k++){
        dist = 0.f;
        // EUCLIDEAN DISTANCE
        if(metric == 0){
          for(j = 0; j < m->col; j++){
            dist += square(m->data[mdc][j] - m->data[k][j]);
          }
          dist = sqrt(dist);
        }
        // MANHATAN DISTANCE
        else if(metric == 1){
          for(j = 0; j < m->col; j++){
            dist += fabs(m->data[mdc][j] - m->data[k][j]);
          }
        }
        // COSINE DISTANCE
        else{
          double d_a, d_b;
          d_a = d_b = 0.f;
          for(j = 0; j < m->col; j++){
            dist += m->data[mdc][j] * m->data[k][j];
            d_a += square(m->data[mdc][j]);
            d_b += square(m->data[k][j]);
          }
          dist = dist/(sqrt(d_a)*sqrt(d_b)); tmprank->data[k][1] = k;
        }
        tmprank->data[k][0] = dist; tmprank->data[k][1] = k;
      }
      */

      /*
      puts("tmpvector wich is rankvector");
      PrintDVector(tmp);
      */
      MatrixSort(tmprank, 0);


      d = 2.;
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

  xfree(threads);
  xfree(args);
  DelMatrix(&tmprank);
  DelDVector(&rankvector);
  DelDVector(&vectinfo);
}

/*
 * MaxDis object selection.
 */
void MaxDis(matrix* m, size_t n, int metric, uivector** selections, size_t nthreads, ssignal *s)
{
  size_t i, j, l, nobj, ntotobj;
  int far_away;
  double far;
  double dis;
  matrix *m1, *m2, *distances;
  uivector *idselection/*, *discounter*/;
  dvector *tmp, *mindists;

  initMatrix(&m2);

  initUIVector(&idselection);

  /* select the faraway compound from centroid */
  dvector *c;
  NewDVector(&c, m->col);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      c->data[j] += m->data[i][j];
    }
  }

  for(j = 0; j < m->col; j++){
    c->data[j] /= (double)m->row;
  }

  far_away = -1;
  far = 0.f;
  for(i = 0; i < m->row; i++){
    double dst = 0.f;
    for(j = 0; j < m->col; j++){
      dst += square(c->data[j] - m->data[i][j]);
    }
    dst = sqrt(dst);
    if(far_away > -1){
      far = dst;
    }
    else{
      if(dst > far){
        far = dst;
        far_away = i;
      }
    }
  }
  DelDVector(&c);

  UIVectorAppend(&idselection, far_away);

  /* OLD IMPLEMENTATION SELECTED THE FIRST COMPOUND AS MDC
   * MDC(m, 1, metric, &idselection, nthreads, s);
   */

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
          EuclideanDistance(m1, m2, &distances, nthreads);
        }
        else if(metric == 1){
          ManhattanDistance(m1, m2, &distances, nthreads);
        }
        else{
          printf("Ciao\n");
          CosineDistance(m1, m2, &distances, nthreads);
          PrintMatrix(distances);
        }

        /*
        puts("--------------------------------------------");
        PrintMatrix(distances);
        */

        /* 3. The next object to be selected is always as distant as possible from already selected molecules */

        NewDVector(&mindists, m1->row);

        /* Select the minumum distances from all distances */
        for(j = 0; j < distances->col; j++){ /*for each molecule remaining in database */
          dis = distances->data[0][j];
          for(i = 1; i < distances->row; i++){
            if(distances->data[i][j] < dis){
              dis = distances->data[i][j];
            }
            else{
              continue;
            }
          }
          mindists->data[j] = dis;
        }


        /*Select the maximum object distant from all minimum distances */

        l = 0;
        for(i = 1; i < mindists->size; i++){
          if(mindists->data[i] > mindists->data[l]){
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

/*
 * Fast implementation but can pose some memory problems with large datasets!
 */
void MaxDis_Fast(matrix* m, size_t n, int metric, uivector** selections, size_t nthreads, ssignal *s)
{
  size_t i, j, indx, nobj;
  double dis;
  dvector *distances;
  uivector *id;
  dvector *c, *mindists;

  NewUIVector(&id, m->row);
  /* Store into the id array the original point positions */
  for(i = 0; i < m->row; i++)
    id->data[i] = i;

  initDVector(&distances);
  /*
   * Calculate a square distance matrix
   * Slow process!
   */
  if(metric == 0){
    EuclideanDistanceCondensed(m, &distances, nthreads);
  }
  else if(metric == 1){
    ManhattanDistanceCondensed(m, &distances, nthreads);
  }
  else{
    CosineDistanceCondensed(m, &distances, nthreads);
  }

  /* select the faraway compound from centroid */
  NewDVector(&c, m->col);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      c->data[j] += m->data[i][j];
    }
  }

  for(j = 0; j < m->col; j++){
    c->data[j] /= (double)m->row;
  }

  int far_away = -1;
  double far = 0.f;
  for(i = 0; i < m->row; i++){
    double dst = 0.f;
    for(j = 0; j < m->col; j++){
      dst += square(c->data[j] - m->data[i][j]);
    }

    dst = sqrt(dst);

    if(far_away > -1){
      far = dst;
    }
    else{
      if(dst > far){
        far = dst;
        far_away = i;
      }
      else{
        continue;
      }
    }
  }
  DelDVector(&c);

  /* Append the far away compound to the final selection */
  UIVectorAppend(selections, far_away);

  /* Remove the selected point from the position list */
  UIVectorRemoveAt(&id, far_away);

  /*
   * The next object to be selected is always as distant as possible
   * from already selected objects. Hence iterate in the distance matrix
   * using the remaining ids.
   */

  /* ntob = 1 because we have already selected the first object, the far away objcet*/
  for(nobj = 1; nobj < n; nobj++){
    if(s != NULL && (*s) == SIGSCIENTIFICRUN){
      /* Select the minumum distance of all remaining objects
       * from the already selected points
       */
      NewDVector(&mindists, id->size);
      for(i = 0; i < id->size; i++){
        size_t ii = id->data[i];
        size_t jj = (*selections)->data[0];
        indx = square_to_condensed_index(ii, jj, m->row);
        dis = distances->data[indx];
        for(j = 1; j < (*selections)->size; j++){
          jj = (*selections)->data[j];
          indx = square_to_condensed_index(ii, jj, m->row);
          if(distances->data[indx] < dis){
            dis = distances->data[indx];
          }
          else{
            continue;
          }
        }
        mindists->data[i] = dis;
      }

      /*
       * From the final smallest distance list select the maximum distant objects
       */
      j = 0;
      for(i = 1; i < mindists->size; i++){
        if(mindists->data[i] > mindists->data[j]){
          j = i;
        }
        else{
          continue;
        }
      }

      /* l is the max min object to select */
      UIVectorAppend(selections, id->data[j]);
      UIVectorRemoveAt(&id, j);
    }
    else{
      break;
    }
    DelDVector(&mindists);
  }
  DelDVector(&distances);
  DelUIVector(&id);
}

void NewHyperGridMap(HyperGridModel **hgm)
{
  (*hgm) = xmalloc(sizeof(HyperGridModel));
  initMatrix(&((*hgm)->gmap));
  initDVector(&((*hgm)->colaverage));
  initDVector(&((*hgm)->colscaling));
  (*hgm)->gsize = 0;
  (*hgm)->bsize = 0;

}

void DelHyperGridMap(HyperGridModel **hgm)
{
  DelMatrix(&((*hgm)->gmap));
  DelDVector(&((*hgm)->colaverage));
  DelDVector(&((*hgm)->colscaling));
  xfree((*hgm));
}

void HyperGridMap(matrix* m, size_t grid_size, hgmbins** bins_id, HyperGridModel **hgm)
{
  size_t i, j;
  matrix *gmap = (*hgm)->gmap;
  dvector *colaverage = (*hgm)->colaverage;
  dvector *colscaling = (*hgm)->colscaling;
  (*hgm)->gsize = grid_size;
  (*hgm)->bsize = 1.f;

  /*Allocate a matrix of min, max, step_size*/
  ResizeMatrix(&gmap, m->col, 3);

  /* Centering and Scaling to unit variance */
  MatrixColAverage(m, &colaverage);
  MatrixColSDEV(m, &colscaling);

  for(j = 0; j < m->col; j++){
    for(i = 0; i < m->row; i++){
      if(FLOAT_EQ(colscaling->data[j], 0.f, EPSILON)){
        m->data[i][j] = 0.f;
      }
      else{
        m->data[i][j] = (m->data[i][j]-colaverage->data[j])/colscaling->data[j];
      }
    }
  }

  /* get the max and min for each column.... */
  for(j = 0; j < m->col; j++){
    MatrixColumnMinMax(m, j, &gmap->data[j][0], &gmap->data[j][1]);
    gmap->data[j][2] = (gmap->data[j][1]-gmap->data[j][0])/(double)grid_size;
    (*hgm)->bsize *= (double)grid_size;
  }

  /*Create two vectors: one for the multiplier, the other for the index id.
   * The formula to get the bin id of each point is the following:
   *  id1*mult1 + id2*mult2 + ... + idN*multN = bin ID
   * if grid_size is 4 mult1 = 4, mult2 = 16, mult3 = 64 ...
   * This approach is ok for small dataset. However for large dataset
   * is better to have an hash long as the number of features for each bin like:
   * 135712736.
   */

  /*DVectorResize(&mult, m->col);

  double mult_ = (double)grid_size;

  mult->data[0] = 1.f;
  for(j = 1; j < m->col; j++){
    mult->data[j] = mult_;
    mult_ *= (double)grid_size;
  }*/

  if(bins_id != NULL){
    /* for each object check what is the bin membership and store in the id bins_id */
    HyperGridMapObjects(m, (*hgm), bins_id);
  }
}


void HyperGridMapObjects(matrix *m, HyperGridModel *hgm, hgmbins **bins_id)
{
  size_t i, j;
  matrix *gmap = hgm->gmap;

  (*bins_id) = xmalloc(sizeof(hgmbins));
  (*bins_id)->nobj = m->row;
  (*bins_id)->hash_size = m->col;
  (*bins_id)->hash = xmalloc(sizeof(size_t*)*m->row);
  for(i = 0; i < (*bins_id)->nobj; i++){
    (*bins_id)->hash[i] = xmalloc(sizeof(size_t)*m->col);
  }

  /*Create two vectors: one for the multiplier, the other for the index id.
   * The formula to get the bin id of each point is the following:
   *  id1*mult1 + id2*mult2 + ... + idN*multN = bin ID
   * if grid_size is 4 mult1 = 4, mult2 = 16, mult3 = 64 ...
   * This approach is ok for small dataset. However for large dataset
   * is better to have an hash long as the number of features for each bin like:
   * 135712736.
   */
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      /*if(FLOAT_EQ(m->data[i][j], gmap->data[j][0], EPSILON)){
        // if x is the minimum then is on 0
        (*bins_id)->hash[i][j] = 0;
      }
      else if(FLOAT_EQ(m->data[i][j], gmap->data[j][1], EPSILON)){
        // if x is the maximum then is on max of the grid position in this axis
        (*bins_id)->hash[i][j] = (hgm->gsize-1);
      }
      else{
        (*bins_id)->hash[i][j] = (size_t)floor((m->data[i][j] - gmap->data[j][0])/gmap->data[j][2]);
      }*/
      if(FLOAT_EQ(gmap->data[j][2], 0.f, EPSILON)){
        (*bins_id)->hash[i][j] = 0;
      }
      else{
        (*bins_id)->hash[i][j] = (size_t)floor((m->data[i][j] - gmap->data[j][0])/gmap->data[j][2]);
      }
    }
  }
}

void PrintHGMBins(hgmbins *bins_id)
{
  size_t i, j;
  for(i = 0; i < bins_id->nobj; i++){
    for(j = 0; j < bins_id->hash_size-1; j++){
      printf("%llu", bins_id->hash[i][j]);
    }
    printf("%llu\n", bins_id->hash[i][bins_id->hash_size-1]);
  }
}

void DelHGMBins(hgmbins **bins_id)
{
  size_t i;
  for(i = 0; i < (*bins_id)->nobj; i++){
    xfree((*bins_id)->hash[i]);
  }
  xfree((*bins_id)->hash);
  xfree(*bins_id);
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

typedef struct
{
  matrix *m;
  uivector *selections;
  dvector *D;
  size_t from, to;
} kmpp_th_args;

void *kmppDistanceWorker(void *arg_)
{
  size_t i, j, k;
  double dist;
  dvector *D_min;
  kmpp_th_args *a = (kmpp_th_args*) arg_;

  for(i = a->from; i < a->to; i++){
    initDVector(&D_min);
    for(k = 0; k < a->selections->size; k++){
      dist = 0.f;
      for(j = 0; j < a->m->col; j++){
        dist += (a->m->data[i][j] - a->m->data[a->selections->data[k]][j])*(a->m->data[i][j] - a->m->data[a->selections->data[k]][j]);
      }
      DVectorAppend(&D_min, sqrt(dist));
    }

    /* get the min value distance of the point x_i from C */
    dist = D_min->data[0];
    for(k = 1; k < D_min->size; k++){
      if(D_min->data[k] < dist){
        dist = D_min->data[k];
      }
      else{
        continue;
      }
    }

    a->D->data[i] = dist;
    DelDVector(&D_min);
  }

  return 0;
}

void KMeansppCenters(matrix *m, size_t n, uivector **selections, int nthreads, ssignal *s)
{
  size_t i, j;
  double dist, tmp, y, A, B;
  dvector *D, *D_square;
  size_t q = n; /*get the number of clusters*/

  pthread_t *threads = xmalloc(sizeof(pthread_t)*nthreads);
  kmpp_th_args *arg = xmalloc(sizeof(kmpp_th_args)*nthreads);

  NewDVector(&D, m->row); /* vettore distanza di tutti i punti. Per ogni punto c'è un valore di distanza */
  NewDVector(&D_square, m->row); /* vettore distanza di tutti i punti al quadrato.*/
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
      /*
      SINGLE THREAD
      for(i = 0; i < m->row; i++){
        initDVector(&D_min);
        for(k = 0; k < (*selections)->size; k++){
          dist = 0.f;
          for(j = 0; j < m->col; j++ ){
            dist += (m->data[i][j] - m->data[(*selections)->data[k]][j])*(m->data[i][j] - m->data[(*selections)->data[k]][j]);
          }
          DVectorAppend(&D_min, sqrt(dist));
        }


        //get the min value distance of the point x_i from C
        dist = D_min->data[0];
        for(k = 1; k < D_min->size; k++){
          if(D_min->data[k] < dist){
            dist = D_min->data[k];
          }
          else{
            continue;
          }
        }

        D->data[i] = dist;
        //DVectorAppend(&D, dist);
        DelDVector(&D_min);
      }
      */
      size_t nobj = ceil(m->row/(double)nthreads);
      size_t from = 0;
      for(i = 0; i < nthreads; i++){
        arg[i].m = m;
        arg[i].selections = (*selections);
        arg[i].D = D;
        arg[i].from = from;
        if(from+nobj > m->row){
          from = m->row;
        }
        else{
          from += nobj;
        }
        arg[i].to = from;
        pthread_create(&threads[i], NULL, kmppDistanceWorker, (void*) &arg[i]);
      }

      for(i = 0; i < nthreads; i++){
        pthread_join(threads[i], NULL);
      }

      /* Step 3 Calculate the square of distances and store in
      * a vector and in a sum (dist)
      */


      dist = 0.f;
      for(i = 0; i < D->size; i++){
        tmp = square(D->data[i]);
        D_square->data[i] = tmp;
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

    }
  }
  DelDVector(&D);
  DelDVector(&D_square);
  xfree(threads);
  xfree(arg);
}

/*This function rank and get the nmaxobj near or far from centroids */
void PruneResults(matrix *m, matrix *centroids, size_t nmaxobj, int type, uivector* clusters, size_t nthreads)
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
    EuclideanDistance(subcentroid, submx, &distmx, nthreads);

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

        if(FLOAT_EQ(var, MISSING, EPSILON)){
          setMatrixValue(distmx, l, 0, MISSING);
        }
        else{
          setMatrixValue(distmx, l, 0, MISSING);
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

/*
 * 1. Random Centroid made by kmeans++ step 1 algorithm
 * 2. Distance object-centroids
 * 3. Group based on minimum distance
 * 4. If no object moved to group stop; else recompute new centroids and go to step 2
*/

int shouldStop(matrix *centroids, matrix *oldcentroids, size_t iterations, size_t max_iterations)
{
  size_t i, j;
  if(iterations > max_iterations){
    return 1;
  }
  else{
    if(centroids->row != oldcentroids->row){
      return 0;
    }
    else{
      for(i = 0; i < centroids->row; i++){
        for(j = 0; j < centroids->col; j++){
          if(FLOAT_EQ(centroids->data[i][j], oldcentroids->data[i][j], EPSILON)){
            continue;
          }
          else{
            return 0;
          }
        }
      }
      return 1;
    }
  }
}

typedef struct{
  matrix *m;
  matrix *centroids;
  int from, to;
  uivector *labels;
} labels_th_arg;

void *getLabelsWorker(void *arg_)
{
  labels_th_arg *a = (labels_th_arg*) arg_;
  size_t i, j, k;
  int c_point; /*closest centroid id*/
  double c_dst; /* closest centroid distance */
  /* parallelize this part! */
  for(i = a->from; i < a->to; i++){
    c_point = -1;
    for(k = 0; k < a->centroids->row; k++){
      double t_dst = 0.f;
      for(j = 0; j < a->m->col; j++){
        t_dst += (a->m->data[i][j]-a->centroids->data[k][j])*(a->m->data[i][j]-a->centroids->data[k][j]);
      }
      t_dst = sqrt(t_dst);
      if(c_point == -1){
        c_point = 0;
        c_dst = t_dst;
      }
      else{
        if(t_dst < c_dst){
          c_point = k;
          c_dst = t_dst;
        }
        else{
          continue;
        }
      }
    }
    a->labels->data[i] = c_point;
  }
  return 0;
}

/* For each element in the dataset, chose the closest centroid.
 * Make that centroid the element's label.
 */
void getLabels_(matrix *m, matrix *centroids, uivector *labels, int nthreads)
{
  size_t i, from;
  pthread_t *threads = xmalloc(sizeof(pthread_t)*nthreads);
  labels_th_arg *arg = xmalloc(sizeof(labels_th_arg)*nthreads);
  int nobj = ceil(m->row/(double)nthreads);
  from = 0;
  for(i = 0; i < nthreads; i++){
    arg[i].m = m;
    arg[i].centroids = centroids;
    arg[i].labels = labels;
    arg[i].from = from;
    if(from+nobj > m->row){
      from = m->row;
    }
    else{
      from += nobj;
    }
    arg[i].to = from;
    pthread_create(&threads[i], NULL, getLabelsWorker, (void*) &arg[i]);
  }

  for(i = 0; i < nthreads; i++){
    pthread_join(threads[i], NULL);
  }

  xfree(threads);
  xfree(arg);
}

/* For each element in the dataset, chose the closest centroid.
 * Make that centroid the element's label.
 * SINGLE THREAD
 */
void getLabels(matrix *m, matrix *centroids, uivector *labels)
{
  size_t i, j, k;
  int c_point; /*closest centroid id*/
  double c_dst; /* closest centroid distance */
  /* parallelize this part! */
  for(i = 0; i < m->row; i++){
    c_point = -1;
    for(k = 0; k < centroids->row; k++){
      double t_dst = 0.f;
      for(j = 0; j < m->col; j++){
        t_dst += (m->data[i][j]-centroids->data[k][j])*(m->data[i][j]-centroids->data[k][j]);
      }
      t_dst = sqrt(t_dst);
      if(c_point == -1){
        c_point = 0;
        c_dst = t_dst;
      }
      else{
        if(t_dst < c_dst){
          c_point = k;
          c_dst = t_dst;
        }
        else{
          continue;
        }
      }
    }
    labels->data[i] = c_point;
  }
}

/* Calculate centroids
 * Each centroid is the geometric mean of the points that
 * have that centroid's label. Important: If a centroid is empty (no points have
 * that centroid's label) you should randomly re-initialize it.
 */
void getCentroids(matrix *m, uivector *cluster_labels, matrix **centroids)
{
  size_t i, j, c_indx;
  matrix *new_centroids;
  uivector *cluster_points;
  NewMatrix(&new_centroids, (*centroids)->row, (*centroids)->col);
  NewUIVector(&cluster_points, (*centroids)->row);
  for(i = 0; i < cluster_labels->size; i++){
    c_indx = cluster_labels->data[i];
    for(j = 0; j < m->col; j++){
      new_centroids->data[c_indx][j] += m->data[i][j];
    }
    cluster_points->data[c_indx]+=1;
  }

  /*Check if there are 0 points in one cluster.
    If yes then pick a random point from mx
    else divide the centroid for the number of points on the cluster*/
  for(i = 0; i < cluster_points->size; i++){
    if(cluster_points->data[i] > 0){
      for(j = 0; j < m->col; j++){
        new_centroids->data[i][j] /= (double)cluster_points->data[i];
      }
    }
    else{
      c_indx = rand() % m->row;
      for(j = 0; j < m->col; j++){
        new_centroids->data[i][j] = m->data[c_indx][j];
      }
    }
  }

  DelUIVector(&cluster_points);
  MatrixCopy(new_centroids, centroids);
  DelMatrix(&new_centroids);
}

void KMeans(matrix* m, size_t nclusters, int initializer, uivector** cluster_labels, matrix **_centroids_, size_t nthreads, ssignal *s)
{
  size_t i, j, it;
  matrix *centroids, *oldcentroids;
  uivector *pre_centroids;

  if(_centroids_ == NULL){
    NewMatrix(&centroids, nclusters, m->col);
  }
  else{
    centroids = (*_centroids_);
    ResizeMatrix(&centroids, nclusters, m->col);
  }

  NewMatrix(&oldcentroids, centroids->row, centroids->col);

  initUIVector(&pre_centroids);

  /* Step 1. Select start centroids */
  if(initializer == 0){ /* Random */
    for(i = 0; i < nclusters; i++){
      //srand(m->col+m->row+nclusters+i);
      srand(time(NULL));
      UIVectorAppend(&pre_centroids, rand() % m->row);
    }
  }
  else if(initializer == 1){ /* KMeansppCenters */
    KMeansppCenters(m, nclusters, &pre_centroids, nthreads, s);
  }
  else if(initializer == 2){ /* MDC */
    MDC(m, nclusters, 0, &pre_centroids, nthreads, s);
  }
  else{ /*if(initializer == 3){  MaxDis */
    MaxDis(m, nclusters, 0, &pre_centroids, nthreads, s);
  }

  /* else personal centroid configuration */

  for(i = 0; i < pre_centroids->size; i++){
    for(j = 0; j < m->col; j++){
      centroids->data[i][j] = m->data[pre_centroids->data[i]][j];
    }
  }

  DelUIVector(&pre_centroids);

  UIVectorResize(cluster_labels, m->row);

  it = 0;
  while(shouldStop(centroids, oldcentroids, it, 100) == 0)
  {
    #ifdef DEBUG
    clock_t t = clock();
    #endif

    MatrixCopy(centroids, &oldcentroids);

    #ifdef DEBUG
    t = clock() - t;
    printf("Matrix copy: %f\n", ((double)t)/CLOCKS_PER_SEC);
    t = clock();
    #endif

    //getLabels(m, centroids, (*cluster_labels));
    getLabels_(m, centroids, (*cluster_labels), nthreads);

    #ifdef DEBUG
    t = clock() - t;
    #endif

    #ifdef DEBUG
    printf("getLabels_: %f\n", ((double)t)/CLOCKS_PER_SEC);
    t = clock();
    #endif

    getCentroids(m, (*cluster_labels), &centroids);

    #ifdef DEBUG
    t = clock() - t;
    printf("getCentroids: %f\n", ((double)t)/CLOCKS_PER_SEC);
    #endif
    it++;
  }

  if(_centroids_ == NULL){
    DelMatrix(&centroids);
  }
  DelMatrix(&oldcentroids);
}

void KMeansRandomGroupsCV(matrix* m, size_t maxnclusters, int initializer, size_t groups, size_t iterations, dvector** ssdist, size_t nthreads, ssignal *s)
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

          KMeans(subm, j, initializer, &clusters, &centroids, nthreads, s);

          initMatrix(&distances);
          EuclideanDistance(centroids, predm, &distances, nthreads);

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
void KMeansJumpMethod(matrix* m, size_t maxnclusters, int initializer, dvector** jumps, size_t nthreads, ssignal *s)
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

    KMeans(m, k, initializer, &clusters, &centroids, nthreads, s);

    dist = MatrixMahalanobisDistance(m, centroids);

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
void HierarchicalClustering(matrix* _m, size_t nclusters, uivector** _clusters, matrix **_centroids_, strvector **dendogram, enum LinkageType linktype, size_t nthreads, ssignal *s)
{
  size_t i, j, k, l, m, min_i, min_j;
  char buffer[MAXCHARSIZE];
  char fl2str[50];
  /* calc the matrix distance */
  matrix *distmx, *distmx_new;
  strvector *pointname, *pointname_new, *clusters, *tokens;
  dvector *clusterdist, *tmp;

  initMatrix(&distmx);
  initStrVector(&pointname);
  initStrVector(&clusters);
  initDVector(&clusterdist);

  if(linktype > 2){
    SquaredEuclideanDistance(_m, _m, &distmx, nthreads);
  }
  else{
    EuclideanDistance(_m, _m, &distmx, nthreads);
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

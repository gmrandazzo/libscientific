/* Implements clustering algorithms for data grouping.
 * Copyright (C) 2016-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>

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
  size_t k;
  size_t j;
  mdc_th_args *arg = (mdc_th_args*) arg_;

  for(k = arg->from; k < arg->to; k++){
    double dist = 0.f;
    /* EUCLIDEAN DISTANCE */
    if(arg->metric == 0){
      for(j = 0; j < arg->m->col; j++){
        dist += square(arg->m->data[arg->mdc][j] - arg->m->data[k][j]);
      }
      dist = sqrt(dist);
    }
    /* MANHATAN DISTANCE */
    else if(arg->metric == 1){
      for(j = 0; j < arg->m->col; j++){
        dist += fabs(arg->m->data[arg->mdc][j] - arg->m->data[k][j]);
      }
    }
    /* COSINE DISTANCE */
    else{
      double d_a = 0.f;
      double d_b = 0.f;
      for(j = 0; j < arg->m->col; j++){
        dist += arg->m->data[arg->mdc][j] * arg->m->data[k][j];
        d_a += square(arg->m->data[arg->mdc][j]);
        d_b += square(arg->m->data[k][j]);
      }
      dist = dist/(sqrt(d_a)*sqrt(d_b));
    }
    arg->tmprank->data[k][0] = dist;
    arg->tmprank->data[k][1] = (double)k;
  }
  return 0;
}

/*
 * Most Dscriptive Compound Selection Algorithm
 */
void MDC(matrix* m,
         size_t n,
         int metric,
         uivector *selections,
         size_t nthreads)
{
  size_t i;
  size_t j;
  size_t k;
  size_t l;
  size_t mdc;
  size_t nmdc;
  size_t th;
  double d;
  double dist;
  dvector *vectinfo;
  dvector *rankvector;
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
    CalculateDistance(m, m, dm, nthreads, EUCLIDEAN);
  }
  else if(metric == 1){
    CalculateDistance(m, m, dm, nthreads, MANHATTAN);
  }
  else{
    CalculateDistance(m, m, dm, nthreads, COSINE);
  }

  for(i = 0; i < m->row; i++){
    for(k = 0; k < m->row; k++){
      tmprank->data[k][0] = dm->data[k][i];
      tmprank->data[k][1] = (double)k;
    }

    MatrixSort(tmprank, 0);

    //Calculate the reciprocal of the rank
    d = 2;
    for(k = 0; k < tmprank->row; k++){
      if(k == i){
        vectinfo->data[k] += 1;
      }
      else{
        j = (size_t)tmprank->data[k][1];
        vectinfo->data[j] += 1/d;
        d += 1;
      }
    }
  }
  DelMatrix(&dm);

  nmdc = 0;
  while(1){
    /* Find the compound with largest value in vectinfo */

    dist = vectinfo->data[0];
    mdc = 0;

    for(i = 0; i < vectinfo->size; i++){
      if(vectinfo->data[i] > dist){
        dist = vectinfo->data[i];
        mdc = i;
      }
      else{
        continue;
      }
    }

    nmdc++;

    /*now mdc is the MDC */

    UIVectorAppend(selections, mdc);
    /*Recalculate the disances of the MDC to all the other compounds 
      * and the reciprocal ranks as aboce subtract these reciprocal 
      * ranks from 1 and store in the rank vector R
      */
    /* MULTITHREAD IMPLEMENTATION */
    size_t step = (size_t)ceil((double)m->row/(double)nthreads);
    size_t from = 0;
    size_t to = step;
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

    MatrixSort(tmprank, 0);


    d = 2.;
    for(i = 0; i < tmprank->row; i++){
      j = (size_t)tmprank->data[i][1];
      if(j == mdc){
        rankvector->data[j] = 0.f;/* set to 0 the mdc in order to unselect this..*/
      }
      else{
        rankvector->data[j] = 1 - (1 / d);
        d += 1;
      }
    }

    /* Multiply values in I by the corresponding values in  R. Store the result in I.*/
    for(i = 0 ; i < vectinfo->size; i++){
      vectinfo->data[i] *= rankvector->data[i];
    }

    /*check that all number in the I Exceded 1 if they do then go to "Find the compound with largest value in vectinfo"
    * else stop....
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
      k = 0; /* objects > 1 */
      l = 0; /* extracted objects */
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

  xfree(threads);
  xfree(args);
  DelMatrix(&tmprank);
  DelDVector(&rankvector);
  DelDVector(&vectinfo);
}

/* Start of MDC_Fast helper types and functions */

typedef struct {
  double dist;
  size_t index;
} mdc_pair;

static int compare_mdc_pair(const void *a, const void *b) {
  const mdc_pair *p1 = (const mdc_pair *)a;
  const mdc_pair *p2 = (const mdc_pair *)b;
  if (p1->dist < p2->dist) return -1;
  if (p1->dist > p2->dist) return 1;
  return 0;
}

typedef struct {
  matrix *m;
  dvector *local_vectinfo;
  size_t from, to;
  int metric;
} mdc_fast_init_args;

static void *MDCFastInitialWorker(void *arg_) {
  mdc_fast_init_args *arg = (mdc_fast_init_args*)arg_;
  size_t i, j, k;
  double d;
  mdc_pair *pairs = xmalloc(sizeof(mdc_pair) * arg->m->row);

  for (i = arg->from; i < arg->to; i++) {
    /* Calculate distances from i to all j */
    for (j = 0; j < arg->m->row; j++) {
      pairs[j].index = j;
      if (i == j) {
        pairs[j].dist = 0.0;
        continue;
      }
      
      double dist = 0.0;
      /* Inline distance calculation for speed */
      if (arg->metric == 0) { /* Euclidean */
        for (k = 0; k < arg->m->col; k++) {
           double diff = arg->m->data[i][k] - arg->m->data[j][k];
           dist += diff * diff;
        }
        pairs[j].dist = sqrt(dist);
      } else if (arg->metric == 1) { /* Manhattan */
        for (k = 0; k < arg->m->col; k++) {
           dist += fabs(arg->m->data[i][k] - arg->m->data[j][k]);
        }
        pairs[j].dist = dist;
      } else { /* Cosine */
        double d_a = 0.0, d_b = 0.0;
        for (k = 0; k < arg->m->col; k++) {
          dist += arg->m->data[i][k] * arg->m->data[j][k];
          d_a += square(arg->m->data[i][k]);
          d_b += square(arg->m->data[j][k]);
        }
        if (d_a > 0 && d_b > 0)
            pairs[j].dist = dist / (sqrt(d_a) * sqrt(d_b));
        else
            pairs[j].dist = 0.0; /* Handle zero vectors */
      }
    }
    
    /* Sort distances */
    qsort(pairs, arg->m->row, sizeof(mdc_pair), compare_mdc_pair);
    
    /* Update local_vectinfo (accumulate votes) */
    d = 2.0;
    for (j = 0; j < arg->m->row; j++) {
       size_t idx = pairs[j].index;
       if (idx == i) {
         arg->local_vectinfo->data[idx] += 1.0; 
       } else {
         arg->local_vectinfo->data[idx] += 1.0 / d;
         d += 1.0;
       }
    }
  }
  
  xfree(pairs);
  return NULL;
}

typedef struct {
  matrix *m;
  mdc_pair *distances; 
  size_t mdc_idx;
  size_t from, to;
  int metric;
} mdc_fast_dist_args;

static void *MDCFastDistWorker(void *arg_) {
  mdc_fast_dist_args *arg = (mdc_fast_dist_args*)arg_;
  size_t j, k;

  for (j = arg->from; j < arg->to; j++) {
      arg->distances[j].index = j;
      if (arg->mdc_idx == j) {
          arg->distances[j].dist = 0.0;
          continue;
      }

      double dist = 0.0;
      size_t i = arg->mdc_idx;
      /* Inline distance calculation */
      if (arg->metric == 0) { /* Euclidean */
        for (k = 0; k < arg->m->col; k++) {
           double diff = arg->m->data[i][k] - arg->m->data[j][k];
           dist += diff * diff;
        }
        arg->distances[j].dist = sqrt(dist);
      } else if (arg->metric == 1) { /* Manhattan */
        for (k = 0; k < arg->m->col; k++) {
           dist += fabs(arg->m->data[i][k] - arg->m->data[j][k]);
        }
        arg->distances[j].dist = dist;
      } else { /* Cosine */
        double d_a = 0.0, d_b = 0.0;
        for (k = 0; k < arg->m->col; k++) {
          dist += arg->m->data[i][k] * arg->m->data[j][k];
          d_a += square(arg->m->data[i][k]);
          d_b += square(arg->m->data[j][k]);
        }
        if (d_a > 0 && d_b > 0)
            arg->distances[j].dist = dist / (sqrt(d_a) * sqrt(d_b));
        else
            arg->distances[j].dist = 0.0;
      }
  }
  return NULL;
}

void MDC_Fast(matrix* m,
              size_t n,
              int metric,
              uivector *selections,
              size_t nthreads)
{
  size_t i, j, k, l;
  size_t mdc, nmdc;
  size_t th;
  double d, dist;
  dvector *vectinfo;
  dvector *rankvector;
  mdc_pair *rank_pairs;

  pthread_t *threads;
  mdc_fast_init_args *init_args;
  mdc_fast_dist_args *dist_args;

  NewDVector(&vectinfo, m->row);
  DVectorSet(vectinfo, 0.0);
  NewDVector(&rankvector, m->row);

  threads = xmalloc(sizeof(pthread_t) * nthreads);
  init_args = xmalloc(sizeof(mdc_fast_init_args) * nthreads);

  /* PHASE 1: Initial Scoring (Parallel & Memory Efficient) */
  size_t step = (size_t)ceil((double)m->row / (double)nthreads);
  size_t from = 0;
  
  for (th = 0; th < nthreads; th++) {
    init_args[th].m = m;
    init_args[th].from = from;
    init_args[th].to = (from + step > m->row) ? m->row : from + step;
    init_args[th].metric = metric;
    NewDVector(&init_args[th].local_vectinfo, m->row);
    DVectorSet(init_args[th].local_vectinfo, 0.0);
    
    pthread_create(&threads[th], NULL, MDCFastInitialWorker, (void*)&init_args[th]);
    
    from = init_args[th].to;
  }

  for (th = 0; th < nthreads; th++) {
    pthread_join(threads[th], NULL);
    /* Accumulate local results into main vectinfo */
    for (i = 0; i < m->row; i++) {
        vectinfo->data[i] += init_args[th].local_vectinfo->data[i];
    }
    DelDVector(&init_args[th].local_vectinfo);
  }
  xfree(init_args);

  /* PHASE 2: Selection Loop */
  dist_args = xmalloc(sizeof(mdc_fast_dist_args) * nthreads);
  rank_pairs = xmalloc(sizeof(mdc_pair) * m->row);
  nmdc = 0;

  while(1) {
    /* Find compound with largest value */
    dist = vectinfo->data[0];
    mdc = 0;
    for (i = 1; i < vectinfo->size; i++) {
      if (vectinfo->data[i] > dist) {
        dist = vectinfo->data[i];
        mdc = i;
      }
    }

    nmdc++;
    UIVectorAppend(selections, mdc);

    /* Break if we have enough */
    if (n > 0 && nmdc >= n) break;
    
    /* Stop if all scores are <= 1 and we want "auto" selection (n=0 logic from original code) */
    if (n == 0) {
        k = 0; l = 0;
        for(i = 0; i < vectinfo->size; i++){
            double val = vectinfo->data[i];
            if(val > 1.0) k++;
            else if(FLOAT_EQ(val, 0.0, EPSILON)) l++;
        }
        /* If there are items with 0 < score <= 1, stop */
        if(k != vectinfo->size - l) break;
    }

    /* Recalculate distances from new MDC to all others (Parallel) */
    step = (size_t)ceil((double)m->row / (double)nthreads);
    from = 0;
    for (th = 0; th < nthreads; th++) {
        dist_args[th].m = m;
        dist_args[th].distances = rank_pairs;
        dist_args[th].mdc_idx = mdc;
        dist_args[th].from = from;
        dist_args[th].to = (from + step > m->row) ? m->row : from + step;
        dist_args[th].metric = metric;
        
        pthread_create(&threads[th], NULL, MDCFastDistWorker, (void*)&dist_args[th]);
        from = dist_args[th].to;
    }
    for (th = 0; th < nthreads; th++) {
        pthread_join(threads[th], NULL);
    }

    /* Sort ranks (Serial - N log N is fast enough for loop iteration) */
    qsort(rank_pairs, m->row, sizeof(mdc_pair), compare_mdc_pair);

    /* Calculate penalty factors */
    d = 2.0;
    for (i = 0; i < m->row; i++) {
        j = rank_pairs[i].index;
        if (j == mdc) {
            rankvector->data[j] = 0.0; /* Select probability 0 for already selected */
        } else {
            rankvector->data[j] = 1.0 - (1.0 / d);
            d += 1.0;
        }
    }

    /* Update scores */
    for (i = 0; i < vectinfo->size; i++) {
        vectinfo->data[i] *= rankvector->data[i];
    }
  }

  xfree(threads);
  xfree(dist_args);
  xfree(rank_pairs);
  DelDVector(&rankvector);
  DelDVector(&vectinfo);
}

/*
 * MaxDis object selection.
 */
void MaxDis(matrix* m,
            size_t n,
            int metric,
            uivector *selections,
            size_t nthreads)
{
  size_t i;
  size_t j;
  size_t nobj;
  size_t ntotobj;
  size_t far_away;
  double far;
  double dis;
  matrix *m1;
  matrix *m2;
  matrix *distances;
  uivector *idselection;
  dvector *tmp;
  dvector *mindists;

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

  far_away = 0;
  far = 0.f;
  for(j = 0; j < m->col; j++){
    far += square(c->data[j] - m->data[far_away][j]);
  }
  far = sqrt(far);

  for(i = 1; i < m->row; i++){
    double dst = 0.f;
    for(j = 0; j < m->col; j++){
      dst += square(c->data[j] - m->data[i][j]);
    }
    dst = sqrt(dst);
    if(dst > far){
      far = dst;
      far_away = i;
    }
  }
  DelDVector(&c);

  UIVectorAppend(idselection, far_away);

  /* OLD IMPLEMENTATION SELECTED THE FIRST COMPOUND AS MDC
   * MDC(m, 1, metric, &idselection, nthreads);
   */

  /* 1. Initialise Subset by transferring to it a componund */
  UIVectorAppend(selections, getUIVectorValue(idselection, 0));

  tmp = getMatrixRow(m, getUIVectorValue(idselection, 0));
  MatrixAppendRow(m2, tmp);
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
    initMatrix(&m1);
    initUIVector(&idselection);

    for(i = 0; i < m->row; i++){
      if(UIVectorHasValue(selections, i) == 1){
        tmp = getMatrixRow(m, i);
        MatrixAppendRow(m1, tmp);
        UIVectorAppend(idselection, i);
        DelDVector(&tmp);
      }
      else{
        continue;
      }
    }

    initMatrix(&distances);

    if(metric == 0){
      CalculateDistance(m1, m2, distances, nthreads, EUCLIDEAN);
    }
    else if(metric == 1){
      CalculateDistance(m1, m2, distances, nthreads,  MANHATTAN);
    }
    else{
      CalculateDistance(m1, m2, distances, nthreads,  COSINE);
    }

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

    size_t l = 0;
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
    MatrixAppendRow(m2, tmp);
    DelDVector(&tmp);


    DelDVector(&mindists);
    DelMatrix(&distances);
    DelUIVector(&idselection);
    DelMatrix(&m1);
  }
  DelMatrix(&m2);
}



/* MaxDis_Parallel Helper Structures and Functions */

typedef struct {
  matrix *m;
  dvector *min_dists;
  unsigned char *is_selected;
  size_t last_selected;
  size_t from, to;
  int metric;
} maxdis_worker_args;

static void *MaxDisUpdateWorker(void *arg_) {
  maxdis_worker_args *arg = (maxdis_worker_args*)arg_;
  size_t i, j, k;
  double d, val_i, val_last;
  
  /* Pointers for faster access */
  double **data = arg->m->data;
  double *dists = arg->min_dists->data;
  unsigned char *sel = arg->is_selected;
  size_t last = arg->last_selected;
  size_t cols = arg->m->col;

  for (i = arg->from; i < arg->to; i++) {
    if (sel[i]) continue;

    d = 0.0;
    if (arg->metric == 0) { /* Euclidean */
      for (k = 0; k < cols; k++) {
        double diff = data[i][k] - data[last][k];
        d += diff * diff;
      }
      d = sqrt(d);
    } else if (arg->metric == 1) { /* Manhattan */
      for (k = 0; k < cols; k++) {
        d += fabs(data[i][k] - data[last][k]);
      }
    } else { /* Cosine Distance (1 - Similarity) */
      double dot = 0.0, norm_i = 0.0, norm_last = 0.0;
      for (k = 0; k < cols; k++) {
        val_i = data[i][k];
        val_last = data[last][k];
        dot += val_i * val_last;
        norm_i += val_i * val_i;
        norm_last += val_last * val_last;
      }
      if (norm_i > 0 && norm_last > 0)
        d = 1.0 - (dot / (sqrt(norm_i) * sqrt(norm_last)));
      else
        d = 1.0; /* Treat zero vector as max distance */
    }

    if (d < dists[i]) {
      dists[i] = d;
    }
  }
  return NULL;
}

void MaxDis_Fast(matrix* m,
                     size_t n,
                     int metric,
                     uivector *selections,
                     size_t nthreads)
{
  size_t i, j, k;
  size_t far_away, best_idx;
  size_t nobj;
  double max_min_dist;
  dvector *min_dists;
  dvector *centroid;
  unsigned char *is_selected;
  pthread_t *threads;
  maxdis_worker_args *args;

  if (n == 0 || m->row == 0) return;

  /* Allocations */
  NewDVector(&min_dists, m->row);
  /* Initialize min_dists to infinity */
  for(i=0; i<m->row; i++) min_dists->data[i] = DBL_MAX;

  is_selected = xmalloc(sizeof(unsigned char) * m->row);
  memset(is_selected, 0, sizeof(unsigned char) * m->row);
  
  threads = xmalloc(sizeof(pthread_t) * nthreads);
  args = xmalloc(sizeof(maxdis_worker_args) * nthreads);

  /* Step 1: Select first point (farthest from centroid) */
  NewDVector(&centroid, m->col);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      centroid->data[j] += m->data[i][j];
    }
  }
  for(j = 0; j < m->col; j++) centroid->data[j] /= (double)m->row;

  far_away = 0;
  double max_d = -1.0;
  
  /* This loop is O(N*D), can be parallelized but fast enough usually. keeping serial for simplicity of Phase 1 */
  for(i = 0; i < m->row; i++){
    double d = 0.0;
    for(j = 0; j < m->col; j++){
      double diff = m->data[i][j] - centroid->data[j];
      d += diff * diff;
    }
    /* We use squared euclidean for centroid check to avoid sqrt */
    if (d > max_d) {
        max_d = d;
        far_away = i;
    }
  }
  DelDVector(&centroid);

  /* Add first selection */
  UIVectorAppend(selections, far_away);
  is_selected[far_away] = 1;
  min_dists->data[far_away] = 0.0;

  /* Step 2: Iteratively select remaining points */
  /* We need to select n points total. We have 1. Loop n-1 times. */
  size_t last_selected = far_away;

  for (nobj = 1; nobj < n && nobj < m->row; nobj++) {
      
      /* Phase 2a: Update min_dists in parallel */
      size_t step = (size_t)ceil((double)m->row / (double)nthreads);
      size_t from = 0;

      for (k = 0; k < nthreads; k++) {
          args[k].m = m;
          args[k].min_dists = min_dists;
          args[k].is_selected = is_selected;
          args[k].last_selected = last_selected;
          args[k].from = from;
          args[k].to = (from + step > m->row) ? m->row : from + step;
          args[k].metric = metric;
          
          pthread_create(&threads[k], NULL, MaxDisUpdateWorker, (void*)&args[k]);
          from = args[k].to;
      }

      for (k = 0; k < nthreads; k++) {
          pthread_join(threads[k], NULL);
      }

      /* Phase 2b: Find point with maximum min_dist (Serial) */
      max_min_dist = -1.0;
      best_idx = (size_t)-1;

      for (i = 0; i < m->row; i++) {
          if (is_selected[i]) continue;
          if (min_dists->data[i] > max_min_dist) {
              max_min_dist = min_dists->data[i];
              best_idx = i;
          }
      }

      if (best_idx == (size_t)-1) break; /* Should not happen if nobj < m->row */

      UIVectorAppend(selections, best_idx);
      is_selected[best_idx] = 1;
      min_dists->data[best_idx] = 0.0;
      last_selected = best_idx;
  }

  xfree(threads);
  xfree(args);
  xfree(is_selected);
  DelDVector(&min_dists);
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
  ResizeMatrix(gmap, m->col, 3);

  /* Centering and Scaling to unit variance */
  MatrixColAverage(m, colaverage);
  MatrixColSDEV(m, colscaling);

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

  if(bins_id != NULL){
    /* for each object check what is the bin membership and store in the id bins_id */
    HyperGridMapObjects(m, (*hgm), bins_id);
  }
}


void HyperGridMapObjects(matrix *m, HyperGridModel *hgm, hgmbins **bins_id)
{
  size_t i;
  size_t j;
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
        if x is the minimum then is on 0
        (*bins_id)->hash[i][j] = 0;
      }
      else if(FLOAT_EQ(m->data[i][j], gmap->data[j][1], EPSILON)){
        if x is the maximum then is on max of the grid position in this axis
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
  size_t i;
  size_t j;
  for(i = 0; i < bins_id->nobj; i++){
    for(j = 0; j < bins_id->hash_size-1; j++){
      printf("%zu", bins_id->hash[i][j]);
    }
    printf("%zu\n", bins_id->hash[i][bins_id->hash_size-1]);
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
  size_t i;
  size_t j;
  size_t k;
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
      DVectorAppend(D_min, sqrt(dist));
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

void KMeansppCenters(matrix *m,
                     size_t n,
                     uivector *selections,
                     int nthreads)
{
  size_t i;
  size_t j;
  double y;
  double A;
  double B;
  dvector *D;
  dvector *D_square;
  size_t q = n; /*get the number of clusters*/

  pthread_t *threads = xmalloc(sizeof(pthread_t)*nthreads);
  kmpp_th_args *arg = xmalloc(sizeof(kmpp_th_args)*nthreads);

  NewDVector(&D, m->row); /* vettore distanza di tutti i punti. Per ogni punto c'è un valore di distanza */
  NewDVector(&D_square, m->row); /* vettore distanza di tutti i punti al quadrato.*/
  /* Step 1 Set Random Function */
  UIVectorAppend(selections, (size_t)randInt(0, (int)m->row)); /* random point from 0 to the max data points */

  /* Step 2 */
  while(q > 1){
    size_t nobj = ceil((double)m->row/(double)nthreads);
    size_t from = 0;
    for(i = 0; i < nthreads; i++){
      arg[i].m = m;
      arg[i].selections = selections;
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
    for(i = 0; i < D->size; i++){
      D_square->data[i] = square(D->data[i]);
    }
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
      y = randDouble(B, A);
      if(A >= y && y > B ){
        if(UIVectorHasValue(selections, i) == 1){
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
  DelDVector(&D);
  DelDVector(&D_square);
  xfree(threads);
  xfree(arg);
}

/*This function rank and get the nmaxobj near or far from centroids */
void PruneResults(matrix *m,
                  matrix *centroids,
                  size_t nmaxobj,
                  int type,
                  uivector* clusters,
                  size_t nthreads)
{
  size_t i;
  size_t j;
  size_t n;
  size_t k;
  size_t l;
  double var;
  matrix *submx;
  matrix *subcentroid;
  matrix *distmx;
  dvector *tmp;
  uivector *clusters_;
  uivector *ids;

  NewUIVector(&clusters_, clusters->size);

  UIVectorSet(clusters_, 0);

  for(n = 0; n < centroids->row; n++){
    initMatrix(&subcentroid);

    tmp = getMatrixRow(centroids, n);
    MatrixAppendRow(subcentroid, tmp);
    DelDVector(&tmp);

    initMatrix(&submx);
    initUIVector(&ids);
    for(i = 0; i < clusters->size; i++){
      if(clusters->data[i] == n+1){
        tmp = getMatrixRow(m, i);
        MatrixAppendRow(submx, tmp);
        DelDVector(&tmp);
        UIVectorAppend(ids, i);
      }
      else{
        continue;
      }
    }

    initMatrix(&distmx);
    CalculateDistance(subcentroid, submx, distmx, nthreads, EUCLIDEAN); 

    /* type == 0 = Near Object else Faar object */
    float (*compare_func)(float, float) = (type == 0) ? fminf : fmaxf;
    for (j = 0; j < nmaxobj; j++) {
      l = 0;
      k = ids->data[0];
      var = distmx->data[0][0];
      for (int i = 1; i < distmx->row; i++) {
        if (compare_func(distmx->data[i][0], var) == distmx->data[i][0]){
          var = distmx->data[i][0];
          k = ids->data[i];
          l = i;
        }
      }

      if(FLOAT_EQ(var, MISSING, EPSILON)){
        distmx->data[l][0] = MISSING;
      }
      else{
        distmx->data[l][0] = MISSING;
        clusters_->data[k] = n+1;
      }
    }
  
    DelMatrix(&distmx);
    DelMatrix(&subcentroid);
    DelUIVector(&ids);
    DelMatrix(&submx);
  }

  for(i = 0; i < clusters->size; i++){
    clusters->data[i] = clusters_->data[i];
  }

  DelUIVector(&clusters_);
}

/*
 * 1. Random Centroid made by kmeans++ step 1 algorithm
 * 2. Distance object-centroids
 * 3. Group based on minimum distance
 * 4. If no object moved to group stop; else recompute new centroids and go to step 2
*/

int shouldStop(matrix *centroids,
               matrix *oldcentroids,
               size_t iterations,
               size_t max_iterations)
{
  size_t i;
  size_t j;
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
  size_t i;
  size_t j;
  size_t k;
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
  size_t i;
  size_t from;
  pthread_t *threads = xmalloc(sizeof(pthread_t)*nthreads);
  labels_th_arg *arg = xmalloc(sizeof(labels_th_arg)*nthreads);
  int nobj = ceil((double)m->row/(double)nthreads);
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
  size_t i;
  size_t j;
  size_t k;
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
  size_t i;
  size_t j;
  size_t c_indx;
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
      c_indx = (size_t)randInt(0, (int)m->row);
      for(j = 0; j < m->col; j++){
        new_centroids->data[i][j] = m->data[c_indx][j];
      }
    }
  }

  DelUIVector(&cluster_points);
  MatrixCopy(new_centroids, centroids);
  DelMatrix(&new_centroids);
}

void KMeans(matrix* m,
            size_t nclusters,
            int initializer,
            uivector *cluster_labels,
            matrix *_centroids_,
            size_t nthreads)
{
  size_t i;
  size_t j;
  size_t it;
  matrix *centroids;
  matrix *oldcentroids;
  uivector *pre_centroids;

  if(_centroids_ == NULL){
    NewMatrix(&centroids, nclusters, m->col);
  }
  else{
    centroids = _centroids_;
    ResizeMatrix(centroids, nclusters, m->col);
  }

  NewMatrix(&oldcentroids, centroids->row, centroids->col);

  initUIVector(&pre_centroids);
  /* Step 1. Select start centroids */
  if(initializer == 0){ /* Random */
    for(i = 0; i < nclusters; i++){
      //srand_(m->col+m->row+nclusters+i);
      UIVectorAppend(pre_centroids, randInt(0, m->row));
    }
  }
  else if(initializer == 1){ /* KMeansppCenters */
    KMeansppCenters(m, nclusters, pre_centroids, nthreads);
  }
  else if(initializer == 2){ /* MDC */
    MDC(m, nclusters, 0, pre_centroids, nthreads);
  }
  else{ /*if(initializer == 3){  MaxDis */
    MaxDis(m, nclusters, 0, pre_centroids, nthreads);
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
    
    //getLabels(m, centroids, cluster_labels);
    getLabels_(m, centroids, cluster_labels, nthreads);
    #ifdef DEBUG
    t = clock() - t;
    #endif

    #ifdef DEBUG
    printf("getLabels_: %f\n", ((double)t)/CLOCKS_PER_SEC);
    t = clock();
    #endif
    getCentroids(m, cluster_labels, &centroids);
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

void KMeansRandomGroupsCV(matrix* m,
                          size_t maxnclusters,
                          int initializer,
                          size_t groups,
                          size_t iterations,
                          dvector *ssdist,
                          size_t nthreads)
{

  size_t i;
  size_t j;
  size_t k;
  size_t g;
  size_t n;
  size_t iterations_;
  int a;
  double mindist;
  matrix *gid;
  matrix *subm;
  matrix *predm;
  matrix *centroids;
  matrix *distances;
  uivector *clusters;

  NewMatrix(&gid, groups, (size_t)ceil((double)m->row/(double)groups));

  srand_(groups*m->row*iterations);

  DVectorResize(ssdist, maxnclusters);

  iterations_ = 0;
  while(iterations_ <  iterations){
    /* Divide in group  all the Dataset */
    MatrixSet(gid, -1);

    /* step 1 generate the random groups */
    k = 0;
    for(i = 0; i <  gid->row; i++){
      for(j = 0; j <  gid->col; j++){
        do{
          n = (size_t)randInt(0, m->row);
        } while(ValInMatrix(gid, n) == 1 && k < (m->row));
        if(k < m->row){
          gid->data[i][j] = n;
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
            if((int)gid->data[i][j] != -1)
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
        if((int)gid->data[g][j] != -1)
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
            a = (int)gid->data[i][j]; /* get the row index */
            if(a != -1){ 
              for(n = 0; n < m->col; n++){
                subm->data[k][n] = m->data[a][n];
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
        a = (int)gid->data[g][j];
        if(a != -1){
          for(n = 0; n < m->col; n++){
            predm->data[k][n] = m->data[a][n];
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

        KMeans(subm, j, initializer, clusters, centroids, nthreads);

        initMatrix(&distances);
        CalculateDistance(centroids, predm, distances, nthreads, EUCLIDEAN);

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
          mindist = distances->data[i][0];
          for(k = 1; k < distances->col; k++){ /*for each centroid */
            if(distances->data[i][k] < mindist){
              mindist = distances->data[i][k];
            }
            else{
              continue;
            }
          }
          ssdist->data[j-1] += mindist;
        }
        DelMatrix(&distances);
        DelUIVector(&clusters);
        DelMatrix(&centroids);
      }


      DelMatrix(&subm);
      DelMatrix(&predm);
    }
    iterations_++;
  }

  /* divide all the value of the ssdist for the number of iteration */
  for(i = 0; i < ssdist->size; i++){
    ssdist->data[i] /= (double)iterations;
  }

  DVectNorm(ssdist, ssdist);
  DelMatrix(&gid);
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
void KMeansJumpMethod(matrix* m,
                      size_t maxnclusters,
                      int initializer,
                      dvector *jumps,
                      size_t nthreads)
{
  size_t k;
  size_t i;
  double dist;
  double y;
  double min;
  double max;
  matrix *centroids;
  matrix *xcenter;
  uivector *clusters;
  dvector *d;

  DVectorResize(jumps, maxnclusters);

  NewMatrix(&xcenter, 1, m->col);

  NewDVector(&d, maxnclusters+1);

  y = (double)m->col/2.f;

  for(k = 2; k <= maxnclusters; k++){
    initUIVector(&clusters);
    initMatrix(&centroids);

    KMeans(m, k, initializer, clusters, centroids, nthreads);
    dist = MatrixMatrixDistance(m, centroids);

    printf("dist %f  %f\n", dist, pow(dist, -y));

    if(FLOAT_EQ(dist, 0.f, EPSILON)){
      d->data[k] = 0.f;
    }
    else{
      d->data[k] = pow(dist, -y);
    }

    DelUIVector(&clusters);
    DelMatrix(&centroids);
  }

  for(i = 1; i < d->size; i++){
    jumps->data[i-1] = d->data[i] - d->data[i-1];
  }

  DVectorMinMax(jumps, &min, &max);

  y = max - min;
  for(i = 0; i < jumps->size; i++){
    jumps->data[i] = (jumps->data[i] - min)/ y;
  }

  DelDVector(&d);
  DelMatrix(&xcenter);
}

/* Helper struct for the optimized tree */
typedef struct {
    int left;
    int right;
    double dist;
    size_t size;
} HClusterNode;

void HierarchicalClustering(matrix* _m,
                            size_t nclusters,
                            uivector *_clusters,
                            strvector *dendogram,
                            enum LinkageType linktype,
                            size_t nthreads)
{
  size_t i, j, k;
  size_t n = _m->row;
  size_t n_active;
  size_t min_i, min_j;
  double min_dist;
  
  matrix *distmx;
  HClusterNode *tree;     /* Stores nodes N to 2N-2 */
  int *node_map;          /* Maps row index to tree node index */
  int *cluster_size;      /* Size of cluster at row i */
  unsigned char *active;  /* 1 if row is active, 0 otherwise */
  
  /* Caching for speed optimization */
  double *row_min_val;    /* Minimum distance for row i */
  int *row_min_idx;       /* Index of minimum distance for row i */

  if (n == 0) return;
  if (nclusters > n) nclusters = n;
  if (nclusters == 0) nclusters = 1;

  /* 1. Initialization */
  NewMatrix(&distmx, n, n);
  
  /* Use Squared Euclidean for Ward, else Euclidean */
  if(linktype == ward_linkage){
    CalculateDistance(_m, _m, distmx, nthreads, SQUARE_EUCLIDEAN);
  } else {
    CalculateDistance(_m, _m, distmx, nthreads, EUCLIDEAN);
  }

  /* Fill diagonal with infinity to ignore self-distances */
  for(i=0; i<n; i++) distmx->data[i][i] = DBL_MAX;

  tree = xmalloc(sizeof(HClusterNode) * (n - 1));
  node_map = xmalloc(sizeof(int) * n);
  cluster_size = xmalloc(sizeof(int) * n);
  active = xmalloc(sizeof(unsigned char) * n);
  row_min_val = xmalloc(sizeof(double) * n);
  row_min_idx = xmalloc(sizeof(int) * n);

  for(i=0; i<n; i++) {
      node_map[i] = i; /* Initial leaves are 0..n-1 */
      cluster_size[i] = 1;
      active[i] = 1;
      
      /* Initialize row mins */
      double min_v = DBL_MAX;
      int min_k = -1;
      for(j=0; j<n; j++) {
          if (i==j) continue;
          if (distmx->data[i][j] < min_v) {
              min_v = distmx->data[i][j];
              min_k = j;
          }
      }
      row_min_val[i] = min_v;
      row_min_idx[i] = min_k;
  }

  /* 2. Main Loop */
  n_active = n;
  for (k = 0; k < n - 1; k++) {
      /* Find global minimum among active rows */
      min_dist = DBL_MAX;
      min_i = 0; 
      min_j = 0;

      for(i = 0; i < n; i++) {
          if(!active[i]) continue;
          if(row_min_val[i] < min_dist) {
              min_dist = row_min_val[i];
              min_i = i;
              min_j = row_min_idx[i];
          }
      }

      /* Ensure min_i < min_j for consistent indexing */
      /* Note: row_min_idx might point to inactive if we didn't update perfectly, 
         but our update logic below ensures validity. */
      
      /* Record Merge */
      tree[k].left = node_map[min_i];
      tree[k].right = node_map[min_j];
      tree[k].dist = min_dist;
      tree[k].size = cluster_size[min_i] + cluster_size[min_j];
      
      /* Update active sets */
      active[min_j] = 0; /* min_j is merged into min_i */
      node_map[min_i] = n + k; /* New node index */
      
      /* Update distances for min_i (the merged row) */
      /* We effectively move the new cluster into row min_i */
      
      for (i = 0; i < n; i++) {
          if (!active[i] || i == min_i) continue;
          
          double d_iu = distmx->data[i][min_i];
          double d_iv = distmx->data[i][min_j];
          double d_new = 0.0;

          if (linktype == single_linkage) {
              d_new = (d_iu < d_iv) ? d_iu : d_iv;
          } else if (linktype == complete_linkage) {
              d_new = (d_iu > d_iv) ? d_iu : d_iv;
          } else if (linktype == average_linkage) {
              /* Standard UPGMA: weighted average */
              /* d_new = (size_i*d_iu + size_j*d_iv) / (size_i + size_j) 
                 Using indices min_i and min_j for u and v */
              /* Wait, i is external. u is min_i, v is min_j */
              double n_u = (double)cluster_size[min_i];
              double n_v = (double)cluster_size[min_j];
              d_new = (n_u * d_iu + n_v * d_iv) / (n_u + n_v);
          } else if (linktype == ward_linkage) {
              /* Standard Ward Update */
              double n_i_ = (double)cluster_size[i];
              double n_u = (double)cluster_size[min_i];
              double n_v = (double)cluster_size[min_j];
              double n_sum = n_i_ + n_u + n_v;
              d_new = ((n_i_ + n_u) * d_iu + (n_i_ + n_v) * d_iv - n_i_ * min_dist) / n_sum;
          } else {
              /* Fallback Average */
              d_new = (d_iu + d_iv) / 2.0; 
          }
          
          /* Update Symmetric Matrix */
          distmx->data[i][min_i] = d_new;
          distmx->data[min_i][i] = d_new;
      }
      
      cluster_size[min_i] += cluster_size[min_j];

      /* Update row_min cache */
      /* 1. Update row min_i: requires full scan of its new values */
      double best_val = DBL_MAX;
      int best_idx = -1;
      for(j=0; j<n; j++) {
          if(!active[j] || min_i == j) continue;
          if(distmx->data[min_i][j] < best_val) {
              best_val = distmx->data[min_i][j];
              best_idx = j;
          }
      }
      row_min_val[min_i] = best_val;
      row_min_idx[min_i] = best_idx;

      /* 2. Update other rows: Check if they pointed to min_j (invalid) or min_i (value changed) */
      for(i=0; i<n; i++) {
          if(!active[i] || i == min_i) continue;
          
          int old_target = row_min_idx[i];
          
          /* If the new distance to min_i is smaller than old min, update */
          if (distmx->data[i][min_i] < row_min_val[i]) {
              row_min_val[i] = distmx->data[i][min_i];
              row_min_idx[i] = min_i;
          } 
          /* If it pointed to min_j (now gone) or min_i (value increased?), we must rescan */
          else if (old_target == min_j || old_target == min_i) {
              double m_v = DBL_MAX;
              int m_k = -1;
              for(j=0; j<n; j++) {
                  if(!active[j] || i==j) continue;
                  if(distmx->data[i][j] < m_v) {
                      m_v = distmx->data[i][j];
                      m_k = j;
                  }
              }
              row_min_val[i] = m_v;
              row_min_idx[i] = m_k;
          }
      }
  }

  /* 3. Extract Flat Clusters */
  UIVectorResize(_clusters, n);
  
  /* We have a tree structure. 
     Root is at tree[n-2] (last merge).
     We need to perform a cut.
     Strategy:
     - Assign Cluster ID 1 to Root.
     - While number of clusters < nclusters:
       - Find cluster with highest merge distance? 
       - Actually, the merges are ordered by distance (agglomerative).
       - The last (nclusters-1) merges form the top of the tree.
       - If we undo the last (nclusters-1) merges, we get nclusters.
  */
  
  /* Array to map Tree Node ID (0..2n-2) to Final Cluster ID (1..nclusters) */
  /* Initialize leaves with 0 */
  int *membership = xmalloc(sizeof(int) * (2*n));
  for(i=0; i<2*n; i++) membership[i] = 0;
  
  /* The last merge created node (2n-2).
     Iterate backwards from the last merge.
     We start with 1 cluster (the root).
     Each step backwards splits one cluster into two.
     We stop when we have enough clusters.
  */
  
  /* Label the top nodes */
  /* The nodes kept are those that were NOT children of the last (nclusters-1) merges? 
     No.
     We process the tree from top (last merge) down.
     Assign label 1 to root (node 2n-2).
     Queue: Nodes to process.
     While queue size < nclusters:
       Take node with highest distance (last merges) from queue?
       Split it: Assign label X to left child, label Y to right child.
       
     Simpler:
     Mark the "cut" nodes.
     The nodes formed by merges [n-2 ... n-nclusters] are the top nodes.
     The children of these top nodes that are NOT in the top set are the roots of our final clusters.
  */
  
  int *is_cluster_root = xmalloc(sizeof(int) * (2*n));
  memset(is_cluster_root, 0, sizeof(int) * (2*n));
  
  /* Valid active nodes at the cut level */
  /* Root is valid */
  is_cluster_root[2*n - 2] = 1;
  
  /* We need to perform (nclusters - 1) splits to get nclusters */
  /* We iterate k from n-2 down to n-nclusters.
     For each k, node (n+k) is a merge. If it is marked as a cluster root, 
     we unmark it and mark its children.
  */
  
  /* Merge indices in tree array are 0..n-2. 
     tree[k] creates node (n+k).
     Last merge is index n-2. 
  */
  
  size_t splits_needed = nclusters - 1;
  for (k = 0; k < splits_needed; k++) {
      int merge_idx = (n - 2) - k; /* Go backwards from last merge */
      int node_id = n + merge_idx;
      
      if (is_cluster_root[node_id]) {
          is_cluster_root[node_id] = 0;
          is_cluster_root[tree[merge_idx].left] = 1;
          is_cluster_root[tree[merge_idx].right] = 1;
      }
  }
  
  /* Now assign IDs to leaves by traversing down from marked roots */
  /* Assign sequential IDs 1..nclusters to the marked roots */
  int current_label = 1;
  for (i = 0; i < 2*n-1; i++) {
      if (is_cluster_root[i]) {
          /* Perform BFS/DFS to mark all leaves under this node with current_label */
          /* Simple stack for DFS */
          int *stack = xmalloc(sizeof(int) * n); 
          int top = 0;
          stack[top++] = i;
          
          while(top > 0) {
              int curr = stack[--top];
              if (curr < n) {
                  /* Leaf */
                  _clusters->data[curr] = current_label;
              } else {
                  /* Internal node - index in tree is curr - n */
                  int t_idx = curr - n;
                  stack[top++] = tree[t_idx].left;
                  stack[top++] = tree[t_idx].right;
              }
          }
          xfree(stack);
          current_label++;
      }
  }
  
  /* 4. Generate Dendrogram Strings (if requested) */
  /* Format: "Name1;Name2;Dist" for each merge step k=0..n-2 */
  if (dendogram != NULL) {
      /* We need to reconstruct the names used in original code logic?
         Original code: merged "Name1;Name2;"
         We can just output indices or reconstruct.
         Since dendogram is just a string list, let's output "LeftID;RightID;Dist"
         Note: The original code accumulated names recursively: "1;2;5;"
         We can do that if we really want, but it's expensive. 
         Let's output simple "ID1;ID2;Dist" where IDs refer to leaves or cluster indices?
         Original format: "clustername\tdistance". 
         Where clustername was "leaf1;leaf2;leaf3...".
         
         Let's replicate that behavior efficiently.
         We maintain an array of strings for each active node.
         Only create them if dendogram != NULL.
      */
      char **node_names = xmalloc(sizeof(char*) * (2*n));
      for(i=0; i<n; i++) {
          size_t len = snprintf(NULL, 0, "%zu", i);
          node_names[i] = xmalloc(len + 2);
          snprintf(node_names[i], len + 2, "%zu;", i); /* Original appended ';' */
      }
      
      for(k=0; k<n-1; k++) {
          int left = tree[k].left;
          int right = tree[k].right;
          int new_id = n + k;
          
          /* Combine names */
          /* Name = "LeftNameRightName" (Original logic did simple concat) */
          size_t len_l = strlen(node_names[left]);
          size_t len_r = strlen(node_names[right]);
          node_names[new_id] = xmalloc(len_l + len_r + 1);
          strcpy(node_names[new_id], node_names[left]);
          strcat(node_names[new_id], node_names[right]);
          
          /* Add to dendogram vector */
          /* Format: Name + Distance */
          size_t len_entry = snprintf(NULL, 0, "%s%f", node_names[new_id], tree[k].dist);
          char *entry = xmalloc(len_entry + 1);
          snprintf(entry, len_entry + 1, "%s%f", node_names[new_id], tree[k].dist);
          StrVectorAppend(dendogram, entry);
          xfree(entry);
      }
      
      for(i=0; i<2*n-1; i++) if(node_names[i]) xfree(node_names[i]);
      xfree(node_names);
  }

  /* Cleanup */
  xfree(is_cluster_root);
  xfree(membership);
  xfree(row_min_idx);
  xfree(row_min_val);
  xfree(active);
  xfree(cluster_size);
  xfree(node_map);
  xfree(tree);
  DelMatrix(&distmx);
}


/* clustering.h
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

#ifndef CLUSTERING_H
#define CLUSTERING_H

#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

typedef struct{
  uivector *left;
  uivector *right;
  double dist;
} NODE;

/* Most Descriptive Compound Method
 * Input:
 *  - m: matrix of data in coordinate
 *  - n: number of object to select.
 *
 * Output:
 *  - selecttions: vector of id selected..
 */
void MDC(matrix* m,
         size_t n,
         int metric,
         uivector *selections,
         size_t nthreads,
         ssignal *s);


/*
 * MaxMin Dissimilarity Method
 *
 * Input:
 * - m: matrix of data in coordinate
 * - n: number of object to select.
 *
 * Output:
 * - selections: vector of id selected
 */
void MaxDis(matrix* m,
            size_t n,
            int metric,
            uivector *selections,
            size_t nthreads,
            ssignal *s);

void MaxDis_Fast(matrix* m,
                 size_t n,
                 int metric,
                 uivector *selections,
                 size_t nthreads,
                 ssignal *s);


/*
 * HyperGrid Map Selections
 *
 * Input:
 * - m: N-Dimensional matrix of objects
 * - grid_size: grid size
 *
 * Output:
 * - bins_id: for each object we define the hyper grid bin of appartenaince
 * - HyperGridModel: Specific model
 *
 * Bins numbering:
 *   ___ ___ ___ ___
 *  |_8_|_9_|_10|_11|
 *  |_4_|_5_|_6_|_7_|
 *  |_0_|_1_|_2_|_3_|
 */

typedef struct{
  matrix *gmap;  /* grid map (min, max, step)*/
  dvector *colaverage;
  dvector *colscaling;
  size_t gsize;  /* grid size*/
  double bsize;  /* total number of bins */
} HyperGridModel;

typedef struct{
  size_t **hash;
  size_t nobj, hash_size;
} hgmbins;

/* Allocate the model */
void NewHyperGridMap(HyperGridModel **hgm);

/* Delete the model */
void DelHyperGridMap(HyperGridModel **hgm);

/*Create the HyperGridMap */
void HyperGridMap(matrix* m, size_t grid_size, hgmbins** bins_id, HyperGridModel **hgm);

/*Extract an object from hypergridmap*/
void HyperGridMapObjects(matrix *m, HyperGridModel *hgm, hgmbins **bins_id);

void PrintHGMBins(hgmbins *bins_id);

void DelHGMBins(hgmbins **bins_id);

/*
 * KMeans++ Centers Method
 *
 * Input:
 * - m: matrix of data in coordinate
 * - n: number of object to select.
 *
 * Output:
 * - selections: vector of id selected
 */
void KMeansppCenters(matrix *m,
                     size_t n,
                     uivector *selections,
                     int nthreads,
                     ssignal *s);

/*
 * KMeans Clustering
 *
 * Imput:
 *  - m: matrix with N-Dimensional coordinate
 *  - nclusters: Number of clusters desired
 *  - initializer: parameter for start kmeans clustering by: 0 Random centroids (Kmeans)
 *                                                           1 David Arthur initialization (Kmeans++)
 *                                                           2 MDC inizialization
 *                                                           3 MaxMinDis inizialization
 * Output:
 * - clusters: vector of size m->row where for each object is defined the cluster membership
 */
void KMeans(matrix* m,
            size_t nclusters,
            int initializer,
            uivector *cluster_labels,
            matrix *_centroids_,
            size_t nthreads,
            ssignal *s);

void PruneResults(matrix *m, matrix *centroids, size_t nmaxobj, int type, uivector* clusters, size_t nthreads);

void KMeansRandomGroupsCV(matrix* m,
                          size_t maxnclusters,
                          int initializer,
                          size_t groups,
                          size_t iterations,
                          dvector *ssdist,
                          size_t nthreads,
                          ssignal *s);

void KMeansJumpMethod(matrix* m,
                      size_t maxnclusters,
                      int initializer,
                      dvector *jumps,
                      size_t nthreads,
                      ssignal *s);

enum LinkageType {
  single_linkage = 0,
  complete_linkage = 1,
  average_linkage = 2,
  ward_linkage = 3
};

void HierarchicalClustering(matrix* _m,
                            size_t nclusters,
                            uivector *_clusters,
                            matrix *_centroids,
                            strvector *dendogram,
                            enum LinkageType linktype,
                            size_t nthreads,
                            ssignal *s);

#endif

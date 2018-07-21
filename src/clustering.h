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
void MDC(matrix *m, size_t n, int metric, uivector **selections, ssignal *s);


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
void MaxDis(matrix* m, size_t n, int metric, uivector** selections, ssignal *s);


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
  dvector *mult; /* multiplier */
  size_t gsize;  /* grid size*/
  size_t bsize;  /* total number of bins */
} HyperGridModel;

/* Allocate the model */
void NewHyperGridMap(HyperGridModel **hgm);

/* Delete the model */
void DelHyperGridMap(HyperGridModel **hgm);

/*Create the HyperGridMap */
void HyperGridMap(matrix* m, size_t grid_size, uivector** bins_id, HyperGridModel **hgm);

/*Extract an object from hypergridmap*/
void HyperGridMapObjects(matrix *m, HyperGridModel *hgm, uivector **bins_id);

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
void KMeansppCenters(matrix *m, size_t n, uivector **selections, ssignal *s);

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
void KMeans(matrix* m, size_t nclusters, int initializer, uivector** clusters, matrix **centroids, ssignal *s);

void PruneResults(matrix *m, matrix *centroids, size_t nmaxobj, int type, uivector* clusters);

void KMeansRandomGroupsCV(matrix* m, size_t maxnclusters, int initializer, size_t groups, size_t iterations, dvector **ssdist, ssignal *s);

void KMeansJumpMethod(matrix* m, size_t maxnclusters, int initializer, dvector** jumps, ssignal *s);

enum LinkageType {
  single_linkage = 0,
  complete_linkage = 1,
  average_linkage = 2,
  ward_linkage = 3
};

void HierarchicalClustering(matrix* _m, size_t nclusters, uivector** _clusters, matrix **_centroids_, strvector **dendogram, enum LinkageType linktype, ssignal *s);

#endif

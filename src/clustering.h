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
 * - m: matrix of data in coordinate
 * - n: number of object to select.
 * 
 * Output:
 * - selections: vector of id selected
 */
void HyperGridMap(matrix* m, size_t step, size_t n, uivector** selections, ssignal *s);

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
  sl = 0,
  cl = 1,
  al = 2,
  wl = 3
};

void HierarchicalClustering(matrix* _m, size_t nclusters, uivector** _clusters, matrix **_centroids_, strvector **dendogram, enum LinkageType linktype, ssignal *s);

#endif

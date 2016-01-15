#ifndef METRICSPACE_H
#define METRICSPACE_H

#include "matrix.h"
#include "vector.h"

void EuclideanDistance(matrix *m1, matrix *m2, matrix** distances);

void SquaredEuclideanDistance(matrix *m1, matrix *m2, matrix **distances);

void ManhattanDistance(matrix *m1, matrix *m2, matrix** distances);

void CosineDistance(matrix *m1, matrix *m2, matrix** distances);

double MahalanobisDistance(matrix* g1, matrix* g2);

#endif

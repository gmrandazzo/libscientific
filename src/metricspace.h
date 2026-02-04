/* Defines metric spaces and distance calculations.
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

#ifndef METRICSPACE_H
#define METRICSPACE_H

#include "matrix.h"
#include "vector.h"

/* Description: calculate a distance between two matrix with the same column
 * number. Multi-thread implementation.
 */
enum cmethod {EUCLIDEAN = 0,
              SQUARE_EUCLIDEAN,
              MANHATTAN,
              COSINE}; 

void CalculateDistance(matrix *m1, 
                       matrix *m2,
                       matrix *distances,
                       size_t nthreads,
                       enum cmethod distance_method);

/* Description: calculate the euclidean distance between two matrix with the same column number.
 *              Single thread implementation
 */
void EuclideanDistance_ST(matrix *m1, matrix *m2, matrix *distances);

/* Description: calculate the square euclidean distance between two matrix with the same column number.
 *              Single thread implementation
 */
void SquaredEuclideanDistance_ST(matrix *m1, matrix *m2, matrix *distances);

/* Description: calculate the manhattan distance between two matrix with the same column number.
 *              Single thread implementation
 */
void ManhattanDistance_ST(matrix *m1, matrix *m2, matrix *distances);

/* Description: calculate the cosine distance between two matrix with the same column number.
 *              Single thread implementation
 */
void CosineDistance_ST(matrix *m1, matrix *m2, matrix *distances);

/* Descriptuon matrix-matrix distance defined as euclidean distance and frobenius norm of it */
double MatrixMatrixDistance(matrix *m1, matrix *m2);

/* Description: convert the matrix mi (mxn) into a covariance distance map (mxn) */
void CovarianceDistanceMap(matrix* mi, matrix *mo);

/*
 * Distance matrix caalculation in a fast and less memory consuming way
 */
/* Description: calculate the index from square to condensed form giving the row and column position*/
size_t square_to_condensed_index(size_t i, size_t j, size_t n);

/* Description: calculate the euclidean distance matrix in a condensed way */
void EuclideanDistanceCondensed(matrix* m, dvector *distances, size_t nthreads);

/* Description: calculate the square euclidean distance in a condensed way */
void SquaredEuclideanDistanceCondensed(matrix *m, dvector *distances, size_t nthreads);

/* Description: calculate the manhattan distance in a condensed way */
void ManhattanDistanceCondensed(matrix *m, dvector *distances, size_t nthreads);

/* Description: calculate the cosine distance in a condensed way */
void CosineDistanceCondensed(matrix *m, dvector *distances, size_t nthreads);

#endif

/* metricspace.h
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

#ifndef METRICSPACE_H
#define METRICSPACE_H

#include "matrix.h"
#include "vector.h"

/* Description: calculate the euclidean distance between two matrix with the same column number */
void EuclideanDistance(matrix *m1, matrix *m2, matrix** distances);

/* Description: calculate the square euclidean distance between two matrix with the same column number */
void SquaredEuclideanDistance(matrix *m1, matrix *m2, matrix **distances);

/* Description: calculate the manhattan distance between two matrix with the same column number */
void ManhattanDistance(matrix *m1, matrix *m2, matrix** distances);

/* Description: calculate the cosine distance between two matrix with the same column number */
void CosineDistance(matrix *m1, matrix *m2, matrix** distances);

/* Description: calculate the malanobis distance between two matrix with the same column number */
double MahalanobisDistance(matrix* g1, matrix* g2);

/* Description: cubic spline interpolation */
void cubic_spline_interpolation(matrix *xy, size_t npoints, matrix **interp_xy);

/* Description: calculate area of a curve.
 * If intervals > 0 will interpolate the curve
 */
double curve_area(matrix *xy, size_t intervals);

/* Description: calculate the ROC curve giving an y_true and an y_score */
void ROC(dvector *y_true, dvector *y_score, matrix **roc, double *auc);

/* Description: calculate the Precision-Recall curve giving an y_true and an y_score */
void PrecisionRecall(dvector *y_true, dvector *y_score,  matrix **pr, double *ap);

#endif

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

void EuclideanDistance(matrix *m1, matrix *m2, matrix** distances);

void SquaredEuclideanDistance(matrix *m1, matrix *m2, matrix **distances);

void ManhattanDistance(matrix *m1, matrix *m2, matrix** distances);

void CosineDistance(matrix *m1, matrix *m2, matrix** distances);

double MahalanobisDistance(matrix* g1, matrix* g2);

#endif

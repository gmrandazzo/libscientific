/* cpca.h
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

#ifndef CPCA_H
#define CPCA_H
#include <stdio.h>
#include "array.h"
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define CPCACONVERGENCE 1e-3

typedef struct {
  array *b_scores;
  array *b_loadings;
  matrix *sscores;
  matrix *sweights;
  dvector *b_scaling;
  dvector *expvar;
  dvector *colaverage;
  dvector *colscaling;
} CPCAMODEL;


void CPCA(array *x, size_t npc, size_t scaling, CPCAMODEL *m);

#endif
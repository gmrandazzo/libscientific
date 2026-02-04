/* Implements Independent Component Analysis (ICA).
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

#ifndef ICA_H
#define ICA_H

#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define ICACONVERGENCE 1e-8

typedef struct{
  matrix *S;
  matrix *W;
  matrix *whitening_matrix;
  dvector *colaverage;
  dvector *colscaling;
} ICAMODEL;

void NewICAModel(ICAMODEL **m);
void DelICAModel(ICAMODEL **m);

void ICA(matrix *mx,
         size_t scaling,
         size_t n_signals,
         ICAMODEL *model);

void ICA_ext(matrix *mx,
         size_t scaling,
         size_t n_signals,
         double alpha,
         double thresh,
         size_t max_iter,
         ICAMODEL *model);

void ICASignalPredictor(matrix *mx,
                        ICAMODEL *model,
                        matrix *p_signals);

void PrintICA(ICAMODEL *m);

#endif
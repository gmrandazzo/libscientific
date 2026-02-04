/* Provides data preprocessing techniques like scaling and centering.
 * Copyright (C) 2022-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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

#ifndef PREPROCESSING_H
#define PREPROCESSING_H

#include "vector.h"
#include "matrix.h"
#include "tensor.h"
#include "list.h"
#include "scientificinfo.h"


/**
 * Matrix Prepreprocess
 * This method will center and scale your "orig" matrix and store the result
 * on "trans".
 *
 * colaverage and colascaling:  are intended to be used for future
 * transformations.
 * type: this will define the scaling type that you want to apply:
 *         -1: no scaling no centering
 *          0: only centering
 *          1: centering and column standard deviation scaling
 *          2: centering and column root mean square scaling
 *          3: centering and column pareto scaling aka sqrt(stdev)
 *          4: centering and column min<->max range scaling
 *          5: centering and level scaling
 *
 * For autoscaling see:
 * Centering scaling and trasfomrations: improving the biological
 * information content of metabolomics dataset
 * A van den Berg
 * BMC Genomics 2006, 7:142  doi:101 186/147-214-7-142
 */
void MatrixPreprocess(matrix *orig,
                      int type,
                      dvector *colaverage,
                      dvector *colscaling,
                      matrix *trans);

/**
 * Matrix Whitening method.
 * This method will transform the original input matrix
 * into a new matrix where the new variables are uncorrelated
 * each others and they will have variance equal to 1.
 */
void MatrixWhitening(matrix *orig,
                    matrix *whitening_matrix,
                    matrix *whiten_matrix);

/**
 * Tensor Prepreprocess
 * This method will center and scale your "orig" tensor and store the result
 * on "trans".
 *
 * colaverage and colascaling:  are intended to be used for future
 * transformations.
 * type: this will define the scaling type that you want to apply:
 *         -1: no scaling no centering
 *          0: only centering
 *          1: centering and column standard deviation scaling
 *          2: centering and column root mean square scaling
 *          3: centering and column pareto scaling aka sqrt(stdev)
 *          4: centering and column min<->max range scaling
 *          5: centering and level scaling
 *
 * For autoscaling see:
 * Centering scaling and trasfomrations: improving the biological
 * information content of metabolomics dataset
 * A van den Berg
 * BMC Genomics 2006, 7:142  doi:101 186/147-214-7-142
 */
void TensorPreprocess(tensor *orig,
                      int type,
                      dvectorlist *colaverages,
                      dvectorlist *colscalings,
                      tensor *trans);

#endif

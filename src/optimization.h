/* Implements optimization algorithms.
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

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "matrix.h"
#include "numeric.h"

double NelderMeadSimplex(double (*func)(), dvector *x0, dvector *step, double xtol, size_t iter, dvector *best);

#endif

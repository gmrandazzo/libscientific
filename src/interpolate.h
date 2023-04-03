/* interpolate.h
*
* Copyright (C) <2018>  Giuseppe Marco Randazzo
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

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "matrix.h"
#include "vector.h"

/* Description: calculate the natural cubic spline interpolation equations 
 * 
 * Algorithm reference:
 * Natural Cubic Spline: Numerical Analysis book Richard L. Burden and J. Douglas Faires
 * pag. 149
 * ISBN-13: 978-0-538-73351-9
 * 
 * N.B.: If two points have different y but share same x the algorithm will fail
 */
void cubic_spline_interpolation(matrix *xy, matrix *S);

/* Description: predict using the  natural cubic spline 
 *              interpolation equations a vector of x
 */
void cubic_spline_predict(dvector *x_, matrix *S, dvector *y_pred);

/* Description: interpolate x and y using the natural cubic 
 *              spline equations and get directly the interpolation.
 */
void interpolate(matrix *xy, size_t npoints, matrix *interp_xy);

#endif

/* modelvalidation.h
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

#ifndef MODELVALIDATION_H
#define MODELVALIDATION_H
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

void YScrambling(matrix *mx, matrix *my,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t block,
                        size_t valtype, size_t rgcv_group, size_t rgcv_iterations,
                        matrix **r2q2scrambling, size_t nthreads, ssignal *s);

/*
 * Description: Calculate the Bootstrap Random group cross validation.
 */
void BootstrapRandomGroupsCV(matrix *mx, matrix *my, size_t xautoscaling, /*Inputs*/
                       size_t yautoscaling, size_t nlv, size_t group, size_t iterations, /*Inputs*/
                      matrix **q2y, matrix **sdep, matrix **bias, matrix **predicted_y, /*Ouputs*/
                      matrix **pred_residuals, size_t nthreads, ssignal *s); /*Ouputs*/

/*
 * Description: Calculate the leave one out cross validation.
 */
void LOOCV(matrix *mx, matrix *my, size_t xautoscaling, size_t yautoscaling, size_t nlv,/*Inputs*/
                        matrix **q2y, matrix **sdep, matrix **bias, /*Ouputs*/
                        matrix **predicted_y, matrix **pred_residuals, /*Ouputs*/
                        size_t ntreads, ssignal *s); /*Inputs*/
#endif

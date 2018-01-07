/* vectorspace.h
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

#ifndef VECTORSPACE_H
#define VECTORSPACE_H

#include <stdio.h>
#include <stdlib.h>

#include "vector.h"

/*
* Definitions:
* a = proportion of 1s that the variables share in the same positions
* b = proportion of 1s in the first variable and 0s in second variable
* in the same positions
* c = proportion of 0s in the first variable and 1s in second variable
* in the same positions
* d = proportion of 0s that both variables share in the same positions.
*/


/* Peirce coefficient */
double SPeir1(double a, double b, double c, double d);
double SPeir2(double a, double b, double c, double d);

/* Doolittle, Pearson coefficient */
double SDoo(double a, double b, double c, double d);

/* Yule (1900), Montgomery and Crittenden */
double SYule1(double a, double b, double c, double d);

#endif

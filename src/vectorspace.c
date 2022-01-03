/* vectorspace.c
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

#include "vectorspace.h"
#include "vector.h"

/*
 * Definitions:
 * a = proportion of 1s that the variables share in the same positions
 * b = proportion of 1s in the first variable and 0s in second variable
 * in the same positions
 * c = proportion of 0s in the first variable and 1s in second variable
 * in the same positions
 * d = proportion of 0s that both variables share in the same positions.
 * p1 = a + b proportion of 1s in the first variable
 * p2 = a + c proportion of 1s in the second variable
 * q1 = c + d proportion of 0s in the first variable
 * q2 = b + d proportion of 0s in the second variable.
 */

/* SPeir1 = (ad − bc)/p1q1 */
double SPeir1(double a, double b, double c, double d)
{
  return (a*d)-(b*c) / ((a+b)*(c+d));
}

/* SPeir2 = (ad − bc)/p2q2 */
double SPeir2(double a, double b, double c, double d)
{
  return (a*d)-(b*c) / ((a+c)*(b+d));
}

/* SDoo = (ad − bc)^2/p1p2q1q2*/
double SDoo(double a, double b, double c, double d)
{
  return square((a*d)-(b*c))/((a+b)*(a+c)*(c+d)*(b+d));
}

/* SYule1 = (ad − bc) / (ad + bc) */
double SYule1(double a, double b, double c, double d)
{
  return square((a*d)-(b*c))/((a+d)+(b+c));
}

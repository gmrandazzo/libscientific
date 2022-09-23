/* ica.c
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

#include <stdio.h>
#include <stdlib.h>

#include "ica.h"

/*
 * See here:
 * https://towardsdatascience.com/independent-component-analysis-ica-in-python-a0ef0db0955e
 */
void ICA(matrix *mx, size_t scaling, size_t npc, ICAMODEL *model, ssignal *s)
{
  // 1. Calculate the PCA extracting
  // 2. Rotation of loadings
  // Central Limit Theorem:
  // S = W*X
  // W is the demelange matrix
  // X = loadings from PCA

}

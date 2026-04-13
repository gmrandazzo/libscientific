/* Implements a generic linked list data structure.
 * Copyright (C) 2017-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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

#ifndef LIST_H
#define LIST_H

#include <stdio.h>
#include <stdlib.h>
#include "vector.h"

/* string vector */
typedef struct{
  dvector **d;
  size_t size;
} dvectorlist;

void initDVectorList(dvectorlist **lst);
void NewDVectorList(dvectorlist **lst, size_t size_);
void DelDVectorList(dvectorlist **lst);
void DVectorListAppend(dvectorlist *lst, dvector *d);


#endif

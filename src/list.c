/* list.c
*
* Copyright (C) <2017>  Giuseppe Marco Randazzo
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

#include "list.h"
#include "numeric.h"
#include "memwrapper.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>


void initDVectorList(dvectorlist **lst){
  (*lst) = xmalloc(sizeof(dvectorlist));
  (*lst)->d = NULL;
  (*lst)->size = 0;
}

void NewDVectorList(dvectorlist **lst, size_t size_)
{
  (*lst) = xmalloc(sizeof(dvectorlist));
  (*lst)->size = size_;
  (*lst)->d = xmalloc(sizeof(dvector*)*size_);
}

void DelDVectorList(dvectorlist **lst)
{
  size_t i;
  for(i = 0; i < (*lst)->size; i++)
    DelDVector(&(*lst)->d[i]);
  xfree((*lst)->d);
  xfree((*lst));
}

void DVectorListAppend(dvectorlist **lst, dvector *d)
{
  size_t i;
  (*lst)->size += 1;
  (*lst)->d = xrealloc((*lst)->d, sizeof(dvector*)*(*lst)->size);
  NewDVector(&(*lst)->d[(*lst)->size-1], d->size);
  for(i = 0; i < d->size; i++){
    (*lst)->d[(*lst)->size-1]->data[i] = d->data[i];
  }
}

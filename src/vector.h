/* vector.h
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

#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <stdlib.h>

#define MAXCHARSIZE 6114

/* string vector */
typedef struct{
  char **data;
  size_t size;
} strvector;

/* Allocate a strvector */
void NewStrVector(strvector **s, size_t size);
void initStrVector(strvector **s);

/* Delete a strvector */
void DelStrVector(strvector **s);

void StrVectorResize(strvector *s, size_t size_);

/* Append a value to a strvector */
void StrVectorAppend(strvector *s, char *str);
void StrVectorAppendInt(strvector *s, int val);
void StrVectorAppendDouble(strvector *s, double val);

void setStr(strvector *s, size_t i, char *str);
char* getStr(strvector *s, size_t i);

/* Append to a strvector an other strvector */
strvector *StrVectorExtend(strvector *d1, strvector *d2);
void PrintStrVector(strvector *s);

/*operations with string*/

/*Trim a string: remove \n*/
char *Trim(char *s);
/*Split a string and fill the splitted string into a strvector*/
void SplitString(char *str, char *sep, strvector *tokens);


/* double vector */
typedef struct{
  double *data;
  size_t size;
} dvector;

/* Allocate a dvector */
void NewDVector(dvector **d, size_t size);
void initDVector(dvector **d);

/* Delete a dvector */
void DelDVector(dvector **d);

/* Resize a dvector */
void DVectorResize(dvector *d, size_t size_);

void PrintDVector(dvector *v);

/* Append a value to a dvector */
void DVectorAppend(dvector *d, double val);

/* Remove a value to a dvector */
void DVectorRemoveAt(dvector *d, size_t indx);

/* Copy a Dvector from dsrc: source to ddst: destination */
void DVectorCopy(dvector *dsrc, dvector *ddst);

/* Append to a dvector an other dvector */
dvector *DVectorExtend(dvector *d1, dvector *d2);

void setDVectorValue(dvector *d, size_t id, double val);
double getDVectorValue(dvector *d, size_t id);

/*check if dvector has a value val. Return 0 if is present, 1 if is not present*/
int DVectorHasValue(dvector *d, double val);

/*Vector operations*/
void DVectorSet(dvector *v, double val);
double DVectorDVectorDotProd(dvector *v1, dvector *v2); /* product between two vector */

double DvectorModule(dvector *v); /* get the Dvector Module */
void DVectNorm(dvector *v, dvector *nv); /* vector normalizing */
void DVectorDVectorDiff(dvector *v1, dvector *v2, dvector *v3);
void DVectorDVectorSum(dvector *v1, dvector *v2, dvector *v3);
void DVectorMinMax(dvector *v, double *min, double *max);
void DVectorMean(dvector *d, double *mean);
void DVectorMedian(dvector *d, double *mean);
void DVectorSDEV(dvector *d, double *sdev);
void DVectorSort(dvector *v);


/* Int Vector */
typedef struct{
  int *data;
  size_t size;
} ivector;

/* Allocate a ivector */
void NewIVector(ivector** d, size_t size_);
void initIVector(ivector **d);

/* Delete a ivector */
void DelIVector(ivector **d);

/* Append a ivector */
void IVectorAppend(ivector *d, int val);

/* Remove a value to a ivector */
void IVectorRemoveAt(ivector *d, size_t indx);

/* Append to a ivector an other ivector */
ivector *IVectorExtend(ivector *d1, ivector *d2);

void setIVectorValue(ivector *d, size_t id, int val);
int getIVectorValue(ivector *d, size_t id);
int IVectorHasValue(ivector *d, int val);

void IVectorSet(ivector *d, int val);


/* size_t VECTOR */
typedef struct{
  size_t *data;
  size_t size;
} uivector;

/* Allocate a uivector */
void NewUIVector(uivector **d, size_t size);
void initUIVector(uivector **d);

/* Delete a uivector */
void DelUIVector(uivector **d);

void UIVectorResize(uivector *d, size_t size_);

void PrintUIVector(uivector *v);

/* Append a uivector */
void UIVectorAppend(uivector *d, size_t val);

/* Remove a value to uivector */
void UIVectorRemoveAt(uivector *d, size_t indx);

/* Append to a uivector an other uivector */
uivector *UIVectorExtend(uivector *d1, uivector *d2);

void setUIVectorValue(uivector *d, size_t id, size_t val);
size_t getUIVectorValue(uivector *d, size_t id);

/*
 * UIVectorHasValue: if the vector "u" has the "id" value
 *                   return 0, else return 1
 */
int UIVectorHasValue(uivector *u, size_t id);
int UIVectorIndexOf(uivector *u, size_t id);

void UIVectorSet(uivector *d, size_t val);
void SortUIVector(uivector *d);

#endif

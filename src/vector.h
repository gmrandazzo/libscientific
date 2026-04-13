/* Defines the vector data structure and operations.
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

#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <stdlib.h>

/**
 * A dynamic array of strings.
 * The strvector struct is a dynamic array of strings represented
 * by a pointer to a character pointer (char data) and its size (size_t size).
 * 
 * - **data** A pointer to a character pointer representing the array of strings.
 * - **size** The size of the array of strings.
 */
typedef struct{
  char **data;
  size_t size;
} strvector;

/**
 * Allocate a strvector with a specific size
 */
void NewStrVector(strvector **s, size_t size);

/**
 * initialize an empty strvector 
 */
void initStrVector(strvector **s);

/**
 *  Delete a strvector
 */
void DelStrVector(strvector **s);

/**
 * Resize a strvector
 */
void StrVectorResize(strvector *s, size_t size_);

/**
 * Append a string value to a strvector
 */
void StrVectorAppend(strvector *s, char *str);

/**
 * Append a integer value to a strvector
 */
void StrVectorAppendInt(strvector *s, int val);

/**
 * Append a double value to a strvector
 */
void StrVectorAppendDouble(strvector *s, double val);

void setStr(strvector *s, size_t i, char *str);
char* getStr(strvector *s, size_t i);

/**
 * Extend a strvector with an other strvector
 */
strvector *StrVectorExtend(strvector *s1, strvector *s2);

/**
 * Print to video a strvector 
 */
void PrintStrVector(strvector *s);

/*operations with string*/

/**
 * Trim a string: remove
 */
char *Trim(char *s);

/**
 * Split a string and fill the splitted string into a strvector
 */
void SplitString(char *str, char *sep, strvector *tokens);


/**
 * A dynamic array of doubles.
 * The dvector struct is a dynamic array of double represented
 * by a pointer to a double and its size (size_t size).
 * 
 * - **data** A pointer to a double type.
 * - **size** The size of the array of strings.
 */
typedef struct{
  double *data;
  size_t size;
} dvector;

/**
 * Allocate a dvector with a specific size
 */
void NewDVector(dvector **d, size_t size);

/**
 * initialize an empty dvector 
 */
void initDVector(dvector **d);

/**
 * Delete a dvector
 */
void DelDVector(dvector **d);

/**
 * Resize a dvector
 */
void DVectorResize(dvector *d, size_t size_);

/**
 * Print a dvector to video
 */
void PrintDVector(dvector *v);

/**
 * Append a value to a dvector
 */
void DVectorAppend(dvector *d, double val);

/**
 * Remove a value to a dvector
 */
void DVectorRemoveAt(dvector *d, size_t indx);

/**
 * Copy a Dvector from dsrc: source to ddst: destination
 */
void DVectorCopy(dvector *dsrc, dvector *ddst);

/*
 * Extend a divector with another dvector
 */
dvector *DVectorExtend(dvector *d1, dvector *d2);

void setDVectorValue(dvector *d, size_t id, double val);
double getDVectorValue(dvector *d, size_t id);

/**
 * Check if dvector has a value val.
 * Return 0 if is present, 1 if is not present
 */
int DVectorHasValue(dvector *d, double val);

/*
 * Vector operations
 */
void DVectorSet(dvector *v, double val);
double DVectorDVectorDotProd(dvector *v1, dvector *v2); /* product between two vector */

double DvectorModule(dvector *v); /* get the Dvector Module */
void DVectNorm(dvector *v, dvector *nv); /* vector normalizing */
void DVectorDVectorDiff(dvector *v1, dvector *v2, dvector *v3);
void DVectorDVectorSum(dvector *v1, dvector *v2, dvector *v3);
void DVectorMinMax(dvector *v, double *min, double *max);
void DVectorMean(dvector *d, double *mean);
void DVectorMedian(dvector *d, double *median);
void DVectorSDEV(dvector *d, double *sdev);
void DVectorSort(dvector *v);


/**
 * A dynamic array of integers.
 * The dvector struct is a dynamic array of integers represented
 * by a pointer to a integers and its size (size_t size).
 * 
 * - **data** A pointer to a double type.
 * - **size** The size of the array of strings.
 */
typedef struct{
  int *data;
  size_t size;
} ivector;

/**
 * Allocate a ivector with a specific size
 */
void NewIVector(ivector** d, size_t size);

/**
 * Initialize an ivector
 */
void initIVector(ivector **d);

/**
 * Delete a ivector
 */
void DelIVector(ivector **d);

/**
 * Print to video an ivector
 */
void PrintIVector(ivector *v);

/**
 * Append a value to an ivector
 */
void IVectorAppend(ivector *d, int val);

/**
 * Remove a value to a ivector
 */
void IVectorRemoveAt(ivector *d, size_t indx);

/**
 * Extend an ivector with another ivector
 */
ivector *IVectorExtend(ivector *d1, ivector *d2);

void setIVectorValue(ivector *d, size_t id, int val);
int getIVectorValue(ivector *d, size_t id);
int IVectorHasValue(ivector *d, int val);

void IVectorSet(ivector *d, int val);


/**
 * A dynamic array of unsigned integers.
 * The dvector struct is a dynamic array of unsigned integers represented
 * by a pointer to a unsigned integers and its size (size_t size).
 * 
 * - **data** A pointer to unsigned integer type.
 * - **size** The size of the array of strings.
 */
typedef struct{
  size_t *data;
  size_t size;
} uivector;

/**
 * Allocate a uivector with a specific size
 */
void NewUIVector(uivector **d, size_t size);

/**
 * Initialize an empty uivector
 */
void initUIVector(uivector **d);

/**
 * Delete a uivector
 */
void DelUIVector(uivector **d);

/**
 * Resize an uivector
 */
void UIVectorResize(uivector *d, size_t size_);

/**
 * Print an uivector 
 */
void PrintUIVector(uivector *v);

/**
 * Append a value to an uivector
 */
void UIVectorAppend(uivector *d, size_t val);

/**
 * Remove a value to an uivector
 */
void UIVectorRemoveAt(uivector *d, size_t indx);

/**
 * Extend an uivector with another uivector
 */
uivector *UIVectorExtend(uivector *d1, uivector *d2);

void setUIVectorValue(uivector *d, size_t id, size_t val);
size_t getUIVectorValue(uivector *d, size_t id);

/**
 * Check if an uivector has a value.
 * if found return 0, else return 1
 */
int UIVectorHasValue(uivector *u, size_t id);
int UIVectorIndexOf(uivector *u, size_t id);

void UIVectorSet(uivector *d, size_t val);
void SortUIVector(uivector *d);

#endif

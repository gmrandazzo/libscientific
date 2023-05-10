/* vector.c
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

#include "vector.h"
#include "numeric.h"
#include "memwrapper.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* strdup is not standard so add this code to copy a string to another place dynamically

 char *strdup(const char *str)
{
    int n = strlen(str) + 1;
    char *dup = malloc(n);
    if(dup)
    {
        strcpy(dup, str);
    }
    return dup;
}

extern char *strdup(const char *s);
*/

void initStrVector(strvector** s){
    (*s) = xmalloc(sizeof(strvector));
    (*s)->data = NULL;
    (*s)->size = 0;
}

void NewStrVector(strvector **s, size_t size)
{
    size_t i;
    (*s) = xmalloc(sizeof(strvector));
    (*s)->size = size;
    (*s)->data = xmalloc(sizeof(char*)*size);
    for(i = 0; i < (*s)->size; i++){
      (*s)->data[i] = xmalloc(sizeof(char));
    }
}

void DelStrVector(strvector **s)
{
  if(s != NULL){
    size_t i;
    for(i = 0; i < (*s)->size; i++)
      xfree((*s)->data[i]);
    xfree((*s)->data);
    xfree((*s));
  }
}

void setStr(strvector *s, size_t i, char *str)
{
  /*memcpy(s->data[i], str, strlen(str)+1);*/
  //   strcpy(s->data[i], str);
  xfree(s->data[i]);
  s->data[i] = strdup(str);
}

char* getStr(strvector *s, size_t i)
{
  return s->data[i];
}

void StrVectorAppend(strvector *s, char *str)
{
  size_t i;
  size_t size = s->size+1;
  strvector *tmp;
  NewStrVector(&tmp, s->size);

  for(i = 0; i < s->size; i++){
    setStr(tmp, i, getStr(s, i));
  }

  s->data = xrealloc(s->data, sizeof(char*)*size);
  s->data[s->size] = xmalloc(sizeof(char));
  s->size += 1;

  for(i = 0; i < tmp->size; i++){
    setStr(s, i, getStr(tmp, i));
  }

  setStr(s, s->size-1, str);

  DelStrVector(&tmp);
}


void StrVectorAppendInt(strvector *s, int val)
{
  size_t size = s->size+1;
  s->data = xrealloc(s->data, sizeof(char*)*size);
  s->size += 1;
  size_t length = snprintf(NULL, 0, "%d", val);
  s->data[s->size-1] = xmalloc(sizeof(char)*length+1);
  snprintf(s->data[s->size-1], length + 1, "%d", val);
}

void StrVectorAppendDouble(strvector *s, double val)
{
  size_t size = s->size+1;
  s->data = xrealloc(s->data, sizeof(char*)*size);
  s->size += 1;
  size_t length = snprintf(NULL, 0, "%f", val);
  s->data[s->size-1] = xmalloc(sizeof(char)*length+1);
  snprintf(s->data[s->size-1], length + 1, "%f", val);
}

strvector *StrVectorExtend(strvector *s1, strvector *s2)
{
  size_t i;
  strvector *sext;
  NewStrVector(&sext, (*s1).size+(*s2).size);
  for(i = 0; i < (*s1).size; i++) {
      (*sext).data[i] = (*s1).data[i];
  }
  for(i = 0; i < (*s2).size; i++) {
      (*sext).data[i+(*s1).size] = (*s2).data[i];
  }
  return sext;
}

void StrVectorResize(strvector *s, size_t size_)
{
  size_t i;
  if(s->size > 0){
    for(i = 0; i < s->size; i++)
      xfree(s->data[i]);
    xfree(s->data);
  }

  s->size = size_;
  s->data = xmalloc(sizeof(char*)*size_);
  for(i = 0; i < s->size; i++){
    s->data[i] = xmalloc(sizeof(char)*256);
    memset(s->data[i], 0, 256);
  }
}

void PrintStrVector(strvector* s)
{
  if(s->size > 0){
    size_t i;
    for(i = 0; i < s->size; i++){
      printf("[%u]: %s\n", (unsigned int)i, getStr(s, i));
    }
  }
}

char *Trim(char *s)
{
  if(!s)
    return NULL;   /* handle NULL string */
  if(!*s)
    return s;      /* handle empty string */
  char *cp1; /* for parsing the whole s */
  char *cp2; /* for shifting & padding */
  /* skip leading spaces, shift remaining chars */
  for (cp1=s; isspace(*cp1); cp1++ ); /* skip leading spaces, via cp1 */
  for (cp2=s; *cp1; cp1++, cp2++) /*shift left remaining chars, via cp2 */
    *cp2 = *cp1;
  *cp2-- = 0; /* mark new end of string for s */

  /* replace trailing spaces with '\0' */
  while ( cp2 > s && isspace(*cp2) )
    *cp2-- = 0; /* pad with '\0's*/
  return s;
}

void SplitString(char *str, char *sep, strvector *tokens)
{
  char *tl=NULL, *saveptr;
  /* Create a buffer of the correct length and copy the string into the buffer */
  char  *buffer = Trim(strdup(str));
  /* Tokenize */
  for (tl = strtok_r(buffer, sep, &saveptr); tl; tl = strtok_r(NULL, sep, &saveptr)){
    StrVectorAppend(tokens, tl);
  }
  xfree(buffer);
}


/* DVECTOR */
void initDVector(dvector **d) {
    (*d) = xmalloc(sizeof(dvector));
    (*d)->data = NULL;
    (*d)->size = 0;
}

void NewDVector(dvector **d, size_t size)
{
    size_t i;
    (*d) = xmalloc(sizeof(dvector));
    (*d)->size = size;
    (*d)->data = xmalloc(sizeof(double)*size);
    for(i = 0; i < (*d)->size; i++)
        (*d)->data[i] = +0.f;
}

void DelDVector(dvector **d )
{
  if(d != NULL){
    xfree((*d)->data);
    xfree((*d));
  }
}

void DVectorResize(dvector *d, size_t size_)
{
  size_t i;
  if(d->size > 0){
    xfree(d->data);
  }

  d->data = xmalloc(sizeof(double)*size_);
  d->size = size_;
  for(i = 0; i < d->size; i++){
    d->data[i] = +0.f;
  }
}

void PrintDVector(dvector *v)
{
  if(v->size > 0){
    size_t i;
    for(i = 0; i < v->size; i++){
      printf("[%u]: %.4f\n", (unsigned int)i, getDVectorValue(v, i));
    }
  }
}

void DVectorAppend(dvector *d, double val)
{
  size_t size = d->size+1;
  d->data = xrealloc(d->data, sizeof(double)*size);
  d->size += 1;
  d->data[size-1] = val;
}

void DVectorRemoveAt(dvector *d, size_t indx)
{
  if(indx < d->size){
    // assuming that sizeOfArray is the count of valid elements in the array
    int elements_to_move = d->size - indx - 1;
    memmove( &d->data[indx], &d->data[indx+1], elements_to_move * sizeof(double));
    d->size -=1;
  }
  else{
    return;
  }
}

void DVectorCopy(dvector* dsrc, dvector *ddst)
{
  size_t i;
  if(ddst->size == 0){
    ddst->size= dsrc->size;
    ddst->data = xmalloc(sizeof(double)*dsrc->size);
    for(i = 0; i < ddst->size; i++){
      ddst->data[i] = +0.f;
    }
  }
  else{
    ddst->size = dsrc->size;
    ddst->data = xrealloc(ddst->data, sizeof(double)*dsrc->size);
    for(i = 0; i < ddst->size; i++){
      setDVectorValue(ddst, i, +0.f);
    }
  }

  for(i = 0; i < ddst->size; i++){
    ddst->data[i] = dsrc->data[i];
  }
}

dvector *DVectorExtend(dvector *d1, dvector *d2)
{
  size_t i;
  dvector *dext;
  NewDVector(&dext, (*d1).size+(*d2).size);
  for(i = 0; i < (*d1).size; i++) {
      (*dext).data[i] = (*d1).data[i];
  }
  for(i = 0; i < (*d2).size; i++) {
      (*dext).data[i+(*d1).size] = (*d2).data[i];
  }
  return dext;
}


void setDVectorValue(dvector *d, size_t id, double val)
{
  if(id < d->size)
    (*d).data[id] = val;
  else{
    fprintf(stdout,"setDVectorValue Error: vector id %d out of range.\n", (int)id);
    fflush(stdout);
    abort();
  }
}

double getDVectorValue(dvector *d, size_t id)
{
  if(id < d->size)
    return d->data[id];
  else{
    fprintf(stdout,"getDVectorValue Error: vector id %d out of range\n.", (int)id);
    fflush(stdout);
    abort();
  }
}

int DVectorHasValue(dvector* d, double val)
{
  size_t i;
  for(i = 0; i < d->size; i++){
    if(FLOAT_EQ(d->data[i], val, EPSILON)){
      return 0;
    }
    else{
      continue;
    }
  }
  return 1;
}

void DVectorSet(dvector *v, double val)
{
  size_t i;
  for(i = 0; i < (*v).size; i++){
    (*v).data[i] = val;
  }
}

double DVectorDVectorDotProd(dvector *v1, dvector *v2)
{
  /* v*t_v = Sum(v[i]*t_v[i])
     the vors must have the same size */
  size_t i;
  double p  = +0.f;
  for(i = 0; i < v1->size; i++){
    if(FLOAT_EQ(v1->data[i], MISSING, 1e-1) ||
       FLOAT_EQ(v2->data[i], MISSING, 1e-1)){
      continue;
    }
    else{
      p += v1->data[i]*v2->data[i];
    }
  }
  return p;
}

double DvectorModule(dvector *v)
{
  size_t i;
  double sum = +0.f;
  for(i = 0; i < v->size; i++){
    if(FLOAT_EQ(v->data[i], MISSING, 1e-1)){
      continue;
    }
    else{
      sum += v->data[i]*v->data[i];
    }
  }
  return sqrt(sum);/* module of the for v */
}

void DVectNorm(dvector *v, dvector *nv)
{
  size_t i;
  double mod;

  mod = DvectorModule(v);

  if((*nv).size != 0 && (*nv).size <= v->size){ /* store the normalized vor value to an other vor named n_v */
    for(i = 0; i < v->size; i++){
      if(FLOAT_EQ(v->data[i], MISSING, 1e-1)){
        nv->data[i] = MISSING;
      }
      else{
        nv->data[i] = v->data[i]/mod;
      }
    }
  }
  else{
    fprintf(stdout,"DVectNorm Error! Problem while normalizing the vector!\n");
    fflush(stdout);
    abort();
  }
}

void DVectorDVectorDiff(dvector *v1, dvector *v2, dvector *v3)
{
  size_t i;
  if(v1->size == v2->size){
    for(i = 0; i < v1->size; i++)
      DVectorAppend(v3, v1->data[i]-v2->data[i]);
  }
  else{
    fprintf(stderr, "Error! Unable to compute DVectorDVectorDiff. Different Vector size.\n");
    abort();
  }
}

void DVectorDVectorSum(dvector *v1, dvector *v2, dvector *v3)
{
  size_t i;
  if(v1->size == v2->size){
    for(i = 0; i < v1->size; i++)
      DVectorAppend(v3, getDVectorValue(v1, i)+getDVectorValue(v2, i));
  }
  else{
    fprintf(stderr, "Error! Unable to compute DVectorDVectorDiff. Different Vector size.\n");
    fflush(stderr);
    abort();
  }
}

void DVectorMinMax(dvector* v, double* min, double* max)
{
  if(v->size > 0){
    size_t i;
    if(min != NULL && max != NULL){
      (*min) = (*max) = getDVectorValue(v, 0);
    }
    else{
      if(min != NULL){
        (*min) = getDVectorValue(v, 0);
      }

      if(max != NULL){
        (*max) = getDVectorValue(v, 0);
      }
    }

    for(i = 0; i < v->size; i++){
      if(min != NULL){
        if(getDVectorValue(v, i) < (*min)){
          (*min) = getDVectorValue(v, i);
        }
      }

      if(max != NULL){
        if(getDVectorValue(v, i) > (*max)){
          (*max) = getDVectorValue(v, i);
        }
      }
    }
  }
  else{
    fprintf(stderr, "Error! Empty DVector\n");
    fflush(stderr);
    abort();
  }
}

int cmp(const void *a, const void *b)
{
    double fa = *(double *) a;
    double fb = *(double *) b;
    return (fa > fb) - (fa < fb);
}

void DVectorMedian(dvector* d, double *median)
{
  (*median) = MISSING;
  qsort(d->data, d->size, sizeof(double), cmp);

  if (d->size%2 == 0){
    (*median) = (d->data[d->size/2] + d->data[(d->size/2) - 1])/2.f;
  }
  else
    (*median) = d->data[d->size/2];
}

void DVectorMean(dvector* d, double* mean)
{
  size_t i;
  (*mean) = +0.f;
  for(i = 0; i < d->size; i++)
    (*mean) += d->data[i];
  (*mean) /= d->size;
}

void DVectorSDEV(dvector* d, double* sdev)
{
  size_t i;
  double mean;
  DVectorMean(d, &mean);

  (*sdev) = +0.f;
  for(i = 0; i < d->size; i++)
    (*sdev) += square(d->data[i] - mean);
  (*sdev) /= d->size;
  (*sdev) = sqrt((*sdev));
}

void DVectorSort(dvector* v)
{
  qsort(v->data, v->size, sizeof(v->data[0]), cmp);
}

/* INT VECTOR */
void initIVector(ivector **d){
  (*d) = xmalloc(sizeof(ivector));
  (*d)->data = NULL;
  (*d)->size = 0;
}

void NewIVector(ivector **d, size_t size)
{
  long int i;
  (*d) = xmalloc(sizeof(ivector));
  (*d)->size = size;
  (*d)->data = xmalloc(sizeof(int)*size);
  for(i = 0; i < (*d)->size; i++)
    (*d)->data[i] = 0;
}

void DelIVector(ivector **d )
{
  if(d != NULL){
    xfree((*d)->data);
    xfree((*d));
  }
}

void PrintIVector(ivector* v)
{
  if(v->size > 0){
    size_t i;
    for(i = 0; i < v->size; i++){
      printf("[%u]: %i\n", (unsigned int)i, (int)getIVectorValue(v, i));
    }
  }
}

void IVectorAppend(ivector *d, int val)
{
  long int size = d->size+1;
  d->data = xrealloc(d->data, sizeof(int)*size);
  d->size += 1;
  d->data[d->size-1] = val;
}

void IVectorRemoveAt(ivector *d, size_t indx)
{
  if(indx < d->size){
    // assuming that sizeOfArray is the count of valid elements in the array
    int elements_to_move = d->size - indx - 1;
    memmove(&d->data[indx], &d->data[indx+1], elements_to_move * sizeof(int));
    d->size -=1;
  }
  else{
    return;
  }
}

ivector *IVectorExtend(ivector *d1, ivector *d2)
{
  int i;
  ivector *dext;
  NewIVector(&dext, (*d1).size+(*d2).size);
  for(i = 0; i < (*d1).size; i++){
    (*dext).data[i] = (*d1).data[i];
  }
  for(i = 0; i < (*d2).size; i++){
    (*dext).data[i+(*d1).size] = (*d2).data[i];
  }
  return dext;
}

void setIVectorValue(ivector* d, size_t id, int val)
{
  if(id < d->size)
    (*d).data[id] = val;
  else{
    fprintf(stdout,"setIVectorValue Error: vector id %d out of range.\n", (int)id);
    fflush(stdout);
  }
}

int getIVectorValue(ivector* d, size_t id)
{
  if(id < d->size)
    return d->data[id];
  else{
    fprintf(stdout,"getIVectorValue Error: vector id %d out of range\n.", (int)id);
    fflush(stdout);
    abort();
    return MISSING;
  }
}

int IVectorHasValue(ivector* d, int val)
{
  size_t i;
  for(i = 0; i < d->size; i++){
    if(d->data[i] == val){
      return 0;
    }
    else{
      continue;
    }
  }
  return 1;
}


void IVectorSet(ivector* d, int val)
{
  size_t i;
  for(i = 0; i < (*d).size; i++){
    (*d).data[i] = val;
  }
}

/* size_t VECTOR */

void initUIVector(uivector **d){
  (*d) = xmalloc(sizeof(uivector));
  (*d)->data = NULL;
  (*d)->size = 0;
}

void NewUIVector(uivector **d, size_t size)
{
  size_t i;
  (*d) = xmalloc(sizeof(uivector));
  (*d)->size = size;
  (*d)->data = xmalloc(sizeof(size_t)*size);
  for(i = 0; i < (*d)->size; i++)
    (*d)->data[i] = 0;
}

void DelUIVector(uivector **d )
{
  if(d != NULL){
    xfree((*d)->data);
    xfree((*d));
  }
}

void UIVectorResize(uivector *d, size_t size_)
{
  size_t i;
  if(d->size > 0){
    xfree(d->data);
  }

  d->data = xmalloc(sizeof(size_t*)*size_);
  d->size = size_;
  for(i = 0; i < d->size; i++)
    setUIVectorValue(d, i, 0);
}

void PrintUIVector(uivector* v)
{
  if(v->size > 0){
    size_t i;
    for(i = 0; i < v->size; i++){
      printf("[%u]: %u\n", (unsigned int)i, (unsigned int)getUIVectorValue(v, i));
    }
  }
}

void UIVectorAppend(uivector *d, size_t val)
{
  size_t size = d->size+1;
  d->data = xrealloc(d->data, sizeof(size_t)*size);
  d->size += 1;
  d->data[d->size-1] = val;
}

void UIVectorRemoveAt(uivector *d, size_t indx)
{
  if(indx < d->size){
    int elements_to_move = d->size - indx - 1;
    memmove( &d->data[indx], &d->data[indx+1], elements_to_move * sizeof(size_t));
    d->size -=1;
  }
  else{
    return;
  }
}

uivector *UIVectorExtend(uivector *d1, uivector *d2)
{
  size_t i;
  uivector *dext;
  NewUIVector(&dext, (*d1).size+(*d2).size);
  for(i = 0; i < (*d1).size; i++){
    (*dext).data[i] = (*d1).data[i];
  }
  for(i = 0; i < (*d2).size; i++){
    (*dext).data[i+(*d1).size] = (*d2).data[i];
  }
  return dext;
}

void setUIVectorValue(uivector* d, size_t id, size_t val)
{
  if(id < d->size)
    (*d).data[id] = val;
  else{
    fprintf(stdout,"setIVectorValue Error: vector id %d out of range.\n", (int)id);
    fflush(stdout);
  }
}

size_t getUIVectorValue(uivector* d, size_t id)
{
  if(id < d->size)
    return d->data[id];
  else{
    fprintf(stdout,"getUIVectorValue Error: vector id %u out of range. Max range: %u\n.", (unsigned int)id, (unsigned int)d->size);
    fflush(stdout);
    abort();
    return MISSING;
  }
}

int UIVectorHasValue(uivector* u, size_t id)
{
  size_t i;
  for(i = 0; i < u->size; i++){
    if(u->data[i] == id){
      return 0;
    }
    else{
      continue;
    }
  }
  return 1;
}

int UIVectorIndexOf(uivector *u, size_t id)
{
  size_t i;
  for(i = 0; i < u->size; i++){
    if(u->data[i] == id){
      return i;
    }
    else{
      continue;
    }
  }
  return -1;
}

void UIVectorSet(uivector* d, size_t val)
{
  size_t i;
  for(i = 0; i < (*d).size; i++){
    (*d).data[i] = val;
  }
}

int intcmp(const void *v1, const void *v2)
{
  return (*(int *)v1 - *(int *)v2);
}

void SortUIVector(uivector* d)
{
  qsort(d->data, d->size, sizeof(d->data[0]), intcmp);
}

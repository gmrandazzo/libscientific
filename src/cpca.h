#ifndef CPCA_H
#define CPCA_H
#include <stdio.h>
#include "array.h"
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define CPCACONVERGENCE 1e-3

typedef struct {
  array *b_scores;
  array *b_loadings;
  matrix *sscores;
  matrix *sweights;
  dvector *b_scaling;
  dvector *expvar;
  dvector *colaverage;
  dvector *colscaling;
} CPCAMODEL;


void CPCA(array *x, size_t npc, size_t scaling, CPCAMODEL *m);

#endif
#ifndef UPCA_H
#define UPCA_H
#include "array.h"
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

typedef struct{
  matrix *scores;
  array *loadings;
  dvector *varexp;
  matrix *colaverage;
  matrix *colscaling;
} UPCAMODEL;

void NewUPCAModel(UPCAMODEL **m);
void DelUPCAModel(UPCAMODEL **m);

int CheckArray(array *X_);
void UPCA(array *X_, size_t npc, size_t autoscaling, UPCAMODEL *m, ssignal *s);

void UPCAScorePredictor(array *X_, UPCAMODEL *model, size_t npc, matrix **pscores);
void UPCAIndVarPredictor(matrix *T, array *P, matrix *colaverage, matrix *colsdev,  size_t npc, array **X);

void PrintUPCAModel(UPCAMODEL *m);

#endif

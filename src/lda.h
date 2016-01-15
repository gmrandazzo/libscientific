#ifndef LDA_H
#define LDA_H
#include <stdio.h>
#include "matrix.h"
#include "vector.h"
#include "array.h"
#include "scientificinfo.h"


typedef struct{
  matrix *inv_cov;
  array *features;
  array *mnpdf;
  matrix *evect;
  matrix *mu;
  matrix *fmean;
  matrix *fsdev;
  dvector *eval;
  dvector *pprob;
  uivector *classid;
  size_t nclass;
  size_t class_start;
  dvector *sens, *spec, *ppv, *npv, *acc;
} LDAMODEL;

void NewLDAModel(LDAMODEL **m);
void DelLDAModel(LDAMODEL **m);
void PrintLDAModel(LDAMODEL *m);

void LDA(matrix *mx, uivector *y, LDAMODEL *lda);
/*prediction
 * OUTPUT:
 *  - predicted features
 *  - probability
 *  - multivariate normal probability distribution of features
 *  - class prediction 
 */
void LDAPrediction(matrix *mx, LDAMODEL *lda, matrix **pfeatures, matrix **probability, matrix **mnpdf, uivector **prediction);

void LDARandomGroupsCV(matrix *mx, uivector *my, size_t group, size_t iterations, dvector **sens, dvector **spec, dvector **ppv, dvector **npv, dvector **acc, size_t nthreads, ssignal *s);

void LDALOOCV(matrix* mx, uivector* my, dvector** sens, dvector** spec, dvector** ppv, dvector** npv, dvector **acc, size_t nthreads, ssignal *s);

#endif

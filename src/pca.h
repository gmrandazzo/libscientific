#ifndef PCA_H
#define PCA_H
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define PCACONVERGENCE 1e-8

typedef struct{
  matrix *scores;
  matrix *loadings;
  matrix *dmodx;
  dvector *varexp;
  dvector *colaverage;
  dvector *colscaling;
} PCAMODEL;

void NewPCAModel(PCAMODEL **m);
void DelPCAModel(PCAMODEL **m);

void calcVarExpressed(double ss, dvector *eval, dvector **varexp);
double calcObjectDistance(matrix *m);

void PCA(matrix *mx, size_t scaling, size_t npc, PCAMODEL *model, ssignal *s);

void PCAScorePredictor(matrix *mx, PCAMODEL *model, size_t npc, matrix **pscores);

void PCAIndVarPredictor(matrix *t, matrix *p, dvector *colaverage, dvector *colsdev,  size_t npc, matrix **mx);

void PCARSquared(matrix *mx, PCAMODEL *model, size_t npc, dvector** r2);

void PCARankValidation(matrix *mx, size_t npc, size_t scaling, size_t group, size_t iterations, dvector **r, ssignal *s);

/*Compute the residual matrix for a specific number of component.
 * mx = matrix of origin model
 * model = PCA model where are stored scores and loadings
 * pc = max component to extract the residual matrix
 * rmx = residual matrix of output. must be initialized with initMatrix(&rmx)
 */
void GetResidualMatrix(matrix *mx, PCAMODEL *model, size_t pc, matrix **rmx);

void PrintPCA(PCAMODEL *m);

#endif

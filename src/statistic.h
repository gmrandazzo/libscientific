#ifndef STATISTIC_H
#define STATISTIC_H

#include "matrix.h"
#include "vector.h"

/*
 * TP = True PositivePredictedValue
 * FN = False Negative
 * 
 * Sensitivity = TP / TP + FN
 */
void Sensitivity(dvector *tp, double thmin, double thmax, double thstep, matrix **s);

/*
 * TP = True Positive
 * FP = False Positive
 * 
 * PPV = TP / TP + FP
 */
void PositivePredictedValue(dvector *dtp, dvector *dtn, double thmin, double thmax, double thstep, matrix **p);


/*
 * This function code the matrix to a sign matrix 
 */
void MatrixCode(matrix *inmx, matrix *outmx);

void BifactorialMatrixExpansion(matrix* inmx, matrix* outmx);

/*
 * This function study the variable effect through the yates algorithm
 */
void YatesVarEffect(matrix *mx, dvector *veff);

#endif

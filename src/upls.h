#ifndef UPLS_H
#define UPLS_H
#include "array.h"
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define UPLSCONVERGENCE 1e-3

typedef struct{
  matrix *xscores;
  array *xloadings;
  array *xweights;
  dvector *xvarexp;
  matrix *yscores;
  array *yloadings;
  dvector *b;
  matrix *xcolaverage;
  matrix *xcolscaling;
  matrix *ycolaverage;
  matrix *ycolscaling;
  dvector *r2x_model;
  dvector *r2x_validation;
  array *r2y_model; /* each order correspond to a principal component; each row correspond to an y dependent variable and each col to an order */
  array *r2y_validation;  /* each order correspond to a principal component; each row correspond to an y dependent variable and each col to an order */
  array *q2y;
  array *sdep; /* Standard Deviation over Prediction */
  array *sdec; /* Standard Deviation over Recalculating */
  array *recalculated_y;
  array *predicted_y;
  array *recalc_residuals;
  array *pred_residuals;
  array *q2y_yscrambling;
  array *sdep_yscrambling;
} UPLSMODEL;

/*
 * Description: Create a new UPLSMODEL
 */
void NewUPLSModel(UPLSMODEL **m);

/*
 * Description: Delete an UPLSMODEL
 */
void DelUPLSModel(UPLSMODEL **m);

int CheckArrays(array *X_, array *Y_);

/*
 *  Description:
 * 
 *  UPLS : Unfold PLS is applied to multilinear data compute a bilinear modeling (bilinear weighs and loadings). 
 *  This means that models describe an equal or higher amount of covariance by using an equal or higher number of model parameters
 * 
 * ----------------------
 *  Variable Explanation
 * ----------------------
 * 
 * xscores->row represent the number of objects that is equal to X_->m[i]->row; 
 * xscores->col represent the number of pc that is less equal to X_->order and depends on the number of principal component computed
 * 
 * xloadings->order represent the number of pc computed.
 * xloadings->m[i]->row represent the number of descriptor that is equal to X_->m[i]->col
 * xloadings->m[i]->col represent the number of order present in the array X_ and is equal to X_->order
 * 
 * xvar is the variance explained for each component for the X matrix
 * 
 * xcolaverage->row and xcolsdev->row are the number of variable for each submatrix in X_, so these are equal to X_->m[i]->col
 * 
 * xcolaverage->col and xcolsdev->col are the number of order present in X_, so these are equal to X_->order
 * 
 * 
 * yscores->row represent the number of objects that is equal to X_->m[i]->row and Y_->m[i]->row; 
 * yscores->col represent the number of pc that is less equal to X_->order and depends on the number of principal component computed
 * 
 * yloadings->order represent the number of pc computed.
 * yloadings->m[i]->row represent the number of descriptor for each matrix in the array Y_->m[i]->col
 * yloadings->m[i]->col represent the number of order present in the array Y_->order
 * 
 * yvar is the variance explained for each component for the Y matrix
 * 
 * 
 * ycolaverage->row and ycolsdev->row are the number of variable for each submatrix in Y_, so these are equal to Y_->m[i]->col
 * 
 * ycolaverage->col and ycolsdev->col are the number of order present in Y_, so these are equal to Y_->order
 * 
 * recalc_residuals & recalculated_y  have size row = number of objects  and col = number of column y * number of principal component and order is the same order of the Y dependent variable. Ex. if we calculate 5 component in 2 y the size will be 10
 */
void UPLS(array *X_, array *Y_, size_t npc, size_t xautoscaling, size_t yautoscaling, UPLSMODEL *m, ssignal *s);

/* Description: Predict Score variable through the UPLSMODEL m.
 * This Function is used by UPLSRSquared and UPLSRSquared_SSErr_SStot
 */
void UPLSScorePredictor(array *X_, UPLSMODEL *m, size_t npc, matrix **ptscores);

/*
 * Description: Predict the dependent variable through UPLSMODEL m.
 * The output will of this function is an array with the size of the array of dependent variable used during the generation of the model.
 */
void UPLSYPredictor(matrix *tscores, UPLSMODEL *m, size_t npc, array **py);

/*
 * Description: Calculation The RSquared for the objects in the model UPLSMODEL m. 
 * This function is used after UPLS model generation and the output is:
 *  ~ dvector r2x where the size correspond to the number of pricipal component calculated
 *  ~ array r2y and sdec where:
 *       - the order of the array correspond to the number of principal component
 *       - the row for each matrix in array correspond to an Y dependent value
 *       - the colunm for each matrix in array correspond to a layer of the Y dependen value, and so is equal to the Y->order if Y is the starting array for the model generation.
 * 
 */
void UPLSRSquared(array *X_, array *Y_, UPLSMODEL *m, size_t npc, dvector **r2x, array **r2y, array **sdec);

/*
 * Description: Calculation of the regression sum of squares and total sum of squares (proportional to the sample variance).
 * This function is used by UPLSCrossValidation and Leave One Out.
 * The output is:
 * 
 *  ~ dvector xss_err and xss_tot where the size correspond to the number of pricipal component calculated
 *  ~ array yss_err and  yss_tot where:
 *       - the order of the array correspond to the number of principal component
 *       - the row for each matrix in array correspond to an Y dependent value
 *       - the colunm for each matrix in array correspond to a layer of the Y dependen value, and so is equal to the Y->order if Y is the starting array for the model generation.
 * 
 */
void UPLSRSquared_SSErr_SStot(array *X_, array *Y_, UPLSMODEL *m, size_t npc, dvector** xss_err, dvector **xss_tot, array** yss_err, array** yss_tot, array** predicted_y);

/*
 * Description: Calculate the consistency of model to prevent Overfitting.
 * Thi funcion is used coupled with UPLSLOOValidation or UPLSCrossValidation in order to exstimate the overfitting of model.
 */
void UPLSYScrambling(array *X_, array *Y_, 
                        size_t xautoscaling, size_t yautoscaling,
                        size_t npc, size_t block,
                        size_t valtype, size_t rgcv_groups, size_t rgcv_iterations,
                        array **q2y, array **sdep, ssignal *s);

/*
 * Description: Random Cross Validation. Is used to extimate the model predictivity.
 * The output is: 
 *  ~ dvector r2x where the size correspond to the number of pricipal component calculated
 *  ~ array q2y and sdep where:
 *       - the order of the array correspond to the number of principal component
 *       - the row for each matrix in array correspond to an Y dependent value
 *       - the colunm for each matrix in array correspond to a layer of the Y dependen value, and so is equal to the Y->order if Y is the starting array for the model generation.
 * 
 */
void UPLSRandomGroupsCV(array *X_, array *Y_, 
                         size_t xautoscaling, size_t yautoscaling, 
                         size_t npc, size_t group, size_t iterations,
                         dvector **r2x, array **q2y, array **sdep, array **predicted_y, array **pred_residuals, ssignal *s);


/*
 * Description: Leave One Out Validation. Is used to extimate the model predictivity.
 * The output is: 
 *  ~ dvector r2x where the size correspond to the number of pricipal component calculated
 *  ~ array q2y and sdep where:
 *       - the order of the array correspond to the number of principal component
 *       - the row for each matrix in array correspond to an Y dependent value
 *       - the colunm for each matrix in array correspond to a layer of the Y dependen value, and so is equal to the Y->order if Y is the starting array for the model generation.
 * 
 */
void UPLSLOOCV(array *X_, array *Y_, 
                         size_t xautoscaling, size_t yautoscaling, 
                         size_t npc,
                         dvector **r2x, array **q2y, array **sdep, array **predicted_y, array **pred_residuals, ssignal *s);


/*This function from Q^2, or R^2 return the number of components to use for make predictions.
 * This cutoff components look at the best r^2 or q^2 and pay attention if r^2 or q^2 get low value 
 * and after high value..
 * 
 * rq2_randomness comes from the PLSModelConsistency that is the q2
 * N.B.: 0 means the first principal component, 1 is the second principal component. etc...
 */
size_t UPLSGetPCModelCutOff(array *rq2);
/* Variable Selection using the Particle Swarm Optimization algorithm
 * 
 * Input:
 *   mx and my matrix for make model; optionally px and py for testset validation.
 *   randomness = a number betwee 0 and 1 that represent the randomness to add on the first model generation.
 *   nrandomnessiter = for compute random models and evaluate also the model consistency...
 * 
 * Output:
 * - varselected vector: vector of varialble selected
 * - map = a matrix of 4 column with the best r2, q2, f-test and number of variables. 
 * - vardistribution = a vector of variables distribution accumulated in all model. Each row is a variable and the value of this row is the number of time that this variable was used.
 */
void UPSLPSOVariableSelection(array *ax, array *ay, array *px, array *py,
                       size_t xautoscaling, size_t yautoscaling, size_t npc, int validation_type, size_t ngroup, size_t niter, size_t population_size, double randomness,
                       uivector **varselected, matrix **map, uivector **vardistribution, ssignal *s);

/*
 * Variable selection using the Genetic Algorithm
 * 
 * Input 
 *   mx and my matrix for make model; optionally px and py for testset validation.
 *   
 * population_size = population of models to calculate; is a number > 2
 * fraction_of_population = fraction of population; is a number between 0 and 1 
 * mutation_rate = rate of mutation; is a number between 0 and 1 
 * crossovertype = 0 for single point crossover, 1 for two point crossover and 2 for uniform crossover
 * nswapping = used by the uniform crossover (crossovertype = 2) and is a number between 0 and 1. It means the % of variable to swapping
 * 
 * Output
 *   varselected vector: vector of varialble selected
 *   rq2 = coefficient of validation q2  or r2
 *  f_distribution = coefficient of distribution
 *  map = a matrix of 4 column with the best r2, q2, f-test and number of variables.
 */
void UPLSGAVariableSelection(array *ax, array *ay, array *px, array *py,
                       size_t xautoscaling, size_t yautoscaling, size_t npc, int validation_type, size_t ngroup, size_t niter, size_t population_size, double fraction_of_population, double mutation_rate, size_t crossovertype, double nswapping, double populationconvergence,
                      uivector **varselected, matrix **map, uivector **vardistribution, ssignal *s);

void PrintUPLSModel(UPLSMODEL *m);

#endif

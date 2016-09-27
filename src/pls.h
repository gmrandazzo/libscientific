/* pls.h
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

#ifndef PLS_H
#define PLS_H
#include "matrix.h"
#include "vector.h"
#include "scientificinfo.h"

#define PLSCONVERGENCE 1e-8

typedef struct{
  matrix *xscores;
  matrix *xloadings;
  matrix *xweights;
  matrix *yscores;
  matrix *yloadings;
  matrix *cweights;
  dvector *b;
  dvector *xvarexp;
  dvector *xcolaverage;
  dvector *xcolscaling;
  dvector *ycolaverage;
  dvector *ycolscaling;
  matrix *r2y_model; /* each column correspond to an y dependent variable and each row correspond to a principal component*/
  matrix *r2y_validation;
  matrix *recalc_residuals;
  matrix *q2y;
  matrix *sdep; /* Standard Deviation over Prediction */
  matrix *sdec; /* Standard Deviation over Recalculating */
  matrix *bias;
  matrix *recalculated_y;
  matrix *predicted_y;
  matrix *pred_residuals;
  matrix *r2q2scrambling;
  matrix *q2_sample_validation;
  matrix *sdep_sample_validation;
  matrix *q2_sample_validation_surface;
  matrix *sdep_sample_validation_surface;
} PLSMODEL;

/*
 * Description: Create a new PLSMODEL
 */
void NewPLSModel(PLSMODEL **m);

/*
 * Description: Delete a PLSMODEL
 */
void DelPLSModel(PLSMODEL **m);

/*
 * Description PLS calculation from P. Geladi algorithm
 */
void PLS(matrix *mx, matrix *my, size_t nlv, size_t xautoscaling, size_t yautoscaling, PLSMODEL *model, ssignal *s);

/*
 * Description: Calculate betas coefficients from a pls model at nlv latent variables
 */
void PLSBetasCoeff(PLSMODEL *model, size_t nlv, dvector **betas);

/*
 * Description: Project a matrix and predict the scores into the new space.
 * This function is used before predict the Y values
 */
void PLSScorePredictor(matrix *mx, PLSMODEL *model, size_t nlv, matrix **xscores);

/*
 * Description: Calculate the Y values.
 * N.B.: This function is dependent of PLSScorePredictor
 */
void PLSYPredictor(matrix *tscore, PLSMODEL *model, size_t nlv, matrix **y);

/*
 * Description: Calculate the R^2 of the model
 */
void PLSRSquared(matrix *mx, matrix *my, PLSMODEL *model, size_t nlv, matrix** r2y, matrix **sdec);

/*
 * Description: Calculate the SS error and the SS total to calculate
 * then the R^2 or other parameters
 */
void PLSRSquared_SSErr_SSTot(matrix *mx, matrix *my, PLSMODEL *model, size_t nlv, dvector** xss_err, dvector **xss_tot, matrix** yss_err, matrix** yss_tot, matrix **pred_y);

/*
 * Description: Calculate the PLS Very important variables
 */
void PLSVIP(PLSMODEL *model, matrix **vip);

/*
 * Description: Generate Random Models by putting randomly y
 * and check if models have correlations.
 */
void PLSYScrambling(matrix *mx, matrix *my,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t block,
                        size_t valtype, size_t rgcv_group, size_t rgcv_iterations,
                        matrix **r2q2scrambling, size_t nthreads, ssignal *s);

/*
 * Description: Calculate the PLS Bootstrap Random group cross validation.
 */
void PLSRandomGroupsCV(matrix *mx, matrix *my, size_t xautoscaling, /*Inputs*/
                       size_t yautoscaling, size_t nlv, size_t group, size_t iterations, /*Inputs*/
                      matrix **q2y, matrix **sdep, matrix **bias, matrix **predicted_y, /*Ouputs*/
                      matrix **pred_residuals, size_t nthreads, ssignal *s); /*Ouputs*/

/*
 * Description: Calculate the PLS leave one out cross validation.
 */
void PLSLOOCV(matrix *mx, matrix *my, size_t xautoscaling, size_t yautoscaling, size_t nlv,/*Inputs*/
                        matrix **q2y, matrix **sdep, matrix **bias, /*Ouputs*/
                        matrix **predicted_y, matrix **pred_residuals, /*Ouputs*/
                        size_t ntreads, ssignal *s); /*Inputs*/


/*
 * Description: Validata Sample Stability reducing or fixing the size to build a model
*/
void PLSStaticSampleValidator(matrix *mx, matrix *my, uivector *obj_class,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t sample_size, size_t niters,
                        size_t rgcv_group, size_t rgcv_iterations, size_t nthreads,
                        matrix **q2_distr, matrix **sdep_distr, uivector **bestid, ssignal *s);

/*
 * Description: Validate Model Stability using a dinamic incremental sample
 */
void PLSDynamicSampleValidator(matrix *mx, matrix *my,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t niters,
                        uivector *obj_class, size_t deltaobj, size_t maxobj,
                        size_t rgcv_group, size_t rgcv_iterations, size_t nthreads,
                        matrix **q2_surface, matrix **sdep_surface, uivector **bestid, ssignal *s);

/*
 * Description: Get the Cutoff based on the grow of r2 or q2
 */
int GetLVCCutoff(matrix* r2q2);

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
void PSLPSOVariableSelection(matrix *mx, matrix *my, matrix *px, matrix *py,
                       size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type, size_t ngroup, size_t niter,
                       size_t population_size, double randomness,
                       uivector **varselected, matrix **map, uivector **vardistribution, size_t nthreads, ssignal *s);

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
void PLSGAVariableSelection(matrix *mx, matrix *my, matrix *px, matrix *py,
                       size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type, size_t ngroup, size_t niter,
                       size_t population_size, double fraction_of_population, double mutation_rate, size_t crossovertype, double nswapping, double populationconvergence,
                       uivector **varselected, matrix **map, uivector **vardistribution, size_t nthreads, ssignal *s);

/* Variable selection using the Spearmann correlation coefficient
 * Input
 * mx and my matrix for make the model; optionally px and py for testset validation
 *
 * nlv = number of principal component
 * validation_type = validation type: 0 = LOO; 1 RGCV
 * ngroup = number of group for RGCV
 * niter = number of iterations for RGCV
 * nrandomnessiter = number of random iteration for Model consistency
 *
 * threshold = the step size to check if variable is in or out
 *
 * Output:
 *  varselected vector: vector of varialble selected
  *  map = a matrix of 4 column with the best r2, q2, treshold and number of variables.
 */
void PLSSpearmannVariableSelection(matrix *mx, matrix *my, matrix *px, matrix *py,
                                   size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type,
                                   size_t ngroup, size_t niter,
                                   double threshold,
                                   uivector **varselected, matrix **map, uivector **vardistribution, size_t nthreads, ssignal *s);

void PrintPLSModel(PLSMODEL *model);

#endif
